import json
import gzip
import math

import os
import subprocess
import glob
import time
import yaml
import copy

import multiprocessing as mp
from queue import Empty
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from openmm.app import (
    PDBFile,
    Modeller,
)

import openfe
from openfe import (
    SmallMoleculeComponent,
    ProteinMembraneComponent,
    ChemicalSystem,
    Transformation,
    AlchemicalNetwork,
)
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

from openff.units import unit

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from cinnabar import FEMap
from cinnabar import plotting as cinnabar_plotting

from openff.units import unit

import MDAnalysis as mda

from modules.utils.rbfe_setup import _prepare_openfe_receptor, _load_protein_membrane_component, _build_protocol_settings
from modules.utils.rbfe_worker import _openfe_worker

from pipeline.task_registry import register_task

from tqdm import tqdm

from backends.openmm_backend import (
    OpenMMBackend,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

from openfe.protocols.openmm_afe import AbsoluteBindingProtocol

from openfe.protocols.openmm_utils.omm_settings import (
    OpenFFPartialChargeSettings,
)

from openfe.protocols.openmm_utils.charge_generation import (
    bulk_assign_partial_charges,
)

def _force_neutralise_protonated_cations(mol):
    """
    Neutralise selected +1 heteroatoms only if they have a removable hydrogen.

    Intended for protonated amines, e.g.

        R3NH+ -> R3N

    This deliberately does not neutralise permanent cations such as
    quaternary ammonium, because those have no removable hydrogen.

    This implementation avoids stale RDKit atom indices by mutating only one
    atom at a time and then restarting the search.
    """

    rw = Chem.RWMol(mol)

    changed_any = False

    while True:
        changed_this_pass = False

        for atom in rw.GetAtoms():
            if atom.GetFormalCharge() != 1:
                continue

            if atom.GetSymbol() not in {
                "N",
                "O",
                "S",
                "P",
            }:
                continue

            atom_idx = atom.GetIdx()

            #
            # Case 1: explicit H atom neighbour.
            #
            h_neighbours = [
                nbr.GetIdx()
                for nbr in atom.GetNeighbors()
                if nbr.GetSymbol() == "H"
            ]

            if h_neighbours:
                h_idx = h_neighbours[0]

                #
                # Set charge on the heavy atom before removing H.
                # After RemoveAtom(), atom indices may shift, so do not
                # fetch atom_idx again in this pass.
                #
                atom.SetFormalCharge(0)
                atom.SetNoImplicit(False)

                rw.RemoveAtom(h_idx)

                changed_any = True
                changed_this_pass = True

                break

            #
            # Case 2: explicit H count rather than explicit H atom.
            #
            if atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(
                    atom.GetNumExplicitHs() - 1
                )
                atom.SetFormalCharge(0)
                atom.SetNoImplicit(False)

                changed_any = True
                changed_this_pass = True

                break

            #
            # Case 3: implicit H.
            #
            if atom.GetNumImplicitHs() > 0:
                atom.SetFormalCharge(0)
                atom.SetNoImplicit(False)

                changed_any = True
                changed_this_pass = True

                break

        if not changed_this_pass:
            break

    neutral = rw.GetMol()

    if changed_any:
        neutral.UpdatePropertyCache(
            strict=False
        )

        Chem.SanitizeMol(
            neutral
        )

    return neutral, changed_any

def _charged_atom_report(mol):
    """
    Return per-atom formal charge diagnostics for non-neutral atoms.
    """
    rows = []

    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()

        if charge == 0:
            continue

        rows.append(
            {
                "atom_idx": atom.GetIdx(),
                "symbol": atom.GetSymbol(),
                "formal_charge": charge,
                "total_num_hs": atom.GetTotalNumHs(),
                "explicit_hs": atom.GetNumExplicitHs(),
                "implicit_hs": atom.GetNumImplicitHs(),
                "degree": atom.GetDegree(),
                "is_aromatic": atom.GetIsAromatic(),
                "smarts": atom.GetSmarts(),
            }
        )

    return rows

def _formal_charge(mol):
    """
    Return total formal charge of an RDKit molecule.
    """
    return int(
        sum(
            atom.GetFormalCharge()
            for atom in mol.GetAtoms()
        )
    )


def _largest_fragment(mol):
    """
    Keep the largest covalent fragment.

    Useful when SDF records contain salts/counterions.
    """
    chooser = rdMolStandardize.LargestFragmentChooser()
    return chooser.choose(mol)


def _neutralise_rdkit_mol(
    mol,
    force_protonated_cations=False,
):
    """
    Neutralise common formal charges using RDKit's Uncharger.

    If requested, additionally neutralise protonated +1 heteroatoms that have
    a removable hydrogen.

    This does not neutralise permanent cations.
    """
    uncharger = rdMolStandardize.Uncharger()

    neutral = uncharger.uncharge(
        mol
    )

    Chem.SanitizeMol(
        neutral
    )

    if (
        _formal_charge(neutral) != 0
        and force_protonated_cations
    ):
        neutral, changed = _force_neutralise_protonated_cations(
            neutral
        )

        if changed:
            Chem.SanitizeMol(
                neutral
            )

    return neutral

def _small_molecule_component_from_rdkit(mol):
    """
    Version-tolerant construction of an OpenFE SmallMoleculeComponent.
    """
    if hasattr(SmallMoleculeComponent, "from_rdkit"):
        return SmallMoleculeComponent.from_rdkit(mol)

    return SmallMoleculeComponent(mol)


def _load_openfe_ligands_from_sdf(
    ligand_sdf,
    ligand_name_property="_Name",
):
    """
    Load RDKit molecules from SDF and return OpenFE SmallMoleculeComponents.

    Ensures every ligand has a stable name because OpenFE uses component names
    in transformation and result labels.
    """
    ligands = []

    supplier = Chem.SDMolSupplier(
        ligand_sdf,
        removeHs=False,
    )

    for idx, mol in enumerate(supplier):
        if mol is None:
            logger.warning(
                f"Skipping unreadable molecule {idx} in {ligand_sdf}"
            )
            continue

        name = None

        if ligand_name_property and mol.HasProp(ligand_name_property):
            name = mol.GetProp(ligand_name_property).strip()

        if not name and mol.HasProp("_Name"):
            name = mol.GetProp("_Name").strip()

        if not name:
            name = f"ligand_{idx:04d}"

        mol.SetProp("_Name", name)

        formal_charge = _formal_charge(
            mol
        )

        if formal_charge != 0:
            raise ValueError(
                f"Ligand {name} has non-zero formal charge "
                f"({formal_charge}). OpenFE ABFE currently rejects charged "
                f"alchemical molecules in the solvation free energy leg. "
                f"Run openfe_neutralise_ligand_sdf first, or use a neutral "
                f"input species."
            )

        ligands.append(
            _small_molecule_component_from_rdkit(mol)
        )

    if not ligands:
        raise ValueError(
            f"No valid ligands loaded from {ligand_sdf}"
        )

    return ligands


def _maybe_charge_ligands_abfe(
    ligands,
    charge_cfg,
):
    """
    Optionally assign a single reproducible set of partial charges per ligand.

    OpenFE's ABFE tutorial recommends assigning a single set of charges for
    reproducibility between repeats.
    """
    if not charge_cfg:
        return ligands

    charge_settings = OpenFFPartialChargeSettings(
        partial_charge_method=charge_cfg.get(
            "partial_charge_method",
            "am1bcc",
        ),
        off_toolkit_backend=charge_cfg.get(
            "off_toolkit_backend",
            "ambertools",
        ),
    )

    logger.info(
        "Assigning ligand partial charges: "
        f"method={charge_settings.partial_charge_method}, "
        f"backend={charge_settings.off_toolkit_backend}"
    )

    return bulk_assign_partial_charges(
        molecules=ligands,
        overwrite=charge_cfg.get("overwrite", False),
        method=charge_settings.partial_charge_method,
        toolkit_backend=charge_settings.off_toolkit_backend,
        generate_n_conformers=charge_settings.number_of_conformers,
        nagl_model=charge_settings.nagl_model,
        processors=charge_cfg.get("processors", 1),
    )


def _build_abfe_protocol_settings(protocol_cfg):
    """
    Build AbsoluteBindingProtocol settings from YAML.

    This deliberately keeps OpenFE's default ABFE lambda schedules unless you
    explicitly extend this function later. ABFE lambda schedules are not the
    same thing as your RBFE 11-window setup.
    """
    settings = AbsoluteBindingProtocol.default_settings()

    repeats = protocol_cfg.get("repeats", None)
    if repeats is not None:
        settings.protocol_repeats = int(repeats)

    platform = protocol_cfg.get("platform", None)
    if platform is not None:
        settings.engine_settings.compute_platform = str(platform)

    gpu_id = protocol_cfg.get("gpu_id", None)
    if gpu_id is not None:
        # OpenMM itself will see CUDA_VISIBLE_DEVICES from the worker.
        # This is logged here for traceability.
        logger.info(
            f"ABFE protocol requested gpu_id={gpu_id}"
        )

    #
    # Bound-state restraint settings.
    #
    if "host_min_distance_nm" in protocol_cfg:
        settings.restraint_settings.host_min_distance = (
            float(protocol_cfg["host_min_distance_nm"])
            * unit.nanometer
        )

    if "host_max_distance_nm" in protocol_cfg:
        settings.restraint_settings.host_max_distance = (
            float(protocol_cfg["host_max_distance_nm"])
            * unit.nanometer
        )

    #
    # Generic simulation lengths.
    #
    production_length_ns = protocol_cfg.get(
        "production_length_ns",
        None,
    )
    equilibration_length_ns = protocol_cfg.get(
        "equilibration_length_ns",
        None,
    )
    equilibration_length_nvt_ns = protocol_cfg.get(
        "equilibration_length_nvt_ns",
        None,
    )

    #
    # Complex leg.
    #
    complex_prod = protocol_cfg.get(
        "complex_production_length_ns",
        production_length_ns,
    )
    if complex_prod is not None:
        settings.complex_simulation_settings.production_length = (
            float(complex_prod) * unit.nanosecond
        )

    complex_equil = protocol_cfg.get(
        "complex_equilibration_length_ns",
        equilibration_length_ns,
    )
    if complex_equil is not None:
        settings.complex_simulation_settings.equilibration_length = (
            float(complex_equil) * unit.nanosecond
        )

    #
    # Solvent leg.
    #
    solvent_prod = protocol_cfg.get(
        "solvent_production_length_ns",
        production_length_ns,
    )
    if solvent_prod is not None:
        settings.solvent_simulation_settings.production_length = (
            float(solvent_prod) * unit.nanosecond
        )

    solvent_equil = protocol_cfg.get(
        "solvent_equilibration_length_ns",
        equilibration_length_ns,
    )
    if solvent_equil is not None:
        settings.solvent_simulation_settings.equilibration_length = (
            float(solvent_equil) * unit.nanosecond
        )

    #
    # Equilibration stage for complex.
    #
    if equilibration_length_ns is not None:
        settings.complex_equil_simulation_settings.equilibration_length = (
            float(equilibration_length_ns) * unit.nanosecond
        )

    if equilibration_length_nvt_ns is not None:
        settings.complex_equil_simulation_settings.equilibration_length_nvt = (
            float(equilibration_length_nvt_ns) * unit.nanosecond
        )

    #
    # Output/checkpoint settings.
    #
    checkpoint_interval_ns = protocol_cfg.get(
        "checkpoint_interval_ns",
        None,
    )
    if checkpoint_interval_ns is not None:
        settings.complex_output_settings.checkpoint_interval = (
            float(checkpoint_interval_ns) * unit.nanosecond
        )
        settings.solvent_output_settings.checkpoint_interval = (
            float(checkpoint_interval_ns) * unit.nanosecond
        )
        settings.complex_equil_output_settings.checkpoint_interval = (
            float(checkpoint_interval_ns) * unit.nanosecond
        )

    positions_write_frequency_ps = protocol_cfg.get(
        "positions_write_frequency_ps",
        None,
    )
    if positions_write_frequency_ps is not None:
        settings.complex_output_settings.positions_write_frequency = (
            float(positions_write_frequency_ps) * unit.picosecond
        )
        settings.solvent_output_settings.positions_write_frequency = (
            float(positions_write_frequency_ps) * unit.picosecond
        )

    if "lambda_windows" in protocol_cfg:
        logger.warning(
            "Ignoring protocol.lambda_windows for ABFE. "
            "OpenFE ABFE uses separate complex/solvent lambda schedules. "
            "Keep defaults unless you explicitly implement custom schedules."
        )

    logger.info(
        "ABFE settings: "
        f"repeats={settings.protocol_repeats}, "
        f"platform={settings.engine_settings.compute_platform}, "
        f"complex_prod={settings.complex_simulation_settings.production_length}, "
        f"solvent_prod={settings.solvent_simulation_settings.production_length}"
    )

    return settings


def _open_json_maybe_gz(path):
    """
    Load JSON from .json or .json.gz.
    """
    path = str(path)

    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return json.load(f)

    with open(path, "r") as f:
        return json.load(f)


def _quantity_to_float_kcal_per_mol(value):
    """
    Convert common OpenFE/gufe JSON quantity representations to kcal/mol float.

    Handles:
      - plain float
      - string-ish float
      - keyed dicts with val/unit
      - nested OpenFF quantity dicts seen in gufe JSON outputs
    """
    if value is None:
        return None

    if isinstance(value, (int, float)):
        return float(value)

    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return value

    if isinstance(value, dict):
        if "val" in value:
            return float(value["val"])

        if "magnitude" in value:
            return float(value["magnitude"])

        if "value" in value:
            return _quantity_to_float_kcal_per_mol(value["value"])

        # gufe/openff keyed quantities can be nested.
        for key in [
            "estimate",
            "uncertainty",
            "m",
        ]:
            if key in value:
                return _quantity_to_float_kcal_per_mol(value[key])

    return value