import math
import os
import glob

from pathlib import Path

import pandas as pd

import openfe
from openfe import (
    ChemicalSystem,
    Transformation,
    AlchemicalNetwork,
)


from openff.units import unit

from rdkit import Chem
from rdkit.Chem import AllChem

from cinnabar import FEMap
from cinnabar import plotting as cinnabar_plotting

from openff.units import unit

import MDAnalysis as mda

from modules.utils.rbfe_setup import _load_protein_membrane_component

from pipeline.task_registry import register_task

from tqdm import tqdm

from openfe.protocols.openmm_afe import AbsoluteBindingProtocol

from openfe.protocols.openmm_utils.omm_settings import (
    OpenFFPartialChargeSettings,
)

from openfe.protocols.openmm_utils.charge_generation import (
    bulk_assign_partial_charges,
)

from modules.utils.abfe_setup import(
    _load_openfe_ligands_from_sdf,
    _maybe_charge_ligands_abfe,
    _build_abfe_protocol_settings,
    _open_json_maybe_gz,
    _quantity_to_float_kcal_per_mol,
    _formal_charge,
    _largest_fragment,
    _neutralise_rdkit_mol,
    _charged_atom_report
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task(
    "openfe_neutralise_ligand_sdf",
    category="Free Energy",
    description="Neutralise ligand SDF formal charges for OpenFE ABFE."
)
def openfe_neutralise_ligand_sdf(self):
    """
    Neutralise formal molecular charges in an SDF before ABFE.

    This is intended for OpenFE ABFE, where charged alchemical molecules
    currently fail during the solvation free energy leg.

    Important:
      - This changes the molecular species.
      - It should be used only if the neutral form is chemically meaningful.
      - Partial charge assignment is still done afterwards by OpenFE/OpenFF.
    """

    cfg = self.config["openfe_neutralise_ligand_sdf"]

    input_sdf = cfg["input_sdf"]

    output_sdf = cfg.get(
        "output_sdf",
        os.path.splitext(input_sdf)[0] + "_neutral.sdf",
    )

    fail_if_charged = cfg.get(
        "fail_if_charged",
        True,
    )

    remove_salts = cfg.get(
        "remove_salts",
        True,
    )

    keep_largest_fragment = cfg.get(
        "keep_largest_fragment",
        True,
    )

    output_dir = os.path.dirname(
        output_sdf
    )

    if output_dir:
        os.makedirs(
            output_dir,
            exist_ok=True,
        )

    supplier = Chem.SDMolSupplier(
        input_sdf,
        removeHs=False,
    )

    writer = Chem.SDWriter(
        output_sdf,
    )

    records = []
    n_written = 0
    n_failed = 0

    for idx, mol in enumerate(
        supplier
    ):
        if mol is None:
            logger.warning(
                f"[{idx}] Could not parse molecule in {input_sdf}"
            )
            n_failed += 1
            continue

        name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else f"mol_{idx:04d}"

        try:
            Chem.SanitizeMol(
                mol
            )

            charge_before = _formal_charge(
                mol
            )

            working = Chem.Mol(
                mol
            )

            if remove_salts and keep_largest_fragment:
                working = _largest_fragment(
                    working
                )

                Chem.SanitizeMol(
                    working
                )

            charge_after_fragment = _formal_charge(
                working
            )

            force_protonated_cations = cfg.get(
                "force_protonated_cations",
                False,
                )

            neutral = _neutralise_rdkit_mol(
                working,
                force_protonated_cations=force_protonated_cations,
            )

            charge_after = _formal_charge(
                neutral
            )

            neutral.SetProp(
                "_Name",
                name,
            )

            #
            # Preserve useful input properties where possible.
            #
            for prop in mol.GetPropNames():
                if prop == "_Name":
                    continue

                try:
                    neutral.SetProp(
                        prop,
                        mol.GetProp(prop),
                    )
                except Exception:
                    pass

            neutral.SetIntProp(
                "formal_charge_before",
                int(charge_before),
            )

            neutral.SetIntProp(
                "formal_charge_after_fragment",
                int(charge_after_fragment),
            )

            neutral.SetIntProp(
                "formal_charge_after_neutralisation",
                int(charge_after),
            )

            records.append(
                {
                    "index": idx,
                    "name": name,
                    "formal_charge_before": charge_before,
                    "formal_charge_after_fragment": charge_after_fragment,
                    "formal_charge_after_neutralisation": charge_after,
                    "status": "ok" if charge_after == 0 else "still_charged",
                }
            )

            if charge_after != 0:
                charged_atoms = _charged_atom_report(
                    neutral
                )

                for atom_info in charged_atoms:
                    logger.warning(
                        f"[{idx}] {name}: charged atom after neutralisation: "
                        f"idx={atom_info['atom_idx']}, "
                        f"symbol={atom_info['symbol']}, "
                        f"charge={atom_info['formal_charge']}, "
                        f"totalHs={atom_info['total_num_hs']}, "
                        f"explicitHs={atom_info['explicit_hs']}, "
                        f"implicitHs={atom_info['implicit_hs']}, "
                        f"degree={atom_info['degree']}, "
                        f"aromatic={atom_info['is_aromatic']}"
                    )

                msg = (
                    f"[{idx}] {name}: still charged after neutralisation. "
                    f"before={charge_before}, "
                    f"after_fragment={charge_after_fragment}, "
                    f"after_neutralisation={charge_after}. "
                    f"If the charged atom has no removable H, this is likely a "
                    f"permanent cation and should not be force-neutralised."
                )

                if fail_if_charged:
                    records.append(
                        {
                            "index": idx,
                            "name": name,
                            "formal_charge_before": charge_before,
                            "formal_charge_after_fragment": charge_after_fragment,
                            "formal_charge_after_neutralisation": charge_after,
                            "status": "failed_still_charged",
                        }
                    )

                    pd.DataFrame(
                        records
                    ).to_csv(
                        os.path.splitext(output_sdf)[0]
                        + "_neutralisation_report.csv",
                        index=False,
                    )

                    writer.close()

                    raise ValueError(
                        msg
                    )

                logger.warning(
                    msg
                )

            writer.write(
                neutral
            )

            n_written += 1

            logger.info(
                f"[{idx}] {name}: formal charge "
                f"{charge_before} -> {charge_after}"
            )

        except Exception as exc:
            n_failed += 1

            if fail_if_charged:
                writer.close()
                raise

            logger.warning(
                f"[{idx}] Failed to neutralise {name}: {exc}"
            )

    writer.close()

    if n_written == 0:
        raise ValueError(
            f"No molecules written to {output_sdf}"
        )

    report_file = os.path.splitext(
        output_sdf
    )[0] + "_neutralisation_report.csv"

    pd.DataFrame(
        records
    ).to_csv(
        report_file,
        index=False,
    )

    self.openfe_ligand_sdf = output_sdf

    logger.info(
        f"Wrote neutralised SDF: {output_sdf}"
    )

    logger.info(
        f"Wrote neutralisation report: {report_file}"
    )

    logger.info(
        f"Neutralisation summary: written={n_written}, failed={n_failed}"
    )

@register_task(
    "openfe_create_abfe_transformations",
    category="Free Energy",
    description="Create OpenFE ABFE transformations."
)
def openfe_create_abfe_transformations(self):
    """
    Create one OpenFE AbsoluteBindingProtocol transformation per ligand.

    ABFE differs from RBFE:
      - no ligand network
      - no atom mapping
      - one transformation per ligand
      - stateA contains ligand + receptor/membrane system
      - stateB contains receptor/membrane system without ligand
    """

    cfg = self.config["openfe_create_abfe_transformations"]
    protocol_cfg = cfg.get(
        "protocol",
        {},
    )

    receptor_pdb = getattr(
        self,
        "openfe_receptor_pdb",
        cfg["receptor_pdb"],
    )

    ligand_sdf = getattr(
        self,
        "openfe_ligand_sdf",
        cfg["ligand_sdf"],
    )

    output_dir = cfg.get(
        "output_dir",
        "openfe_abfe",
    )

    os.makedirs(
        output_dir,
        exist_ok=True,
    )

    transforms_dir = os.path.join(
        output_dir,
        "transformations",
    )

    os.makedirs(
        transforms_dir,
        exist_ok=True,
    )

    logger.info(
        f"Loading ABFE receptor/membrane component from {receptor_pdb}"
    )

    protein = _load_protein_membrane_component(
        receptor_pdb,
    )

    ligands = _load_openfe_ligands_from_sdf(
        ligand_sdf=ligand_sdf,
        ligand_name_property=cfg.get(
            "ligand_name_property",
            "_Name",
        ),
    )

    logger.info(
        f"Loaded {len(ligands)} ligand(s) from {ligand_sdf}"
    )

    if cfg.get("charge_ligands", False):
        ligands = _maybe_charge_ligands_abfe(
            ligands,
            cfg.get("charge", {}),
        )

    settings = _build_abfe_protocol_settings(
        protocol_cfg,
    )

    protocol = AbsoluteBindingProtocol(
        settings=settings,
    )

    transformations = []

    for idx, ligand in enumerate(
        ligands,
        start=1,
    ):
        ligand_name = ligand.name or f"ligand_{idx:04d}"

        logger.info(
            f"[{idx}/{len(ligands)}] Creating ABFE transformation "
            f"for {ligand_name}"
        )

        #
        # For a ProteinMembraneComponent/SolvatedPDBComponent, the receptor
        # component already contains explicit membrane/water/box, so we
        # mirror the setup in rbfe so avoid duplication

        solvent = openfe.SolventComponent()

        system_a = ChemicalSystem(
            {
                "ligand": ligand,
                "protein": protein,
                "solvent": solvent,
            },
            name=f"{ligand_name}_complex",
        )

        system_b = ChemicalSystem(
            {
                "protein": protein,
                "solvent": solvent,
            },
            name=f"{ligand_name}_apo",
        )

        transformation = Transformation(
            stateA=system_a,
            stateB=system_b,
            mapping=None,
            protocol=protocol,
            name=f"abfe_{ligand_name}",
        )

        transformations.append(
            transformation,
        )

    alchemical_network = AlchemicalNetwork(
        transformations,
    )

    self.alchemical_network = alchemical_network

    written = []

    for i, transformation in enumerate(
        alchemical_network.edges,
    ):
        safe_name = (
            transformation.name
            .replace("/", "_")
            .replace(" ", "_")
            .replace(":", "_")
        )

        out_json = os.path.join(
            transforms_dir,
            f"{i:03d}_{safe_name}.json",
        )

        transformation.dump(
            out_json,
        )

        written.append(out_json)

    logger.info(
        f"Wrote {len(written)} ABFE transformation JSON files "
        f"to {transforms_dir}"
    )

    self.openfe_abfe_transformations_dir = transforms_dir

@register_task(
    "openfe_gather_abfe_results",
    category="Free Energy",
    description="Gather OpenFE ABFE result JSON files."
)
def openfe_gather_abfe_results(self):
    """
    Gather ABFE quickrun outputs.

    OpenFE quickrun writes a JSON containing estimate and uncertainty.
    This task extracts those fields into a CSV.
    """

    cfg = self.config["openfe_gather_abfe_results"]

    results_dir = getattr(
        self,
        "openfe_results_dir",
        cfg["results_dir"],
    )

    output_file = cfg.get(
        "output_file",
        os.path.join(
            results_dir,
            "abfe_results.csv",
        ),
    )

    json_files = sorted(
        glob.glob(
            os.path.join(
                results_dir,
                "*.json",
            )
        )
    )

    json_files += sorted(
        glob.glob(
            os.path.join(
                results_dir,
                "*.json.gz",
            )
        )
    )

    if not json_files:
        raise ValueError(
            f"No ABFE result JSON files found in {results_dir}"
        )

    records = []

    for result_json in json_files:
        stem = Path(result_json).stem

        if stem.endswith(".json"):
            stem = Path(stem).stem

        try:
            result = _open_json_maybe_gz(
                result_json,
            )
        except Exception as exc:
            logger.warning(
                f"Could not read ABFE result JSON {result_json}: {exc}"
            )
            continue

        estimate = _quantity_to_float_kcal_per_mol(
            result.get("estimate")
        )

        uncertainty = _quantity_to_float_kcal_per_mol(
            result.get("uncertainty")
        )

        ligand = stem

        if ligand.startswith("000_abfe_"):
            ligand = ligand.replace("000_abfe_", "", 1)
        elif "_abfe_" in ligand:
            ligand = ligand.split("_abfe_", 1)[1]
        elif ligand.startswith("abfe_"):
            ligand = ligand.replace("abfe_", "", 1)

        records.append(
            {
                "ligand": ligand,
                "result_file": result_json,
                "dG_kcal_mol": estimate,
                "uncertainty_kcal_mol": uncertainty,
            }
        )

    if not records:
        raise ValueError(
            f"No readable ABFE results found in {results_dir}"
        )

    df = pd.DataFrame(records)

    #
    # Optional approximate conversion to pKi-like value.
    # pK = -dG / (RT ln 10)
    # At 298.15 K, RT ln 10 is about 1.364 kcal/mol.
    #
    temperature_k = cfg.get(
        "temperature_k",
        298.15,
    )

    r_kcal_mol_k = 0.00198720425864083
    rtln10 = r_kcal_mol_k * float(temperature_k) * math.log(10.0)

    df["approx_pK_from_dG"] = (
        -df["dG_kcal_mol"] / rtln10
    )

    df = df.sort_values(
        "dG_kcal_mol",
        ascending=True,
    )

    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(
            output_dir,
            exist_ok=True,
        )

    df.to_csv(
        output_file,
        index=False,
    )

    logger.info(
        f"Wrote ABFE results: {output_file}"
    )

    logger.info(
        "ABFE result summary:\n"
        + df[
            [
                "ligand",
                "dG_kcal_mol",
                "uncertainty_kcal_mol",
                "approx_pK_from_dG",
            ]
        ].to_string(index=False)
    )

    self.openfe_abfe_results_table = output_file