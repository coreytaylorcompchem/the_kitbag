import os
import json
import sys
import platform
import re

from pathlib import Path
import json
import numpy as np

import openfe
from openfe import (
    SmallMoleculeComponent,
    ProteinMembraneComponent,
    ChemicalSystem,
    Transformation,
    AlchemicalNetwork,
)

from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

from rdkit import Chem
from rdkit.Chem import AllChem

import MDAnalysis as mda

from pipeline.task_registry import register_task

from tqdm import tqdm

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def _load_box_vectors_from_state_json(json_file):
    with open(json_file) as f:
        data = json.load(f)

    return np.asarray(data["box_vectors"])

@register_task(
    "openfe_prepare_receptor",
    category="Free Energy",
    description="Prepare membrane receptor for OpenFE."
)
def openfe_prepare_receptor(self):

    cfg = self.config["openfe_prepare_receptor"]

    pdb_file = cfg["equilibrated_pdb"]
    ligand_resname = cfg.get("ligand_resname", "LIG")

    output_dir = cfg.get("output_dir", "openfe")
    os.makedirs(output_dir, exist_ok=True)

    u = mda.Universe(pdb_file)

    receptor = u.select_atoms(
        f"not resname {ligand_resname}"
    )

    receptor_pdb = os.path.join(
        output_dir,
        "receptor_membrane.pdb"
    )

    receptor.write(receptor_pdb)

    logger.info(
        f"Wrote receptor file: {receptor_pdb}"
    )

    with open(receptor_pdb) as f:
        first_lines = [next(f) for _ in range(20)]

    logger.info(
        f"CRYST1 present: "
        f"{any('CRYST1' in x for x in first_lines)}"
    )

    logger.info(
        f"Saved receptor membrane system "
        f"to {receptor_pdb}"
    )

    self.openfe_receptor_pdb = receptor_pdb

@register_task(
    "openfe_create_network",
    category="Free Energy",
    description="Create OpenFE ligand network for membrane RBFE."
)
def openfe_create_network(self):

    cfg = self.config["openfe_create_network"]

    receptor_pdb = cfg["receptor_pdb"]
    ligand_sdf = cfg["ligand_sdf"]

    output_dir = cfg.get(
        "output_dir",
        "openfe"
    )

    os.makedirs(output_dir, exist_ok=True)

    receptor = (
        ProteinMembraneComponent
        .from_pdb_file(
            receptor_pdb
        )
    )

    supplier = Chem.SDMolSupplier(
        ligand_sdf,
        removeHs=False
    )

    ligands = []

    for i, mol in enumerate(supplier):

        if mol is None:
            logger.warning(
                f"Failed to read molecule {i}"
            )
            continue

        ligands.append(
            SmallMoleculeComponent(
                mol
            )
        )

    ligands = []

    for mol in Chem.SDMolSupplier(
        ligand_sdf,
        removeHs=False
    ):

        if mol is None:
            continue

        ligands.append(
            SmallMoleculeComponent(
                mol
            )
        )

    mapper = openfe.setup.KartografAtomMapper()

    scorer = (
        openfe.lomap_scorers
        .default_lomap_score
    )

    network = (
        openfe.ligand_network_planning
        .generate_lomap_network(
            ligands=ligands,
            mappers=[mapper],
            scorer=scorer,
        )
    )

    network_json = os.path.join(
        output_dir,
        "ligand_network.json"
    )

    with open(network_json, "w") as f:
        f.write(network.to_json())

    logger.info(
        f"Wrote OpenFE network: "
        f"{network_json}"
    )

    self.openfe_network = network
    self.openfe_receptor = receptor

    logger.info(
        f"Generated network with "
        f"{len(network.nodes)} ligands and "
        f"{len(network.edges)} edges"
    )

@register_task(
    "openfe_create_alchemical_network",
    category="Free Energy",
    description="Create OpenFE RBFE alchemical network."
)
def openfe_create_alchemical_network(self):

    cfg = self.config["openfe_create_alchemical_network"]

    receptor_pdb = cfg["receptor_pdb"]

    output_dir = cfg.get(
        "output_dir",
        "openfe"
    )

    os.makedirs(output_dir, exist_ok=True)

    ligand_network = self.openfe_network

    protein = openfe.ProteinMembraneComponent.from_pdb_file(
        receptor_pdb
    )

    solvent = openfe.SolventComponent()

    transformations = []

    for mapping in ligand_network.edges:

        for leg in ["solvent", "complex"]:

            if leg == "solvent":

                sysA_dict = {
                    "ligand": mapping.componentA,
                    "solvent": solvent,
                }

                sysB_dict = {
                    "ligand": mapping.componentB,
                    "solvent": solvent,
                }

            else:

                sysA_dict = {
                    "ligand": mapping.componentA,
                    "solvent": solvent,
                    "protein": protein,
                }

                sysB_dict = {
                    "ligand": mapping.componentB,
                    "solvent": solvent,
                    "protein": protein,
                }

            system_a = ChemicalSystem(
                sysA_dict,
                name=f"{mapping.componentA.name}_{leg}"
            )

            system_b = ChemicalSystem(
                sysB_dict,
                name=f"{mapping.componentB.name}_{leg}"
            )

            settings = (
                RelativeHybridTopologyProtocol
                .default_settings()
            )

            settings.protocol_repeats = 1

            protocol = RelativeHybridTopologyProtocol(
                settings=settings
            )

            transformation = Transformation(
                stateA=system_a,
                stateB=system_b,
                mapping=mapping,
                protocol=protocol,
                name=(
                    f"rbfe_"
                    f"{system_a.name}_"
                    f"{system_b.name}"
                )
            )

            transformations.append(
                transformation
            )

    alchemical_network = AlchemicalNetwork(
        transformations
    )

    self.alchemical_network = (
        alchemical_network
    )

    #
    # Write transformations
    #

    transforms_dir = os.path.join(
        output_dir,
        "transformations"
    )

    os.makedirs(
        transforms_dir,
        exist_ok=True
    )

    for i, transformation in enumerate(
        alchemical_network.edges
    ):

        transformation.dump(
            os.path.join(
                transforms_dir,
                f"{i:03d}_{transformation.name}.json"
            )
        )

    logger.info(
        f"Wrote "
        f"{len(alchemical_network.edges)} "
        f"transformations"
    )