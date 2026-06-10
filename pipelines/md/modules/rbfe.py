import os
import json
import sys
import platform
import re
import subprocess
import glob
import time
import yaml

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

from openff.units import unit

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

def _get_openfe_progress(workdir):

    matches = glob.glob(
        os.path.join(
            workdir,
            "**",
            "simulation_real_time_analysis.yaml"
        ),
        recursive=True
    )

    if not matches:
        return None

    try:

        with open(matches[0]) as f:
            data = yaml.safe_load(f)

        if not data:
            return None

        latest = data[-1]

        return latest.get(
            "percent_complete"
        )

    except Exception:

        return None

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

    # mapper = openfe.setup.KartografAtomMapper()

    # scorer = (
    #     openfe.lomap_scorers
    #     .default_lomap_score
    # )

    # network = (
    #     openfe.ligand_network_planning
    #     .generate_lomap_network(
    #         ligands=ligands,
    #         mappers=[mapper],
    #         scorer=scorer,
    #     )
    # )
    mapper = openfe.setup.KartografAtomMapper()

    scorer = (
        openfe.lomap_scorers
        .default_lomap_score
    )

    network_cfg = cfg.get(
        "network",
        {}
    )

    network_method = network_cfg.get(
        "method",
        "lomap"
    )

    if network_method == "lomap":

        network = (
            openfe.ligand_network_planning
            .generate_lomap_network(
                ligands=ligands,
                mappers=[mapper],
                scorer=scorer,
            )
        )

    elif network_method == "mst":

        network = (
            openfe.ligand_network_planning
            .generate_minimal_spanning_network(
                ligands=ligands,
                mappers=[mapper],
                scorer=scorer,
            )
        )

    elif network_method == "star":

        center_idx = network_cfg.get(
            "central_ligand",
            0
        )

        network = (
            openfe.ligand_network_planning
            .generate_radial_network(
                ligands=ligands,
                central_ligand=ligands[center_idx],
                mappers=[mapper],
                scorer=scorer,
            )
        )

    else:

        raise ValueError(
            f"Unknown network method "
            f"{network_method}"
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
    protocol_cfg = cfg.get(
        "protocol",
        {}
    )

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
                    # "solvent": solvent,
                    "protein": protein,
                }

                sysB_dict = {
                    "ligand": mapping.componentB,
                    # "solvent": solvent,
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

            settings.protocol_repeats = (
                protocol_cfg.get(
                    "repeats",
                    1
                )
            )

            if "production_length_ns" in protocol_cfg:

                settings.simulation_settings.production_length = (
                    protocol_cfg["production_length_ns"]
                    * unit.nanosecond
                )

            if "equilibration_length_ns" in protocol_cfg:

                settings.simulation_settings.equilibration_length = (
                    protocol_cfg["equilibration_length_ns"]
                    * unit.nanosecond
                )
            
            if "lambda_windows" in protocol_cfg:

                n = protocol_cfg["lambda_windows"]

                settings.lambda_settings.lambda_windows = n

            # setup will blindly apply barostat settings if we don't set it explicitly for the complex legs.

            if leg == "complex":

                settings.integrator_settings.barostat = (
                    "MonteCarloMembraneBarostat"
                )

            logger.debug(
                f"Production length default: "
                f"{settings.simulation_settings.production_length}"
            )

            logger.debug(
                f"Equilibration length default: "
                f"{settings.simulation_settings.equilibration_length}"
            )

            protocol = RelativeHybridTopologyProtocol(
                settings=settings
            )

            logger.debug(
                f"{leg}: "
                f"{list(sysA_dict.keys())}"
            )

            logger.debug(
                f"{leg}: "
                f"repeats="
                f"{settings.protocol_repeats}"
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

@register_task(
    "openfe_run_network",
    category="Free Energy",
    description="Run all OpenFE transformations sequentially."
)
def openfe_run_network(self):

    cfg = self.config["openfe_run_network"]

    transforms_dir = cfg["transformations_dir"]

    output_dir = cfg.get(
        "output_dir",
        "openfe_results"
    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    json_files = sorted(
        glob.glob(
            os.path.join(
                transforms_dir,
                "*.json"
            )
        )
    )

    total = len(json_files)

    logger.info(
        f"Found {total} transformations"
    )

    failed = []

    for i, transform_json in enumerate(
        json_files,
        start=1
    ):

        stem = Path(transform_json).stem

        result_json = os.path.join(
            output_dir,
            f"{stem}.json"
        )

        working_dir = os.path.join(
            output_dir,
            stem
        )

        #
        # Skip completed calculations
        #

        if os.path.exists(result_json):

            logger.info(
                f"[{i}/{total}] "
                f"Skipping completed {stem}"
            )

            continue

        logger.info(
            f"[{i}/{total}] "
            f"Running {stem}"
        )

        cmd = [
            "openfe",
            "quickrun",
            transform_json,
            "-o",
            result_json,
            "-d",
            working_dir,
        ]

        proc = subprocess.Popen(cmd)

        last_report = -1

        while proc.poll() is None:

            progress = _get_openfe_progress(
                working_dir
            )

            if (
                progress is not None
                and int(progress) != last_report
            ):

                logger.info(
                    f"[{i}/{total}] "
                    f"{stem}: "
                    f"{progress:.1f}%"
                )

                last_report = int(progress)

            time.sleep(300)

        if proc.returncode != 0:

            logger.error(
                f"[{i}/{total}] "
                f"FAILED {stem}"
            )

            failed.append(stem)

            continue

        logger.info(
            f"[{i}/{total}] "
            f"Completed {stem}"
        )

    logger.info(
        f"Completed "
        f"{total - len(failed)} / {total} "
        f"transformations"
    )

    # deal with failed jobs

    if failed:

        logger.warning(
            f"{len(failed)} transformations failed"
        )

        for name in failed:

            logger.warning(
                f"FAILED: {name}"
            )

    failed_file = os.path.join(
        output_dir,
        "failed_transformations.txt"
    )

    with open(failed_file, "w") as f:

        for name in failed:
            f.write(f"{name}\n")

    self.openfe_results_dir = output_dir

@register_task(
    "openfe_gather_results",
    category="Free Energy",
    description="Gather OpenFE results."
)
def openfe_gather_results(self):

    cfg = self.config["openfe_gather_results"]

    results_dir = cfg["results_dir"]

    output_file = cfg.get(
        "output_file",
        "final_results.tsv"
    )

    report = cfg.get(
        "report",
        "dg"
    )

    cmd = [
        "openfe",
        "gather",
        results_dir,
        "--report",
        report,
        "-o",
        output_file,
    ]

    subprocess.run(
        cmd,
        check=True
    )

    logger.info(
        f"Wrote gathered results: "
        f"{output_file}"
    )

    self.openfe_results_table = output_file