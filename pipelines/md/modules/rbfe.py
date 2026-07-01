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

@register_task(
    "openfe_prepare_receptor",
    category="Free Energy",
    description="Prepare membrane receptor for OpenFE."
)
def openfe_prepare_receptor(self):

    cfg = self.config["openfe_prepare_receptor"]

    if "equilibrated_pdb" in cfg:

        u = mda.Universe(
            cfg["equilibrated_pdb"]
        )

    else:

        universes = []

        for key in [
            "protein_pdb",
            "membrane_pdb",
            "water_pdb",
        ]:

            if key in cfg:

                universes.append(
                    mda.Universe(
                        cfg[key]
                    )
                )

        if not universes:

            raise ValueError(
                "No receptor input files provided"
            )

        u = mda.Merge(
            *[
                x.atoms
                for x in universes
            ]
        )

        u.dimensions = universes[0].dimensions.copy()

    ligand_resname = cfg.get("ligand_resname", "LIG")

    output_dir = cfg.get("output_dir", "openfe")
    os.makedirs(output_dir, exist_ok=True)

    receptor = u.select_atoms(
        f"not resname {ligand_resname}"
    )

    receptor_pdb = os.path.join(
        output_dir,
        "receptor_membrane.pdb"
    )

    receptor.write(receptor_pdb)

    fixed_receptor_pdb = os.path.join(
        output_dir,
        "receptor_membrane_fixed.pdb"
    )

    _prepare_openfe_receptor(
        receptor_pdb,
        fixed_receptor_pdb
    )

    self.openfe_receptor_pdb = (
        fixed_receptor_pdb
    )

    logger.info(
        f"Wrote receptor file: {fixed_receptor_pdb}"
    )

    with open(receptor_pdb) as f:
        first_lines = [next(f) for _ in range(20)]

    logger.info(
        f"CRYST1 present: "
        f"{any('CRYST1' in x for x in first_lines)}"
    )

    logger.info(
        f"Saved receptor membrane system "
        f"to {fixed_receptor_pdb}"
    )

@register_task(
    "openfe_create_network",
    category="Free Energy",
    description="Create OpenFE ligand network for membrane RBFE."
)
def openfe_create_network(self):

    cfg = self.config["openfe_create_network"]

    receptor_pdb = getattr(
        self,
        "openfe_receptor_pdb",
        cfg["receptor_pdb"]
    )
    ligand_sdf = cfg["ligand_sdf"]

    output_dir = cfg.get(
        "output_dir",
        "openfe"
    )

    os.makedirs(output_dir, exist_ok=True)

    receptor = _load_protein_membrane_component(
        receptor_pdb
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

    network_cfg = cfg.get(
        "network",
        {}
    )

    max_path_length = network_cfg.get(
        "max_path_length",
        6
    )

    radial = network_cfg.get(
        "radial",
        False
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
                max_path_length=max_path_length,
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

    receptor_pdb = getattr(
        self,
        "openfe_receptor_pdb",
        cfg["receptor_pdb"]
    )

    output_dir = cfg.get(
        "output_dir",
        "openfe"
    )

    os.makedirs(output_dir, exist_ok=True)

    ligand_network = self.openfe_network

    protein = _load_protein_membrane_component(
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

            settings = _build_protocol_settings(
                protocol_cfg,
                leg,
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

    # # set active GPU

    # env = os.environ.copy()

    # gpu_id = (
    #     self.config
    #     .get("openfe_create_alchemical_network", {})
    #     .get("protocol", {})
    #     .get("gpu_id")
    # )

    # logger.debug(
    #     f"gpu_id from config = {gpu_id!r}"
    # )

    # if gpu_id is not None:

    #     env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    #     logger.debug(
    #         f"Using GPU {gpu_id}"
    #     )

    # json_files = sorted(
    #     glob.glob(
    #         os.path.join(
    #             transforms_dir,
    #             "*.json"
    #         )
    #     )
    # )

    # total = len(json_files)

    # logger.info(
    #     f"Found {total} transformations"
    # )

    gpus = cfg.get(
        "gpus",
        [0],
    )

    jobs_per_gpu = cfg.get(
        "jobs_per_gpu",
        1,
    )

    json_files = sorted(
        glob.glob(
            os.path.join(
                transforms_dir,
                "*.json",
            )
        )
    )

    logger.info(
        f"Found {len(json_files)} transformations"
    )

    mp.set_start_method("spawn", force=True)

    manager = mp.Manager()

    queue = manager.Queue()

    failed = manager.list()

    for job in json_files:

        queue.put(job)

    workers = []

    for gpu in gpus:

        for _ in range(jobs_per_gpu):

            p = mp.Process(
                target=_openfe_worker,
                args=(
                    gpu,
                    queue,
                    failed,
                    output_dir,
                ),
            )

            p.start()

            workers.append(p)

            logger.info(
                f"Started worker on GPU {gpu}"
            )

    for p in workers:
        p.join()

    failed = list(failed)

    completed = len(json_files) - len(failed)

    logger.info(
        f"Completed {completed} / {len(json_files)} transformations"
    )

    # for i, transform_json in enumerate(
    #     json_files,
    #     start=1
    # ):

    #     stem = Path(transform_json).stem

    #     result_json = os.path.join(
    #         output_dir,
    #         f"{stem}.json"
    #     )

    #     working_dir = os.path.join(
    #         output_dir,
    #         stem
    #     )

    #     #
    #     # Skip completed calculations
    #     #

    #     if os.path.exists(result_json):

    #         logger.info(
    #             f"[{i}/{total}] "
    #             f"Skipping completed {stem}"
    #         )

    #         continue

    #     logger.info(
    #         f"[GPU {env['CUDA_VISIBLE_DEVICES']}] "
    #         f"Running {stem}"
    #     )

    #     cmd = [
    #         "openfe",
    #         "quickrun",
    #         transform_json,
    #         "-o",
    #         result_json,
    #         "-d",
    #         working_dir,
    #     ]

    #     logger.info(
    #         "CUDA_VISIBLE_DEVICES="
    #         f"{env.get('CUDA_VISIBLE_DEVICES')}"
    #     )

    #     proc = subprocess.Popen(
    #         cmd,
    #         env=env
    #     )

    #     last_report = -1

    #     while proc.poll() is None:

    #         progress = _get_openfe_progress(
    #             working_dir
    #         )

    #         if (
    #             progress is not None
    #             and int(progress) != last_report
    #         ):

    #             logger.info(
    #                 f"[{i}/{total}] "
    #                 f"{stem}: "
    #                 f"{progress:.1f}%"
    #             )

    #             last_report = int(progress)

    #         time.sleep(300)

    #     if proc.returncode != 0:

    #         logger.error(
    #             f"[{i}/{total}] "
    #             f"FAILED {stem}"
    #         )

    #         failed.append(stem)

    #         continue

    #     logger.info(
    #         f"[{i}/{total}] "
    #         f"Completed {stem}"
    #     )

    # logger.info(
    #     f"Completed "
    #     f"{total - len(failed)} / {total} "
    #     f"transformations"
    # )

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

@register_task(
    "openfe_compare_exp_predicted",
    category="Free Energy",
    description="Compare OpenFE RBFE predictions against experiment using Cinnabar."
)
def openfe_compare_exp_predicted(self):

    cfg = self.config["openfe_compare_exp_predicted"]

    #
    # files
    #

    ddg_file = cfg["ddg_file"]
    experimental_file = cfg["experimental_file"]

    #
    # column names
    #

    ligand_column = cfg.get(
        "ligand_column",
        "ligand"
    )

    experimental_column = cfg.get(
        "experimental_column",
        "DG"
    )

    experimental_uncertainty_column = cfg.get(
        "experimental_uncertainty_column",
        None
    )

    output_plot = cfg.get(
        "output_plot",
        "openfe_vs_experiment.png"
    )

    output_csv = cfg.get(
        "output_csv",
        "predicted_absolute_dg.csv"
    )

    #
    # load data
    #

    ddg = pd.read_csv(
        ddg_file,
        sep="\t",
    )

    exp = pd.read_csv(
        experimental_file
    )

    femap = FEMap()

    #
    # add experimental measurements
    #

    for _, row in exp.iterrows():

        uncertainty = (
            row[experimental_uncertainty_column]
            if experimental_uncertainty_column
            else 0.3
        )

        femap.add_experimental_measurement(
            label=row[ligand_column],
            value=row[experimental_column]
                * unit.kilocalorie_per_mole,
            uncertainty=uncertainty
                * unit.kilocalorie_per_mole,
            source="Experiment",
        )

    #
    # add OpenFE RBFE edges
    #

    for _, row in ddg.iterrows():

        femap.add_relative_calculation(
            labelA=row["ligand_i"],
            labelB=row["ligand_j"],
            value=row["DDG(i->j) (kcal/mol)"]
                * unit.kilocalorie_per_mole,
            uncertainty=row["uncertainty (kcal/mol)"]
                * unit.kilocalorie_per_mole,
            source="OpenFE",
        )

    #
    # convert to legacy graph
    #

    graph = femap.to_legacy_graph()

    #
    # save reconstructed absolute values
    #

    absolute = pd.DataFrame(
        [
            {
                "ligand": ligand,
                "DG (kcal/mol)": data["calc_DG"],
                "uncertainty (kcal/mol)": data["calc_dDG"],
            }
            for ligand, data in graph.nodes(data=True)
            if "calc_DG" in data
        ]
    )

    absolute.to_csv(
        output_csv,
        index=False,
    )

    #
    # make comparison plot
    #

    cinnabar_plotting.plot_DGs(
        graph,
        method_name="OpenFE",
        target_name="Experiment",
        title="Experimental vs OpenFE",
        filename=output_plot,
    )

    schrodinger_column = cfg.get(
        "schrodinger_column",
        None,
    )

    if schrodinger_column:

        graph_sch = copy.deepcopy(graph)

        for ligand, data in graph_sch.nodes(data=True):

            match = exp.loc[
                exp[ligand_column] == ligand
            ]

            if match.empty:
                continue

            value = match.iloc[0][schrodinger_column]

            if pd.isna(value):
                continue

            #
            # Replace the OpenFE prediction
            # with the Schrödinger prediction.
            #

            data["calc_DG"] = float(value)

            #
            # Use experimental uncertainty if you
            # have one, otherwise zero.
            #

            data["calc_dDG"] = 0.0

        output_plot_sch = cfg.get(
            "output_plot_schrodinger",
            "schrodinger_vs_experiment.png",
        )

        cinnabar_plotting.plot_DGs(
            graph_sch,
            method_name="Schrodinger FEP+",
            target_name="Experiment",
            title="Experimental vs Schrodinger",
            filename=output_plot_sch,
        )

        logger.info(
            f"Wrote {output_plot_sch}"
        )

    logger.info(
        f"Wrote {output_plot}"
    )

    logger.info(
        f"Wrote {output_csv}"
    )

    self.openfe_absolute_results = output_csv