import numpy as np

from openmm.app import (
    PDBFile,
    Modeller,
)

from openfe import (
    ProteinMembraneComponent,
)

from openff.units import unit

from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

from pipeline.task_registry import register_task

from tqdm import tqdm

from backends.openmm_backend import (
    OpenMMBackend,
)

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def _prepare_openfe_receptor(pdb_in, pdb_out):

    pdb = PDBFile(pdb_in)
    
    box = pdb.topology.getPeriodicBoxVectors()

    modeller = Modeller(
        pdb.topology,
        pdb.positions
    )

    backend = OpenMMBackend({})
    modeller = backend.cap_internal_chain_breaks(modeller)

    # Need this after capping otherwise box vectors get lost.
    
    if box is not None:
        modeller.topology.setPeriodicBoxVectors(box)

    with open(pdb_out, "w") as f:
        PDBFile.writeFile(
            modeller.topology,
            modeller.positions,
            f,
            keepIds=True,
        )

    return pdb_out

def _check_protein_connectivity(pdb_file):

    pdb = PDBFile(pdb_file)
    topology = pdb.topology

    logger.debug(
        f"Topology contains "
        f"{topology.getNumChains()} chains"
    )

    for chain in topology.chains():

        residues = list(chain.residues())

        logger.debug(
            f"Chain {chain.id}: "
            f"{len(residues)} residues"
        )

        for r1, r2 in zip(
            residues[:-1],
            residues[1:]
        ):

            c_atom = None
            n_atom = None

            for atom in r1.atoms():

                if atom.name == "C":
                    c_atom = atom
                    break

            for atom in r2.atoms():

                if atom.name == "N":
                    n_atom = atom
                    break

            if (
                c_atom is None or
                n_atom is None
            ):
                continue

            pos1 = pdb.positions[
                c_atom.index
            ]

            pos2 = pdb.positions[
                n_atom.index
            ]

            dist = np.linalg.norm(
                np.array(
                    [
                        pos1.x - pos2.x,
                        pos1.y - pos2.y,
                        pos1.z - pos2.z,
                    ]
                )
            )

            if dist > 0.25:

                logger.warning(
                    f"Large peptide gap: "
                    f"chain={chain.id} "
                    f"{r1.name}{r1.id} -> "
                    f"{r2.name}{r2.id} "
                    f"distance={dist:.3f} nm"
                )

def _has_valid_cryst1(pdb_file):

    try:
        pdb = PDBFile(pdb_file)
        box = pdb.topology.getPeriodicBoxVectors()

        if box is None:
            return False

        return True

    except Exception:

        return False

def _load_protein_membrane_component(
    pdb_file,
):

    _check_protein_connectivity(
        pdb_file
    )

    has_box = _has_valid_cryst1(
        pdb_file
    )

    logger.info(
        f"Periodic box present: {has_box}"
    )

    if has_box:

        logger.info(
            "Using periodic box vectors "
            "from PDB"
        )

        return (
            ProteinMembraneComponent
            .from_pdb_file(
                pdb_file
            )
        )

    logger.warning(
        "No periodic box detected. "
        "Inferring box vectors."
    )

    return (
        ProteinMembraneComponent
        .from_pdb_file(
            pdb_file,
            infer_box_vectors=True,
        )
    )

def _build_protocol_settings(
    protocol_cfg,
    leg,
):

    settings = (
        RelativeHybridTopologyProtocol
        .default_settings()
    )

    #
    # repeats
    #

    settings.protocol_repeats = (
        protocol_cfg.get(
            "repeats",
            1
        )
    )

    #
    # simulation lengths
    #

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

    #
    # lambda schedule
    #

    if "lambda_windows" in protocol_cfg:

        n = protocol_cfg["lambda_windows"]

        settings.lambda_settings.lambda_windows = n
        settings.simulation_settings.n_replicas = n

    #
    # replica exchange frequency
    #

    if "swap_frequency" in protocol_cfg:

        settings.simulation_settings.n_replicas_exchange_attempts = (
            protocol_cfg["swap_frequency"]
        )

    #
    # output frequency
    #

    if "checkpoint_interval_ns" in protocol_cfg:

        logger.debug(
            settings.output_settings
        )

        settings.output_settings.checkpoint_interval = (
            protocol_cfg["checkpoint_interval_ns"]
            * unit.nanosecond
        )

    if "positions_write_frequency_ps" in protocol_cfg:

        settings.output_settings.positions_write_frequency = (
            protocol_cfg["positions_write_frequency_ps"]
            * unit.picosecond
        )

    #
    # platform
    #

    if "platform" in protocol_cfg:

        settings.engine_settings.compute_platform = (
            protocol_cfg["platform"]
        )

        logger.debug(f"engine_settings: {settings.engine_settings}")
        logger.debug(f"dir(engine_settings): {dir(settings.engine_settings)}")
        logger.debug(f"engine_settings__dict__: {settings.engine_settings.__dict__}")
        logger.debug(f"gpu_device_index: {settings.engine_settings.gpu_device_index}")

    #
    # membrane systems only
    #

    if leg == "complex":

        settings.integrator_settings.barostat = (
            "MonteCarloMembraneBarostat"
        )

    return settings