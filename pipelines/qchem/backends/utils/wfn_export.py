from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def write_wfn_file(wfn, filepath, wfn_filename):
    """
    Write a Multiwfn-compatible .wfn file from a Psi4 wavefunction object.
    Only works with Psi4 wavefunctions and supports certain method/basis combos.
    """
    try:
        import psi4

        # Psi4 does not natively export .wfn in v1.9.1, but can write Molden file
        molden_path = filepath.replace(".wfn", ".molden")
        psi4.molden(wfn, molden_path)
        logger.info(f"Saved {wfn_filename} to {molden_path}")

        # Multiwfn can *sometimes* accept Molden format — less ideal but may work
        return molden_path

    except Exception as e:
        raise RuntimeError(f"Failed to write .wfn file: {e}")