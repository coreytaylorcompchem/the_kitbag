import os
import io
import sys
import psi4
from backends.base import BaseBackend

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class Psi4Backend(BaseBackend):

    def single_point(self, xyz_file, config, log_path=None):
        with open(xyz_file) as f:
            mol_xyz = f.read()

        mol = psi4.geometry(mol_xyz)

        method = config.get('method', 'b3lyp')
        basis = config.get('basis', 'def2-svp')

        ram = config.get('ram')
        ncpu = config.get('ncpu')

        psi4_opts = {
            'basis': basis,
            'scf_type': config.get('scf_type', 'pk'),
            'reference': config.get('reference', 'rhf'),
        }

        if ram:
            psi4_opts['memory'] = f"{ram} MB"
        if ncpu:
            psi4_opts['num_threads'] = ncpu

        # Ensure output directory exists
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            psi4.core.set_output_file(log_path, False)  # redirect Psi4 log
        else:
            psi4.core.set_output_file("psi4.log", False)

        logger.info(f"Running single-point calculation with {method}/{basis}")
        energy, wfn = psi4.energy(f"{method}/{basis}", molecule=mol, return_wfn=True)

        return energy, wfn.molecule().save_string_xyz()

    def optimise(self, xyz_file, config, log_path=None):
        with open(xyz_file) as f:
            mol_xyz = f.read()

        mol = psi4.geometry(mol_xyz)

        method = config.get('method', 'b3lyp')
        basis = config.get('basis', 'def2-svp')
        maxiter = config.get('maxiter', 50)

        ram = config.get('ram')
        ncpu = config.get('ncpu')

        psi4_opts = {
            'basis': basis,
            'scf_type': config.get('scf_type', 'pk'),
            'reference': config.get('reference', 'rhf'),
        }

        if ram:
            psi4_opts['memory'] = f"{ram} MB"
        if ncpu:
            psi4_opts['num_threads'] = ncpu

        # Prepare output redirection
        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            psi4.core.set_output_file(log_path, False)
        else:
            psi4.core.set_output_file("psi4_opt.log", False)

        # Capture Python stdout/stderr to ensure all prints are logged
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()
        sys_stdout_backup = sys.stdout
        sys_stderr_backup = sys.stderr
        sys.stdout = stdout_capture
        sys.stderr = stderr_capture

        try:
            # print(f"[Psi4] Running geometry optimization with {method}/{basis}")
            energy, wfn = psi4.optimize(f"{method}/{basis}", molecule=mol, return_wfn=True)
        finally:
            sys.stdout = sys_stdout_backup
            sys.stderr = sys_stderr_backup

        # Append captured stdout and stderr to the log file
        with open(log_path, "a") as log_file:
            log_file.write("\n[Python stdout]\n")
            log_file.write(stdout_capture.getvalue())
            log_file.write("\n[Python stderr]\n")
            log_file.write(stderr_capture.getvalue())

        # Return energy and optimized geometry
        return energy, wfn.molecule().save_string_xyz()

