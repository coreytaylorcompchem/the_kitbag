import os
import io
import sys
import shutil
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

        if ram:
            psi4.set_memory(f"{ram} MB")
        if ncpu:
            psi4.set_num_threads(ncpu)

        psi4_opts = {
            'basis': basis,
            'scf_type': config.get('scf_type', 'pk'),
            'reference': config.get('reference', 'rhf'),
        }

        psi4.set_options(psi4_opts)

        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            psi4.core.set_output_file(log_path, False)
        else:
            psi4.core.set_output_file("psi4.log", False)

        logger.info(f"Running single-point calculation with {method}/{basis}")
        energy, wfn = psi4.energy(f"{method}/{basis}", molecule=mol, return_wfn=True)

        # Save geometry
        output_cfg = config.get('output', {})
        geometry_filename = output_cfg.get('geometry')
        if geometry_filename:
            output_dir = os.path.dirname(log_path) if log_path else "."
            geometry_path = os.path.join(output_dir, geometry_filename)
            with open(geometry_path, "w") as xyz_out:
                xyz_out.write(wfn.molecule().save_string_xyz())
            logger.info(f"Saved final geometry to: {geometry_path}")

        return {
            "energy": energy,
            "wfn": wfn,
            "method": method,
            "basis": basis
        }


    def optimise(self, xyz_file, config, log_path=None):
        with open(xyz_file) as f:
            mol_xyz = f.read()
        mol = psi4.geometry(mol_xyz)

        method = config.get('method', 'b3lyp')
        basis = config.get('basis', 'def2-svp')
        maxiter = config.get('maxiter', 50)
        ram = config.get('ram')
        ncpu = config.get('ncpu')

        if ram:
            psi4.set_memory(f"{ram} MB")
        if ncpu:
            psi4.set_num_threads(ncpu)

        psi4_opts = {
            'basis': basis,
            'scf_type': config.get('scf_type', 'pk'),
            'reference': config.get('reference', 'rhf'),
            'maxiter': maxiter,
        }

        psi4.set_options(psi4_opts)

        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            psi4.core.set_output_file(log_path, False)
        else:
            psi4.core.set_output_file("psi4_opt.log", False)

        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()
        sys_stdout_backup = sys.stdout
        sys_stderr_backup = sys.stderr
        sys.stdout = stdout_capture
        sys.stderr = stderr_capture

        try:
            energy, wfn = psi4.optimize(f"{method}/{basis}", molecule=mol, return_wfn=True)
        finally:
            sys.stdout = sys_stdout_backup
            sys.stderr = sys_stderr_backup

        with open(log_path, "a") as log_file:
            log_file.write("\n[Python stdout]\n")
            log_file.write(stdout_capture.getvalue())
            log_file.write("\n[Python stderr]\n")
            log_file.write(stderr_capture.getvalue())

        # Save geometry 
        output_cfg = config.get('output', {})
        geometry_filename = output_cfg.get('geometry')
        if geometry_filename:
            output_dir = os.path.dirname(log_path) if log_path else "."
            geometry_path = os.path.join(output_dir, geometry_filename)
            with open(geometry_path, "w") as xyz_out:
                xyz_out.write(wfn.molecule().save_string_xyz())
            logger.info(f"Saved final geometry to: {geometry_path}")

        return {
            "energy": energy,
            "wfn": wfn,
            "method": method,
            "basis": basis
        }


    def mesp(self, xyz_or_wfn, config, mesp_cube_path, density_cube_path, log_path=None):
        # Ensure output dir exists
        output_dir = os.path.dirname(os.path.abspath(mesp_cube_path))
        os.makedirs(output_dir, exist_ok=True)

        # Prepare wavefunction
        if hasattr(xyz_or_wfn, 'molecule'):
            wfn = xyz_or_wfn
            mol = wfn.molecule()

            # Try to pull method and basis from wfn if not provided
            method = config.get('method', None)
            basis = config.get('basis', None)
            if method is None or basis is None:
                try:
                    method = method or wfn.functional().name()
                except Exception:
                    method = method or 'unknown_method'
                try:
                    basis = basis or wfn.basisset().name()
                except Exception:
                    basis = basis or 'unknown_basis'
        else:
            with open(xyz_or_wfn) as f:
                mol_xyz = f.read()
            mol = psi4.geometry(mol_xyz)

            method = config.get('method', 'b3lyp')
            basis = config.get('basis', 'def2-svp')
            ram = config.get('ram')
            ncpu = config.get('ncpu')

            if ram:
                psi4.set_memory(f"{ram} MB")
            if ncpu:
                psi4.set_num_threads(ncpu)

            psi4.set_options({
                'basis': basis,
                'scf_type': config.get('scf_type', 'pk'),
                'reference': config.get('reference', 'rhf'),
            })

            if log_path:
                os.makedirs(os.path.dirname(log_path), exist_ok=True)
                psi4.core.set_output_file(log_path, False)
            else:
                psi4.core.set_output_file("psi4_mesp.log", False)

            logger.info(f"Running energy calculation to obtain wavefunction for MESP with {method}/{basis}")
            energy, wfn = psi4.energy(f"{method}/{basis}", molecule=mol, return_wfn=True)

        # Log method/basis used
        logger.info(f"Running cubeprop with method: {method}, basis: {basis}")

        # Set cubeprop options to only do ESP and DENSITY
        psi4.set_options({
            "CUBEPROP_TASKS": ["ESP", "DENSITY"],
            "CUBEPROP_FILEPATH": output_dir
        })

        prev_cwd = os.getcwd()
        os.chdir(output_dir)

        try:
            psi4.cubeprop(wfn)
        except Exception as e:
            logger.error(f"Error during cubeprop: {e}")
        finally:
            os.chdir(prev_cwd)

        # Move generated files
        esp_default = os.path.join(output_dir, "esp.cube")
        density_default = os.path.join(output_dir, "density.cube")

        for default_name, target_path in [(esp_default, mesp_cube_path), (density_default, density_cube_path)]:
            if os.path.exists(default_name):
                if os.path.abspath(default_name) != os.path.abspath(target_path):
                    shutil.move(default_name, target_path)
                logger.info(f"Saved cube file: {target_path}")
            else:
                logger.warning(f"Expected cube file {default_name} not found!")

        if log_path:
            psi4.core.set_output_file("stdout")

        return None

