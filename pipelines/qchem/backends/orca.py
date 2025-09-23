import os
import tempfile
import shutil
# from orca import OrcaCalculation, InputFile
from backends.base import BaseBackend

from pipeline.logger import setup_logger
logger = setup_logger(__name__, debug_mode=True, simple_format=True)


class OrcaBackend(BaseBackend):

    def _write_xyz_to_tmp(self, xyz_str):
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.xyz', mode='w')
        tmp.write(xyz_str)
        tmp.close()
        return tmp.name

    def _make_input_file(self, xyz_path, config, calc_type='opt'):

        method = config.get('method', 'B3LYP')
        basis = config.get('basis', 'def2-SVP')
        ncpu = config.get('ncpu', 2)
        solvent = config.get('solvent', None)

        inp = InputFile()
        inp.set_xyz_from_file(xyz_path)

        # Route section
        route = f"! {method} {basis} {calc_type} TightSCF Grid4"
        if solvent:
            # Add implicit solvent model (e.g., CPCM)
            route += f" CPCM({solvent})"
        inp.route = route

        # Additional settings
        inp.settings["%pal"] = f"nprocs {ncpu}"
        inp.settings["%output"] = 'Print[P_MOs] 1'

        return inp

    def optimise(self, xyz_file, config, log_path=None):
        logger.info("Running ORCA optimization...")

        inp = self._make_input_file(xyz_file, config, calc_type='opt')
        calc = OrcaCalculation(inp)

        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            calc.set_output_file(log_path)

        calc.run()

        energy = calc.extract_energy()
        final_geom = calc.get_final_structure().to_xyz_string()

        # Optional save geometry
        output_cfg = config.get('output', {})
        geometry_filename = output_cfg.get('geometry')
        if geometry_filename:
            output_dir = os.path.dirname(log_path) if log_path else "."
            geometry_path = os.path.join(output_dir, geometry_filename)
            with open(geometry_path, "w") as xyz_out:
                xyz_out.write(final_geom)
            logger.info(f"Saved final geometry to: {geometry_path}")

        return {
            "energy": energy,
            "method": config.get('method', 'B3LYP'),
            "basis": config.get('basis', 'def2-SVP'),
        }

    def single_point(self, xyz_file, config, log_path=None):
        logger.info("Running ORCA single point...")

        inp = self._make_input_file(xyz_file, config, calc_type='sp')
        calc = OrcaCalculation(inp)

        if log_path:
            os.makedirs(os.path.dirname(log_path), exist_ok=True)
            calc.set_output_file(log_path)

        calc.run()

        energy = calc.extract_energy()
        return {
            "energy": energy,
            "method": config.get('method', 'B3LYP'),
            "basis": config.get('basis', 'def2-SVP'),
        }