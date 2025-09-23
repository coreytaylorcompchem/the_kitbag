import os
import subprocess
import shutil
from backends.base import BaseBackend
from pipeline.logger import setup_logger

# For accessing package data
import pkg_resources

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

class MultiwfnBackend(BaseBackend):
    def __init__(self, config=None):
        self.config = config
        self.executable_path = "/home/corey/local_software/Multiwfn_3.8_dev_bin_Linux/Multiwfn"
        if not os.path.exists(self.executable_path):
            raise FileNotFoundError(f"Multiwfn executable not found at: {self.executable_path}")

    def run_multiwfn(self, input_file, output_dir, ncpu=4, log_path=None, script_input=None):
        # Ensure input file is in the output directory
        input_name = os.path.basename(input_file)
        input_path_in_output = os.path.join(output_dir, input_name)
        
        ncpu = self.config.get('qtaim', {}).get('ncpu', 4)

        if not os.path.samefile(input_file, input_path_in_output):
            shutil.copyfile(input_file, input_path_in_output)
            
        # Write custom settings.ini to boost the number of cores for intensive calculations
        self._generate_settings_ini(ncpu, output_dir)
        # Build the mwfn command
        command = [self.executable_path, input_name]
        logger.debug(f"Executing Multiwfn: {' '.join(command)} in {output_dir}")

        process = subprocess.Popen(
            command,
            cwd=output_dir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(script_input)

        if log_path:
            with open(log_path, "w") as log_file:
                log_file.write(stdout)
                log_file.write("\n[stderr]\n")
                log_file.write(stderr)
        
        if process.returncode != 0:
            raise RuntimeError(f"Multiwfn failed with exit code {process.returncode}")

    def _generate_settings_ini(self, ncpu: int, output_dir: str):
        """
        Copy the template settings.ini file and replace the nthreads line
        with the user-specified number of threads (ncpu).
        """
        try:
            template_path = pkg_resources.resource_filename(
                'backends.data', 'settings.ini'
            )
        except ModuleNotFoundError:
            raise FileNotFoundError("settings.ini template not found in backends.data")

        with open(template_path, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.strip().startswith("nthreads="):
                if "//" in line:
                    comment = line.split("//", 1)[1]
                    lines[i] = f"  nthreads=  {ncpu} //{comment}\n"
                else:
                    lines[i] = f"  nthreads=  {ncpu}\n"
                break
        
        dest_path = os.path.join(output_dir, "settings.ini")
        with open(dest_path, 'w') as f:
            f.writelines(lines)

        logger.debug(f"Generated settings.ini with nthreads={ncpu} at {dest_path}")
