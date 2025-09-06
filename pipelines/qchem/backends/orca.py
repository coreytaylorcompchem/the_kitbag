from orca_interface import Orca 
from backends.base import BaseBackend

class OrcaBackend(BaseBackend):
    def run(self, xyz_file: str, task: str):
        orca = Orca()
        job = orca.create_job(xyz_file, method=self.config['method'], basis=self.config['basis'])

        if task == 'optimize':
            job.optimize()
        elif task == 'single_point':
            job.single_point()
