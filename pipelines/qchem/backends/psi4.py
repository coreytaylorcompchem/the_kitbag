import psi4
from backends.base import BaseBackend

class Psi4Backend(BaseBackend):
    def run(self, xyz_file, task):
        if task == "optimize" or task == "optimise":
            return self.optimise(xyz_file, self.config)
        elif task == "single_point":
            return self.single_point(xyz_file, self.config)
        else:
            raise NotImplementedError(f"Task '{task}' not implemented in Psi4 backend.")
    
    def optimise(self, xyz_file, config):
        with open(xyz_file) as f:
            mol_xyz = f.read()

        mol = psi4.geometry(mol_xyz)

        method = config.get('method', 'b3lyp')
        basis = config.get('basis', 'def2-svp')
        maxiter = config.get('maxiter', 50)

        psi4.set_options({
            'basis': basis,
            'maxiter': maxiter,
            'scf_type': config.get('scf_type', 'pk'),
            'reference': config.get('reference', 'rhf'),
        })

        print(f"[Psi4] Optimizing with {method}/{basis}")
        energy, wfn = psi4.optimize(f"{method}/{basis}", molecule=mol, return_wfn=True)

        return wfn.molecule().save_string_xyz()

