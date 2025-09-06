import os
import importlib
import yaml
from backends.psi4 import Psi4Backend
from backends.xtb import XTBBackend
from modules import load_all_tasks
from workflows import load_all_workflows
# from backends.orca import OrcaBackend

from task_registry import get_task, get_unregistered_tasks_used
from workflows import get_workflow

load_all_tasks()
load_all_workflows()

# Helpers for importing all available modules and workflows

def auto_import_modules(package_dir='modules'):
    for filename in os.listdir(package_dir):
        if filename.endswith('.py') and not filename.startswith('_'):
            module_name = filename[:-3]
            importlib.import_module(f"{package_dir}.{module_name}")
def auto_import_workflows(package_dir='workflows'):
    for filename in os.listdir(package_dir):
        if filename.endswith('.py') and not filename.startswith('_'):
            module_name = filename[:-3]
            importlib.import_module(f"{package_dir}.{module_name}")

# The actual pipeline

class QuantumPipeline:
    def __init__(self, xyz_file: str, config_file: str):
        self.xyz_file = xyz_file
        self.config = self.load_config(config_file)
        self.backend = self.init_backend()

    def load_config(self, file_path):
        with open(file_path, 'r') as f:
            return yaml.safe_load(f)

    def init_backend(self):
        backend_name = self.config['backend'].lower()
        if backend_name == 'psi4':
            return Psi4Backend(self.config)
        elif backend_name == 'xtb':
            return XTBBackend(self.config)
        # elif backend_name == 'orca':
        #     return OrcaBackend(self.config)
        else:
            raise ValueError(f"Unsupported backend: {backend_name}")

    def run(self):
        task = self.config['task'].lower()

        if task in ['optimize', 'single_point']:
            self.backend.run(self.xyz_file, task)

        elif task == 'workflow':
            wf_name = self.config.get('workflow_name', 'default')
            wf_func = get_workflow(wf_name)
            if not wf_func:
                raise ValueError(f"No workflow named '{wf_name}' found.")
            
            print(f"[Pipeline] Running workflow '{wf_name}'...")
            wf_func(self.backend, self.xyz_file, self.config)

            # warn about unregistered tasks
            from task_registry import get_unregistered_tasks_used
            unreg = get_unregistered_tasks_used()
            if unreg:
                print(f"\n[Warning] The following tasks were called but not registered:")
                for task in unreg:
                    print(f" - {task}")
                raise RuntimeError("Workflow tried to use some unregistered tasks. Fix or register them.")

        else:
            task_func = get_task(task)
            if task_func:
                task_func(self.backend, self.xyz_file, self.config)
            else:
                raise ValueError(f"Unsupported task: {task}")

if __name__ == '__main__':
    import sys
    xyz_file, config_file = sys.argv[1], sys.argv[2]
    pipeline = QuantumPipeline(xyz_file, config_file)
    pipeline.run()

