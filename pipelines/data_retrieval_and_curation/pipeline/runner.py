import yaml
from pipeline.task_registry import TASK_REGISTRY
from backends.chembl import ChemblBackend
from backends.pubchem import PubChemBackend  # Example placeholder

BACKENDS = {
    "chembl": ChemblBackend,
    "pubchem": PubChemBackend,
    # Add more backends as needed
}

class PipelineRunner:
    def __init__(self, config_path):
        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        self.backend_instances = {
            name: BACKENDS[name]() for name in self.config['backends']
        }

    def run(self):
        # Step 1: Fetch data
        data = []
        for backend_name, compound_list in self.config['compounds'].items():
            backend = self.backend_instances[backend_name]
            for compound_id in compound_list:
                activities = backend.fetch_bioactivity(compound_id)
                # Add source label for traceability
                for act in activities:
                    act['source'] = backend_name
                data.append(activities)

        # Step 2: Run Tasks
        task_data = data
        for task_cfg in self.config['tasks']:
            task_cls = TASK_REGISTRY[task_cfg['name']]
            task = task_cls(task_cfg.get('params', {}))
            task_data = task.run(task_data)

        return task_data
