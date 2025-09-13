import yaml
from backends.rdkit_backend import RDKitBackend
from registry.workflow_registry import get_workflow

def get_backend(name):
    if name == 'rdkit':
        return RDKitBackend()
    # Add mordred, deepchem backends as needed
    raise ValueError(f"Unknown backend: {name}")

def main(config_file):
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    backend = get_backend(config['backend'])
    workflow_class = get_workflow(config['workflow'])
    workflow = workflow_class(config, backend)
    workflow.setup()
    results = workflow.run(config['input_file'])

    print(f"Workflow complete. Processed {len(results)} molecules.")

if __name__ == "__main__":
    main("configs/example_workflow.yaml")
