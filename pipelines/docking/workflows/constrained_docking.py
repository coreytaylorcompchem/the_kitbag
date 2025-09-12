from workflows import register_workflow

@register_workflow("constrained_docking", description="Dock with core constraints.")
def run(config_path: str):
    print(f"Running constrained docking with config: {config_path}")
    # Constrained docking code here!