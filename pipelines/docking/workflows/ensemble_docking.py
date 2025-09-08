from workflows import register_workflow

@register_workflow("ensemble_docking", description="Ensemble docking with Gnina.")
def run(config_path: str):
    print(f"Running ensemble docking with config: {config_path}")
    # Your actual constrained docking logic...