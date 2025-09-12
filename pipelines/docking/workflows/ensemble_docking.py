from workflows import register_workflow

@register_workflow("ensemble_docking", description="Preparation, dock and score emsemble of proteins.")
def run(config_path: str):
    print(f"Running ensemble docking with config: {config_path}")
    # Ensemble docking code here!