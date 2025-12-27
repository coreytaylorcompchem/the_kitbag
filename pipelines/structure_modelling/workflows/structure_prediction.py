from workflows import register_workflow
from pipeline.task_registry import get_task
from pathlib import Path

@register_workflow(
    "structure_prediction",
    description="Parallel AI-based structure prediction"
)
def run(config):
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    for step in config["workflow"]:
        task = get_task(step)
        task(None, config)
