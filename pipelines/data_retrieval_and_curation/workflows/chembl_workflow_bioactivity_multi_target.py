from workflows import register_workflow
from pipeline.parallel_runner import ParallelWorkflowRunner
from pipeline.task_registry import get_task


def workflow_function(config):
    steps = config.get("workflow", [])
    data = None

    for step in steps:
        task_func = get_task(step)
        if not task_func:
            raise ValueError(f"Task '{step}' not found in registry.")
        data = task_func(config, data)

    return data


@register_workflow(
    "chembl_multi_target",
    description="Retrieve, standardise and collate bioactivities for multiple targets - CHEMBL."
)
def run_chembl_multi_target_parallel_workflow(config):
    runner = ParallelWorkflowRunner(
        workflow_func=workflow_function,
        config=config,
        input_key="uniprot_ids",
        output_key="uniprot_id",
        filename_pattern="{uniprot_id}_bioactivity.csv",
        combined_filename=config.get("output", {}).get("filename", "combined_bioactivity.csv"),
        output_dir=config.get("output", {}).get("directory", "outputs/parallel")
    )

    return runner.run()
