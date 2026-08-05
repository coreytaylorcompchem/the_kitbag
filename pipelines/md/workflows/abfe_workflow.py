from workflows import register_workflow

from pipeline.task_registry import (
    get_task,
    list_tasks,
)

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True,
)


class ABFEWorkflow:
    """
    Lightweight workflow context object.

    Tasks receive this instance as `self`,
    allowing them to share state through attributes:

        self.openfe_network
        self.alchemical_network
        self.receptor_pdb

    etc.
    """

    def __init__(self, config):
        self.config = config


@register_workflow(
    "abfe",
    description="Prepare and build OpenFE RBFE calculations."
)
def run_abfe_workflow(config):

    logger.info(
        f"Running workflow: "
        f"{config.get('workflow_name', 'abfe')}"
    )

    logger.info(
        f"Workflow steps: "
        f"{config.get('workflow', [])}"
    )

    workflow = ABFEWorkflow(config)

    context = {}

    for task_name in config.get("workflow", []):

        task_func = get_task(task_name)

        if task_func is None:

            logger.warning(
                f"Task '{task_name}' not found."
            )

            logger.warning(
                f"Available tasks: "
                f"{list_tasks()}"
            )

            continue

        logger.info(
            f"Running task "
            f"'{task_name}'..."
        )

        result = task_func(workflow)

        context[task_name] = result

        logger.info(
            f"Task '{task_name}' completed."
        )

    return context