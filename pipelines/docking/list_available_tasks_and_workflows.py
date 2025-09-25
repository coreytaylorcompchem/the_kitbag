from collections import defaultdict

from modules import load_all_tasks
from workflows import load_all_workflows, list_workflows, get_workflow_metadata
from backends import discover_backends, list_backends
from pipeline.docking_task_registry import (
    list_tasks, get_task_metadata, finalize_task_registration
)

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    spartan_format=True
)

def main():
    load_all_tasks()
    load_all_workflows()
    discover_backends()
    finalize_task_registration() 

    logger.info("\nAvailable Tasks:")

    # organise tasks by category
    tasks_by_category = defaultdict(list)
    for task_name in list_tasks():
        meta = get_task_metadata(task_name) or {}
        category = meta.get('category', 'Uncategorized')
        tasks_by_category[category].append((task_name, meta))

    # print tasks grouped by category
    for category in sorted(tasks_by_category):
        logger.info(f"  [{category}]:")
        for task_name, meta in sorted(tasks_by_category[category]):
            desc = meta.get('description', '')
            supported = ', '.join(meta.get('supported_backends', []))
            logger.info(f"    - {task_name}: {desc}")
            logger.info(f"       ↳ Backends: {supported or 'None'}")

    logger.info("\nAvailable Workflows:")
    for wf in list_workflows():
        meta = get_workflow_metadata(wf) or {}
        desc = meta.get('description', '')
        logger.info(f" - {wf}: {desc}")

if __name__ == "__main__":
    main()
