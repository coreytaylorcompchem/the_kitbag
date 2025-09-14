from tasks import load_all_tasks
from workflows import list_workflows, get_workflow_metadata, load_all_workflows
#from backends import list_backends # later
from pipeline.task_registry import list_tasks, get_task_metadata

import sys
import os

from pipeline.logger import setup_logger

logger = setup_logger(
    __name__,
    debug_mode=False,
    simple_format=True
)

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

def print_pipeline_capabilities():
    load_all_tasks()
    load_all_workflows()

   # logger.info("Available Backends:")
   # for backend in list_backends():
   #     logger.info(f" - {backend}")

    logger.info("\nAvailable Tasks:")
    for task in list_tasks():
        meta = get_task_metadata(task) or {}
        desc = meta.get('description', '')
        supported = ', '.join(meta.get('supported_backends', []))
        logger.info(f" - {task}: {desc}")
       # logger.info(f"   ↳ Backends: {supported or 'None'}")

    logger.info("\nAvailable Workflows:")
    for wf in list_workflows():
        meta = get_workflow_metadata(wf) or {}
        desc = meta.get('description', '')
        logger.info(f" - {wf}: {desc}")

if __name__ == "__main__":
    print_pipeline_capabilities()

