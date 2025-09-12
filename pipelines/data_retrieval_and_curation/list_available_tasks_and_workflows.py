from tasks import load_all_tasks
from workflows import list_workflows, get_workflow_metadata, load_all_workflows
#from backends import list_backends # later
from pipeline.task_registry import list_tasks, get_task_metadata

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

def print_pipeline_capabilities():
    load_all_tasks()
    load_all_workflows()

   # print("Available Backends:")
   # for backend in list_backends():
   #     print(f" - {backend}")

    print("\nAvailable Tasks:")
    for task in list_tasks():
        meta = get_task_metadata(task) or {}
        desc = meta.get('description', '')
        supported = ', '.join(meta.get('supported_backends', []))
        print(f" - {task}: {desc}")
       # print(f"   ↳ Backends: {supported or 'None'}")

    print("\nAvailable Workflows:")
    for wf in list_workflows():
        meta = get_workflow_metadata(wf) or {}
        desc = meta.get('description', '')
        print(f" - {wf}: {desc}")

if __name__ == "__main__":
    print_pipeline_capabilities()

