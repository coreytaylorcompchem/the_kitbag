from modules import load_all_tasks
from workflows import load_all_workflows, list_workflows, get_workflow_metadata
from backends import list_backends
from task_registry import list_tasks, get_task_metadata

load_all_tasks()
load_all_workflows()

print("Available Backends:")
for backend in list_backends():
    print(f" - {backend}")

print("\nAvailable Tasks:")
for task in list_tasks():
    meta = get_task_metadata(task)
    desc = meta.get('description', '')
    supported = ', '.join(meta.get('supported_backends', []))
    print(f" - {task}: {desc}")
    print(f"   ↳ Backends: {supported or 'None'}")

print("\nAvailable Workflows:")
for wf in list_workflows():
    meta = get_workflow_metadata(wf)
    desc = meta.get('description', '')
    print(f" - {wf}: {desc}")