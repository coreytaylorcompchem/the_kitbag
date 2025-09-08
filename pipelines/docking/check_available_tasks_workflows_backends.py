from modules import load_all_tasks
from workflows import load_all_workflows, list_workflows, get_workflow_metadata
from backends import load_all_backends, list_backends
from docking_task_registry import list_tasks, get_task_metadata

def main():
    load_all_tasks()
    load_all_workflows()
    load_all_backends()

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

if __name__ == "__main__":
    main()
