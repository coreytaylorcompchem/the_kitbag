from workflows.base_workflow import BaseWorkflow
from registry.wokrkflow_registry import register_workflow
from registry.task_registry import get_task

@register_workflow('early_hit_discovery')
class EarlyHitWorkflow(BaseWorkflow):
    def setup(self):
        for task_conf in self.config['tasks']:
            TaskClass = get_task(task_conf['name'])
            task = TaskClass(task_conf['params'], backend=self.backend)
            self.tasks.append(task)
