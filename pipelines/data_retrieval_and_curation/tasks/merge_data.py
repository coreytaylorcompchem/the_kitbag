from .base import BaseTask
from pipeline.task_registry import register_task
from models.activity import Activity

@register_task("merge_bioactivity")
class MergeBioactivityTask(BaseTask):
    def run(self, data):
        combined = []
        for dataset in data:
            combined.extend(dataset)
        return combined
