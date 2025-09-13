from modules.base_task import BaseTask
from registry.task_registry import register_task

@register_task('physchem_filtering')
class PhyschemFiltering(BaseTask):
    def run(self, molecules):
        filtered = []
        for mol in molecules:
            if mol is None:
                continue
            if self.backend.filter_molecule(mol, mw_cutoff=self.config.get("mw_cutoff", 500)):
                filtered.append(mol)
        return filtered
