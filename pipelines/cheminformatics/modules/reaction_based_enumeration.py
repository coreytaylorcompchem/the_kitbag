from rdkit import Chem
from rdkit.Chem import BRICS
import pandas as pd
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger
# from workflows import register_workflow

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("reaction_based_enumeration", category="Library generation", description="Generate library from fragment with reaction SMARTS (WIP).")
def combinatorial_enumerator(config: dict, data: dict = None) -> dict:
    
    return None
