import matplotlib as plt
import numpy as np
import pandas as pd

from typing import TYPE_CHECKING, ClassVar, Dict, List, Literal, Optional, Tuple

import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import rdDepictor
from matplotlib.colors import ListedColormap
from prolif.plotting.utils import separated_interaction_colors

def _get_color_mapper():
    interactions = [None] + list(separated_interaction_colors.keys())
    return {name: i for i, name in enumerate(interactions)}

def _get_inv_color_mapper():
    cmap = _get_color_mapper()
    return {v: k for k, v in cmap.items()}

def _bit_to_color_value(s: pd.Series) -> pd.Series:
    """Replaces a bit value with it's corresponding color value"""

    color_mapper = _get_color_mapper()
    
    interaction = s.name[-1]
    return s.apply(
        lambda v: (
            color_mapper[interaction] if v else color_mapper[None]
        ),
    )