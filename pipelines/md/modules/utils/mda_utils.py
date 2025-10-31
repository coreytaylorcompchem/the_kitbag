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

    COLORS: ClassVar[Dict[Optional[str], str]] = {
            None: "white",
            **separated_interaction_colors,
    }

    color_mapper = {
        interaction: value for value, interaction in enumerate(COLORS)
    }
    return color_mapper

def _get_inv_color_mapper():
    
    color_mapper = _get_color_mapper()
    
    inv_color_mapper = {
                value: interaction for interaction, value in color_mapper.items()
    }

    return inv_color_mapper

def _bit_to_color_value(s: pd.Series) -> pd.Series:
    """Replaces a bit value with it's corresponding color value"""

    color_mapper = _get_color_mapper()
    
    interaction = s.name[-1]
    return s.apply(
        lambda v: (
            color_mapper[interaction] if v else color_mapper[None]
        ),
    )