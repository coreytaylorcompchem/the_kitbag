from rdkit import Chem
from rdkit.Chem import Draw
import os

from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def save_fragment_grid(fragments_counts, title, out_dir, filename="fragments.png"):
    """
    fragments_counts: list of (smiles, count) tuples
    """
    mols = []
    legends = []
    for i, (frag, count) in enumerate(fragments_counts):
        mol = Chem.MolFromSmiles(frag)
        if mol is None:
            continue
        mols.append(mol)
        legend = f"{frag}\nCount: {count}"
        legends.append(legend)

    if not mols:
        logger.warning(f"No valid fragments found for {title}")
        return

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=5,
        subImgSize=(200, 200),
        legends=legends,
        useSVG=False,
        returnPNG=True
    )

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    with open(out_path, "wb") as f:
        f.write(img)

    logger.info(f"{title} image saved to {out_path}")