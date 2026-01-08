from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED

from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)


def compute_descriptors_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "csp3": Descriptors.FractionCSP3(mol),
        "QED": QED.qed(mol),
    }


@register_task(
    "compute_rdkit_descriptors",
    category="Annotation",
    description="Compute RDKit molecular descriptors per ligand.",
)
def compute_rdkit_descriptors(backend, ligand, config, **kwargs):
    """
    Ligand-level RDKit descriptors.
    Runs once per ligand and attaches values directly to ligand dict.
    """

    smiles = ligand.get("smiles")
    if not smiles:
        raise ValueError(f"Ligand '{ligand.get('name')}' missing SMILES.")

    try:
        desc = compute_descriptors_from_smiles(smiles)
        if desc is None:
            logger.warning(f"[RDKit] Failed to compute descriptors for {ligand['name']}")
            return

        ligand.update(desc)

        logger.debug(
            f"[RDKit] Descriptors computed for {ligand['name']}: "
            f"{', '.join(desc.keys())}"
        )

    except Exception as e:
        logger.error(f"[RDKit] Failed for {ligand['name']}: {e}")
