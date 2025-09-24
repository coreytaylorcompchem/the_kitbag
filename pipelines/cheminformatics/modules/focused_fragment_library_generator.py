from pipeline.task_registry import register_task
from rdkit import Chem
from rdkit.Chem import BRICS, MACCSkeys
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Recap import RecapDecompose
from rdkit.Chem import SDWriter
from collections import Counter, defaultdict

import pandas as pd
from pathlib import Path
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

@register_task("focused_fragment_library_generator", category="Library generation", description="Generate a focused library by fragment frequency.")
def focused_fragment_library_generator(config, data=None):
    input_file = config.get("input_file")
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = pd.read_csv(input_path)
    if "smiles" not in df.columns:
        raise ValueError("Input CSV must contain a 'smiles' column.")
    
    fragment_type = config.get("focused_fragment_library", {}).get("fragment_type", "BRICS").upper()
    frequency_threshold = config.get("focused_fragment_library", {}).get("frequency_threshold", 0.3)
    max_fragments = config.get("focused_fragment_library", {}).get("max_fragments", None)
    group_by_target = config.get("focused_fragment_library", {}).get("group_by_target", False)

    output_fragments = config.get("focused_fragment_library", {}).get("output_fragments", True)

    if group_by_target and "target" not in df.columns:
        logger.warning("group_by_target is True but 'target' column not found. Proceeding without grouping.")
        group_by_target = False

    grouped = [("all", df)] if not group_by_target else df.groupby("target")

    results = []

    for group_name, group_df in grouped:
        frag_counter = Counter()
        total_mols = 0

        for smi in group_df["smiles"].dropna().unique():
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            total_mols += 1
            fragments = get_fragments(mol, method=fragment_type)
            frag_counter.update(fragments)

        if total_mols == 0:
            logger.warning(f"No valid molecules for group '{group_name}'")
            continue

        frequent_frags = {
            frag: count / total_mols
            for frag, count in frag_counter.items()
            if (count / total_mols) >= frequency_threshold
        }

        sorted_frags = sorted(frequent_frags.items(), key=lambda x: x[1], reverse=True)
        if max_fragments:
            sorted_frags = sorted_frags[:max_fragments]

        for frag, freq in sorted_frags:
            results.append({
                "group": group_name,
                "fragment_smiles": frag,
                "frequency": round(freq, 4),
                "count": frag_counter[frag],
                "total_molecules": total_mols
            })

    result_df = pd.DataFrame(results)
    logger.info(f"Extracted {len(result_df)} frequent fragments")

    # Optionally write to SDF file
    if output_fragments:
        output_dir = Path(config.get("output", {}).get("directory", "."))
        output_dir.mkdir(parents=True, exist_ok=True)
        write_sdfs = config.get("focused_fragment_library", {}).get("output_sdf_per_fragment", True)
        output_sdf_dir = output_dir / "matched_molecules"
        if write_sdfs:
            output_sdf_dir.mkdir(parents=True, exist_ok=True)

        for group_name, group_df in grouped:
            frag_counter = Counter()
            total_mols = 0
            frag_to_mols = defaultdict(list)  # Track mols per fragment

            for smi in group_df["smiles"].dropna().unique():
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                total_mols += 1
                fragments = get_fragments(mol, method=fragment_type)

                unique_fragments = set(fragments)
                frag_counter.update(unique_fragments)
                for frag in unique_fragments:
                    frag_to_mols[frag].append(mol)

            if total_mols == 0:
                logger.warning(f"No valid molecules for group '{group_name}'")
                continue

            frequent_frags = {
                frag: count / total_mols
                for frag, count in frag_counter.items()
                if (count / total_mols) >= frequency_threshold
            }

            sorted_frags = sorted(frequent_frags.items(), key=lambda x: x[1], reverse=True)
            if max_fragments:
                sorted_frags = sorted_frags[:max_fragments]

            for frag, freq in sorted_frags:
                results.append({
                    "group": group_name,
                    "fragment_smiles": frag,
                    "frequency": round(freq, 4),
                    "count": frag_counter[frag],
                    "total_molecules": total_mols
                })

                if write_sdfs:
                    frag_safe = frag.replace("/", "_").replace("*", "x").replace("\\", "_")
                    sdf_path = output_sdf_dir / f"{group_name}__{frag_safe}.sdf"
                    writer = SDWriter(str(sdf_path))
                    for mol in frag_to_mols[frag]:
                        if mol:
                            mol.SetProp("Fragment", frag)
                            mol.SetProp("Group", str(group_name))
                            writer.write(mol)
                    writer.close()
        # Write csv
        out_file = output_dir / config.get("output", {}).get("filename", "focused_library.csv")
        result_df.to_csv(out_file, index=False)
        logger.info(f"Fragments written to: {out_file}")

    return (input_path.stem, result_df)


def get_fragments(mol, method="BRICS", merge_brics=False):
    """
    Return a list of fragment SMILES from the molecule using the selected method.
    """
    if method == "BRICS":
        try:
            frags = BRICS.BRICSDecompose(mol)
            if merge_brics:
                merged = Chem.FragmentOnBRICSBonds(mol)
                smiles = Chem.MolToSmiles(merged)
                return [smiles]
            return list(frags)
        except Exception:
            return []

    elif method == "RECAP":
        try:
            tree = RecapDecompose(mol)
            return list(tree.GetLeaves().keys()) if tree else []
        except Exception:
            return []

    elif method == "MURCKO":
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            return [Chem.MolToSmiles(scaffold)] if scaffold else []
        except Exception:
            return []

    elif method == "MACCS":
        keys = MACCSkeys.GenMACCSKeys(mol)
        return [f"MACCS_{i}" for i, bit in enumerate(keys) if bit]
    else:
        raise ValueError(f"Unsupported fragment type: {method}")
