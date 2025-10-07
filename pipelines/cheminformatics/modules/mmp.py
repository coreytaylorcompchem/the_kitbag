import pandas as pd
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from pipeline.task_registry import register_task
from pipeline.logger import setup_logger

logger = setup_logger(__name__, debug_mode=False, simple_format=True)

def write_mmpdb_inputs(df, output_prefix, smiles_col="smiles", id_col="id", props_cols=None):
    output_prefix = Path(output_prefix)

    if id_col not in df.columns:
        df[id_col] = df.index.astype(str)

    smi_path = output_prefix.with_suffix(".smi")
    with open(smi_path, "w") as f:
        # f.write("smiles\tid\n")  # header required for mmpdb
        for _, row in df.iterrows():
            smi = row.get(smiles_col)
            mol_id = row.get(id_col)
            if pd.notna(smi) and pd.notna(mol_id):
                f.write(f"{smi}\t{mol_id}\n")
    logger.info(f"[✔] SMILES file written: {smi_path}")

    prop_path = None
    if props_cols:
        missing_cols = [col for col in props_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing property columns in DataFrame: {missing_cols}")

        props_df = df[[id_col] + props_cols].copy()
        props_df.rename(columns={id_col: "id"}, inplace=True)
        prop_path = output_prefix.with_suffix(".props")
        props_df.to_csv(prop_path, index=False, sep="\t", encoding="utf-8", lineterminator="\n")
        logger.info(f"[✔] Properties file written: {prop_path}")
    else:
        logger.info("[ℹ] No props_cols provided — skipping props file.")

    return smi_path, prop_path

def run_transform(smiles, mmpdb_file):
    cmd = [
        "mmpdb", "transform",
        str(mmpdb_file),
        "--smiles", smiles,
        "--property", "property"
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        lines = result.stdout.strip().splitlines()
        # Remove header line if present
        if lines and lines[0].startswith("smiles"):
            lines = lines[1:]
        return "\n".join(lines)
    except subprocess.CalledProcessError as e:
        logger.warning(f"[!] Transform command failed for SMILES {smiles}")
        logger.debug(e.stderr)
        return None

@register_task("mmp_analysis", category="Analysis", description="Matched Molecular Pair (MMP) analysis via mmpdb")
def mmp_analysis(config, data=None):
    input_file = config.get("input_file")
    activity_col = config.get("activity_col", "pActivity")
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out_filename = config.get("output", {}).get("filename", Path(input_file).stem + "_mmp.csv")

    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", activity_col])

    # Rename activity column to "property" as required by mmpdb
    df = df.rename(columns={activity_col: "property"})

    smi_path, prop_csv_path = write_mmpdb_inputs(
        df,
        output_prefix=output_dir / "temp_input",
        smiles_col="smiles",
        id_col="id",
        props_cols=["property"]
    )

    fragdb = output_dir / "temp.fragdb"
    mmpdb_file = output_dir / "temp.mmpdb"
    transform_out = output_dir / out_filename

    cmd1 = ["mmpdb", "fragment", str(smi_path), "-o", str(fragdb)]
    cmd2 = ["mmpdb", "index", str(fragdb), "--properties", str(prop_csv_path), "-o", str(mmpdb_file)]

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd1)}")
    subprocess.run(cmd1, check=True)

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd2)}")
    subprocess.run(cmd2, check=True)

    smiles_list = df["smiles"].tolist()

    results = []
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(run_transform, smi, mmpdb_file): smi for smi in smiles_list}
        for future in as_completed(futures):
            res = future.result()
            if res:
                results.append(res)
            else:
                logger.warning(f"[!] Transform failed for SMILES: {futures[future]}")

    # Write combined output with one header and all data rows
    with open(transform_out, "w") as fout:
        fout.write("\n".join(results))
        fout.write("\n")

    result_df = pd.read_csv(transform_out)
    logger.info(f"MMP: got {len(result_df)} transform rows")

    return (Path(input_file).stem, result_df)
