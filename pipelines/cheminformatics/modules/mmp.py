import csv

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

def run_transform(mmpdb_file, smiles, property_name):
    cmd = [
        "mmpdb", "transform",
        str(mmpdb_file),
        "--smiles", smiles,
        "--property", property_name
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        output_lines = result.stdout.strip().splitlines()
        # Add property name to each line to identify property in merged output
        # But mmpdb output usually includes headers already, so careful to only add to data lines
        # Let's just return raw output and handle header outside
        return property_name, smiles, output_lines
    except subprocess.CalledProcessError as e:
        logger.warning(f"[!] Transform failed for SMILES {smiles} with property {property_name}")
        logger.debug(e.stderr)
        return property_name, smiles, None

@register_task("mmp_analysis", category="Analysis", description="Matched Molecular Pair (MMP) analysis via mmpdb")
def mmp_analysis(config, data=None):
    import csv
    import subprocess
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import pandas as pd
    from pathlib import Path
    import logging

    logger = logging.getLogger(__name__)

    input_file = config.get("input_file")
    activity_col = config.get("activity_col", "pActivity")
    output_dir = Path(config.get("output", {}).get("directory", "outputs/mmp"))
    output_dir.mkdir(parents=True, exist_ok=True)
    out_filename = config.get("output", {}).get("filename", Path(input_file).stem + "_mmp.tsv")

    df = pd.read_csv(input_file)
    df = df.dropna(subset=["smiles", activity_col])

    # Recode toxic and reactive flags to numeric 0/1 if present
    if 'toxic_flag' in df.columns:
        df['toxic_flag'] = df['toxic_flag'].map({'N': 0, 'Y': 1})
    if 'reactive_flag' in df.columns:
        df['reactive_flag'] = df['reactive_flag'].map({'N': 0, 'Y': 1})

    # Rename activity column to "property" for mmpdb as required
    df = df.rename(columns={activity_col: "property"})

    props_cols = ["property", "mw", "logp", "hbd", "hba", "rotatable_bonds", "tpsa", "qed", "stereocenters"]
    if 'toxic_flag' in df.columns:
        props_cols.append("toxic_flag")
    if 'reactive_flag' in df.columns:
        props_cols.append("reactive_flag")

    smi_path, prop_csv_path = write_mmpdb_inputs(
        df,
        output_prefix=output_dir / "temp_input",
        smiles_col="smiles",
        id_col="id",
        props_cols=props_cols
    )

    fragdb = output_dir / "temp.fragdb"
    mmpdb_file = output_dir / "temp.mmpdb"
    transform_out = output_dir / out_filename

    # Run fragment and index steps
    cmd1 = ["mmpdb", "fragment", str(smi_path), "-o", str(fragdb)]
    cmd2 = ["mmpdb", "index", str(fragdb), "--properties", str(prop_csv_path), "-o", str(mmpdb_file)]

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd1)}")
    subprocess.run(cmd1, check=True)

    logger.info(f"[mmp_analysis] Running: {' '.join(cmd2)}")
    subprocess.run(cmd2, check=True)

    smiles_list = df["smiles"].tolist()

    all_results = []
    header_line = None

    with ProcessPoolExecutor() as executor:
        futures = []
        for prop in props_cols:
            for smi in smiles_list:
                futures.append(executor.submit(run_transform, mmpdb_file, smi, prop))

        for future in as_completed(futures):
            prop, smi, output_lines = future.result()
            if output_lines:
                if header_line is None:
                    header_line = output_lines[0]

                for line in output_lines[1:]:
                    if line.strip():
                        all_results.append({
                            "property": prop,
                            "smiles": smi,
                            "data": line
                        })
            else:
                logger.warning(f"[!] No output for SMILES {smi} with property {prop}")

    if header_line is None:
        raise RuntimeError("No MMP transform results generated.")

    columns = ["property", "smiles"] + header_line.strip().split("\t")

    rows = []
    for res in all_results:
        row_fields = res["data"].strip().split("\t")
        row = [res["property"], res["smiles"]] + row_fields
        rows.append(row)

    with open(transform_out, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)

    result_df = pd.read_csv(transform_out, sep="\t")
    logger.info(f"MMP: got {len(result_df)} transform rows")

    return (Path(input_file).stem, result_df)


