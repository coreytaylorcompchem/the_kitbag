import csv

def parse_residue_list(value):
    if value is None:
        return []

    value = str(value).strip()

    if not value:
        return []

    residues = []

    for part in value.replace(",", ";").split(";"):
        part = part.strip()

        if not part:
            continue

        residues.append(int(part))

    return residues


def parse_csv_constraints(row):
    chain1 = row.get("constraint_chain1")
    chain2 = row.get("constraint_chain2")
    residues1 = parse_residue_list(row.get("constraint_residues1"))
    residues2 = parse_residue_list(row.get("constraint_residues2"))

    if not chain1 or not chain2 or not residues1 or not residues2:
        return []

    return [
        {
            "type": "chain_contact",
            "chain1": str(chain1).strip(),
            "residues1": residues1,
            "chain2": str(chain2).strip(),
            "residues2": residues2,
        }
    ]

def load_sequences(csv_path):
    entries = []

    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)

        if not reader.fieldnames:
            raise ValueError(f"CSV has no header: {csv_path}")

        fieldnames = [
            field.strip()
            for field in reader.fieldnames
        ]

        reader.fieldnames = fieldnames

        if "id" not in reader.fieldnames:
            raise ValueError(
                f"CSV missing required 'id' column. "
                f"Detected columns: {reader.fieldnames}"
            )

        for row_idx, row in enumerate(reader, start=2):
            entry_id = row.get("id")

            if entry_id:
                entry_id = entry_id.strip()

            if not entry_id:
                raise ValueError(
                    f"CSV row {row_idx} missing 'id'. "
                    f"Detected row keys: {list(row.keys())}. "
                    f"Row contents: {row}"
                )

            proteins = {}
            ligands = []
            templates = []

            # # --- proteins ---
            # if row.get("chain_A"):
            #     proteins["A"] = row["chain_A"].strip()

            # if row.get("chain_B"):
            #     proteins["B"] = row["chain_B"].strip()

            for key, value in row.items():
                if not key.startswith("chain_"):
                    continue

                if value is None:
                    continue

                value = value.strip()

                if not value:
                    continue

                chain_id = key.replace("chain_", "")
                proteins[chain_id] = value

            # --- ligands ---
            smiles_field = row.get("ligand_SMILES")
            
            if smiles_field:
                ligands = [
                    s.strip()
                    for s in smiles_field.split(";")
                    if s.strip()
                ]

            template_field = row.get("template_pdbs")
            
            if template_field:
                templates = [
                    t.strip()
                    for t in template_field.split(";")
                    if t.strip()
                ]

            # strict validation
            if not proteins and not ligands:
                raise ValueError(
                    f"Entry '{entry_id}' has no proteins or ligands"
                )

            constraints = parse_csv_constraints(row)

            if constraints:
                print(
                    f"[CSV] {entry_id}: loaded {len(constraints)} constraint(s)"
                )

            entries.append({
                "id": entry_id,
                "proteins": proteins,
                "ligands": ligands,
                "templates": templates,
                "constraints": constraints,
                "msas": {}
            })


    return entries