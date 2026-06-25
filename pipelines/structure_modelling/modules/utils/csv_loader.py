import csv

def load_sequences(csv_path):
    entries = []

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)

        for row in reader:
            entry_id = row.get("id")

            if not entry_id:
                raise ValueError("CSV row missing 'id'")

            proteins = {}
            ligands = []

            # --- proteins ---
            if row.get("chain_A"):
                proteins["A"] = row["chain_A"].strip()

            if row.get("chain_B"):
                proteins["B"] = row["chain_B"].strip()

            # --- ligands ---
            smiles_field = row.get("ligand_SMILES")

            if smiles_field:
                ligands = [
                    s.strip()
                    for s in smiles_field.split(";")
                    if s.strip()
                ]

            # ✅ strict validation
            if not proteins and not ligands:
                raise ValueError(
                    f"Entry '{entry_id}' has no proteins or ligands"
                )

            entries.append({
                "id": entry_id,
                "proteins": proteins,
                "ligands": ligands
            })

    return entries