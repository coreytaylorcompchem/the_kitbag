import csv

def load_sequences(csv_path):
    entries = []

    with open(csv_path) as f:
        reader = csv.DictReader(f)

        for row in reader:
            seqs = {}
            for col, val in row.items():
                if col.startswith("chain_") and val:
                    chain_id = col.split("_")[1]
                    seqs[chain_id] = val

            entries.append({
                "id": row["id"],
                "sequences": seqs
            })

    return entries