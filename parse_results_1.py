import os
import csv

# Pasta onde estão os arquivos
INPUT_DIR = "results"
OUTPUT_CSV = "results_experimento_1.csv"

def parse_file(path):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                data[key] = value
    return data

def main():
    rows = []
    all_keys = set()

    for filename in os.listdir(INPUT_DIR):
        filepath = os.path.join(INPUT_DIR, filename)

        if not os.path.isfile(filepath):
            continue
        
        record = parse_file(filepath)
        record["dataset"] = filename

        rows.append(record)
        all_keys.update(record.keys())

    # Garantir que a coluna 'dataset' seja a primeira
    all_keys = ["dataset"] + sorted(k for k in all_keys if k != "dataset")

    # Escrever CSV
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(rows)

    print(f"CSV gerado com sucesso: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
