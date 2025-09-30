import os
import csv

RESULTS_DIR = "./results"
DATA_DIR = "./data"
OUTPUT_CSV = "results.csv"

def parse_log_file(filepath):
    """Lê um arquivo de log e retorna um dicionário {chave: valor}"""
    data = {}
    with open(filepath, "r") as f:
        for line in f:
            if ":" in line:
                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()
                try:
                    value = float(value)
                except ValueError:
                    pass  # mantém string se não for número
                data[key] = value
    return data

def count_codes(filepath):
    """Conta quantos códigos (linhas não vazias) existem no arquivo .txt"""
    with open(filepath, "r") as f:
        return sum(1 for line in f if line.strip())

def main():
    all_data = {}
    keys = []

    # Itera sobre os arquivos no diretório ./results
    for filename in os.listdir(RESULTS_DIR):
        if filename.endswith("_results.log"):
            dataset = filename.replace("_results.log", "")
            filepath = os.path.join(RESULTS_DIR, filename)

            parsed = parse_log_file(filepath)

            # conta códigos do dataset correspondente
            txt_path = os.path.join(DATA_DIR, dataset + ".txt")
            if os.path.exists(txt_path):
                code_count = count_codes(txt_path)
            else:
                code_count = 0  # se não existir .txt correspondente

            parsed = {"code_count": code_count, **parsed}
            all_data[dataset] = parsed

            # mantém ordem das chaves na primeira vez
            if not keys:
                keys = list(parsed.keys())
            else:
                for k in parsed.keys():
                    if k not in keys:
                        keys.append(k)

    # Escreve CSV
    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["dataset"] + keys
        writer.writerow(header)
        for dataset, values in all_data.items():
            row = [dataset] + [values.get(k, "") for k in keys]
            writer.writerow(row)

    print(f"CSV gerado em: {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
