import os
import csv
import re

# Pasta onde estão os arquivos .txt
input_folder = "results"
# Nome do arquivo de saída CSV
output_csv = "results.csv"

# Expressão regular para capturar os pares chave: valor
pattern = re.compile(r"^(.*?):\s+([0-9eE\.\-]+)$")

# Lista para armazenar os dados extraídos
data_rows = []
headers_set = set()

# Lê todos os arquivos .txt da pasta
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        filepath = os.path.join(input_folder, filename)
        with open(filepath, "r") as f:
            lines = f.readlines()

        row_data = {"filename": filename.split("_results")[0]}
        
        for line in lines:
            match = pattern.match(line.strip())
            if match:
                key = match.group(1).strip()
                value = float(match.group(2))
                row_data[key] = value
                headers_set.add(key)

        data_rows.append(row_data)

# Ordenar headers para consistência (exceto filename, que vai primeiro)
headers = ["filename"] + sorted(headers_set)

# Escreve o arquivo CSV
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=headers)
    writer.writeheader()
    for row in data_rows:
        writer.writerow(row)

print(f"Arquivo CSV '{output_csv}' gerado com sucesso.")
