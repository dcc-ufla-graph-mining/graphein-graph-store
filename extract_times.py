import os
import re
import csv

def extract_last_increment(filepath):
    """Extrai o valor do campo 'Increment' na última linha do arquivo."""
    increment_value = None
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            match = re.search(r'\b([\d.]+)\s*MiB', line)
            if match:
                increment_value = float(match.group(1))
    return increment_value

def main():
    times_dir = './times'
    output_csv = 'memory_summary.csv'

    # Coleta todos os arquivos *_v1.txt e *_v2.txt
    v1_files = [f for f in os.listdir(times_dir) if f.endswith('_v1.txt')]
    v2_files = [f for f in os.listdir(times_dir) if f.endswith('_v2.txt')]

    datasets = {}
    
    # Lê valores v1
    for f in v1_files:
        dataset = f.replace('_v1.txt', '')
        path = os.path.join(times_dir, f)
        datasets[dataset] = {'v1': extract_last_increment(path), 'v2': None}
    
    # Lê valores v2
    for f in v2_files:
        dataset = f.replace('_v2.txt', '')
        path = os.path.join(times_dir, f)
        if dataset not in datasets:
            datasets[dataset] = {'v1': None, 'v2': None}
        datasets[dataset]['v2'] = extract_last_increment(path)

    # Escreve CSV
    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['dataset', 'v1', 'v2', 'compression_ratio'])
        for dataset, values in sorted(datasets.items()):
            writer.writerow([dataset, values['v1'], values['v2'],(float(values['v1'])/float(values['v2']))])

    print(f'Arquivo gerado: {output_csv}')

if __name__ == '__main__':
    main()
