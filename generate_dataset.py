import pandas as pd
import wget
import re

DSET_URL = {
    "SO4": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/SO4.txt",
    "PO4": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/PO4.txt",
    "NAG": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/NAG.txt",
    "HEM": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/HEM.txt",
    "BME": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/BME.txt",
    "EDO": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/EDO.txt",
    "PLP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/PLP.txt",
}

import re
import pandas as pd
from Bio import SeqIO

all_pdbs = []

def parse_dset(filename: str) -> pd.DataFrame:
    with open(f'{filename}.txt') as fasta_file:
        identifiers = []
        chains = []
        sequences = []
        interactions = []
        lengths = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            identifiers.append(seq_record.id[:-1])
            chains.append(seq_record.id[-1])
            lengths.append(len(seq_record.seq))

            parsed_sequence_and_interactions = re.split(
                ';',
                re.sub(r"\+|-", lambda match: ';' + match.group(), str(seq_record.seq), count=1),
                maxsplit=1
            )
            sequences.append(parsed_sequence_and_interactions[0])
            interactions.append(parsed_sequence_and_interactions[1])

    df = pd.DataFrame({
        "PDB": identifiers,
        "chain": chains,
        "sequence": sequences,
        "interacting_residues": interactions,
        "length": lengths,
        "interactor": filename.split(".")[0]
    })

    # Adiciona os PDBs ao conjunto global
    all_pdbs.extend(df["PDB"].tolist())
    
    return df

# Exemplo de uso com múltiplos arquivos
for k, v in DSET_URL.items():
    print(f"Downloading: {k}")
    wget.download(v)
    parse_dset(k)

# Remover duplicatas e salvar todos os PDBs em um único arquivo
unique_pdbs = sorted(set(all_pdbs))

with open("all_pdbs.txt", "w") as f:
    for pdb in unique_pdbs:
        f.write(pdb + "\n")
