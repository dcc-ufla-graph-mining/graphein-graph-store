import pandas as pd
import wget
import re

datasets = {
    "ligand": {
        "SO4": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/SO4.txt",
        "PO4": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/PO4.txt",
        "NAG": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/NAG.txt",
        "HEM": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/HEM.txt",
        "BME": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/BME.txt",
        "EDO": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/EDO.txt",
        "PLP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/PLP.txt",
    },

    "nucleotides": {
        "ATP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/ATP.txt",
        "ADP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/ADP.txt",
        "GTP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/GTP.txt",
        "GDP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/GDP.txt",
        "NAD": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/NAD.txt",
        "FAD": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/FAD.txt",
        "FMN": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/FMN.txt",
        "UDP": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/UDP.txt",  
    },

    "metals": {
        "FE": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/FE.txt",
        "MG": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/MG.txt",
        "CA": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/CA.txt",
        "MN": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/MN.txt",
        "ZN": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/ZN.txt",
        "CO": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/CO.txt",
        "NI": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/NI.txt",
    },

    "nicleic": {
        "rna": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/rna.txt",
        "dna": "https://webs.iiitd.edu.in/raghava/ccpdb/datasets/dna.txt",
    },
}



import re
import pandas as pd
from Bio import SeqIO
import os

def parse_dset(filename: str, database) -> pd.DataFrame:
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

    os.remove(f'{filename}.txt')

    df = pd.DataFrame({
        "PDB": identifiers,
        "chain": chains,
        "sequence": sequences,
        "interacting_residues": interactions,
        "length": lengths,
        "interactor": filename.split(".")[0]
    })

    this_pdbs = df["PDB"].tolist()

    this_pdbs = sorted(set(this_pdbs))

    with open(f"data/{database}_{filename}.txt", "a") as f:
        for pdb in this_pdbs:
            f.write(pdb + "\n")
    
    return this_pdbs

def main():
    for dataset, sub in datasets.items():
        dat = list()
        for k, v in sub.items():
            print(f"Downloading: {k}")
            wget.download(v)
            dat.extend(parse_dset(k, dataset))

        unique_pdbs = sorted(set(dat))

        with open(f"data/{dataset}.txt", "w") as f:
            for pdb in unique_pdbs:
                f.write(pdb + "\n")


if __name__ == "__main__":
    main()