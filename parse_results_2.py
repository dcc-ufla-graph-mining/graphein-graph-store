# Creating a parser that preserves separate "memory size" values for each insertion
import os
import re
from pathlib import Path
import pandas as pd

SAMPLE_LOG = """Result log for ligand_PLP
number of nodes: 293709
number of edges initially (with one func): 109551        
memory size initially: 265.63 MB
Number of edges after 1 insertion (aromatic_sulphur): 109551            
Memory size: 265.63 MB
Number of edges after 2 insertion (bb_carbonyl_carbonyl): 147964            
Memory size: 295.91 MB
Number of edges after 3 insertion (cation_pi): 147964            
Memory size: 295.91 MB
Number of edges after 4 insertion (delaunay): 3675535            
Memory size: 3333.72 MB
Number of edges after 5 insertion (disulfide): 3675535            
Memory size: 3333.72 MB
Number of edges after 6 insertion (hbond): 3675594            
Memory size: 3333.77 MB
Number of edges after 7 insertion (hydrophobic): 3827322            
Memory size: 3444.65 MB
Number of edges after 8 insertion (ionic): 3927645            
Memory size: 3514.41 MB
Number of edges after 9 insertion (knn): 3977176            
Memory size: 3553.12 MB
Number of edges after 10 insertion (peptide_bond): 3977204            
Memory size: 3553.14 MB
"""

def sanitize_key(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r'[^\w]+', '_', s)  # replace non-alphanumeric with underscore
    s = re.sub(r'__+', '_', s)
    s = s.strip('_')
    return s

def unique_key(data: dict, key: str) -> str:
    """If key exists, append suffix _2, _3... to make it unique."""
    if key not in data:
        return key
    i = 2
    while f"{key}_{i}" in data:
        i += 1
    return f"{key}_{i}"

def parse_log_text(text: str) -> dict:
    data = {}
    last_insertion = None  # will be tuple (index, name_sanitized) or "initial"
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        # Result log for ...
        m = re.match(r'(?i)^result log for\s+(.+)', line)
        if m:
            data['result_log_for'] = m.group(1).strip()
            continue

        # number of nodes
        m = re.match(r'(?i)^number of nodes[:\s]+([\d,]+)', line)
        if m:
            data['number_of_nodes'] = int(m.group(1).replace(',', ''))
            continue

        # number of edges initially
        m = re.match(r'(?i)^number of edges initially[^\:]*:\s*([\d,]+)', line)
        if m:
            data['number_of_edges_initially'] = int(m.group(1).replace(',', ''))
            last_insertion = 'initial'
            continue

        # memory size initially
        m = re.match(r'(?i)^memory size initially[:\s]+([\d\.]+)\s*mb', line)
        if m:
            data['memory_size_initial'] = float(m.group(1))
            last_insertion = None
            continue

        # Number of edges after N insertion (optional name)
        m = re.match(r'(?i)^number of edges after\s+(\d+)\s+insertion(?:\s*\((.*?)\))?\s*:\s*([\d,]+)', line)
        if m:
            idx = int(m.group(1))
            name = m.group(2) or ''
            name_s = sanitize_key(name) if name else ''
            key_edges = f"number_of_edges_after_{idx}" + (f"_{name_s}" if name_s else "")
            data[key_edges] = int(m.group(3).replace(',', ''))
            last_insertion = (idx, name_s)  # set context so next Memory size lines attach here
            continue

        # Memory size (general) - attach to last_insertion if present, else store generically
        m = re.match(r'(?i)^memory size[:\s]+([\d\.]+)\s*mb', line)
        if m:
            value = float(m.group(1))
            if last_insertion == 'initial':
                key = 'memory_size_initial'
            elif isinstance(last_insertion, tuple):
                idx, name_s = last_insertion
                key = f"memory_size_after_{idx}" + (f"_{name_s}" if name_s else "")
            else:
                # generic memory size (no context) -> make unique
                key = unique_key(data, 'memory_size')
            data[key] = value
            # keep last_insertion as is (there may be multiple lines referencing same insertion)
            continue

        # fallback - generic key: value pattern
        m = re.match(r'(.+?):\s*(.+)', line)
        if m:
            k = sanitize_key(m.group(1))
            v = m.group(2).strip()
            # try convert to int/float if possible (strip MB and commas)
            if v.lower().endswith('mb'):
                try:
                    v_num = float(v[:-2].strip())
                    val = v_num
                except Exception:
                    val = v
            else:
                v_clean = v.replace(',', '')
                try:
                    if '.' in v_clean:
                        val = float(v_clean)
                    else:
                        val = int(v_clean)
                except Exception:
                    val = v
            k_unique = unique_key(data, k)
            data[k_unique] = val
            continue

    return data

def parse_results_dir(results_dir: str = "./results"):
    p = Path(results_dir)
    if not p.exists():
        p.mkdir(parents=True, exist_ok=True)
        # write sample file so you can see the parser output immediately
        (p / "ligand_PLP.txt").write_text(SAMPLE_LOG, encoding="utf-8")

    parsed_list = []
    for filepath in sorted(p.iterdir()):
        if not filepath.is_file():
            continue
        text = filepath.read_text(encoding="utf-8")
        parsed = parse_log_text(text)
        parsed['filename'] = filepath.name
        parsed_list.append(parsed)

    df = pd.DataFrame(parsed_list)
    # put filename first column if present
    cols = list(df.columns)
    if 'filename' in cols:
        cols = ['filename'] + [c for c in cols if c != 'filename']
        df = df[cols]

    # save CSV to /mnt/data so you can download
    out_path = Path("./results.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved CSV to: {out_path}")

# Run the parser on ./results (creates a sample file if none exist)
parse_results_dir("./results")

