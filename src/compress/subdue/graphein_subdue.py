from graphein.protein.config import ProteinGraphConfig
from graphein.protein.edges.atomic import add_atomic_edges
from graphein.protein.graphs import construct_graph
from Subdue import nx_subdue

params_to_change = {"granularity": "atom", "edge_construction_functions": [add_atomic_edges]}

config = ProteinGraphConfig(**params_to_change)
print(config.dict())

# List of PDB codes
pdb_codes = ["1CRN", "4HHB", "2MNR"]  # Example PDB codes


# Function to construct graphs from PDB codes
protein_graphs = {}
for pdb_code in pdb_codes:
    graph = construct_graph(config=config, pdb_code=pdb_code)
    protein_graphs[pdb_code] = graph  # Store graph

print(protein_graphs)

g = protein_graphs["1CRN"]
print(g.nodes)
print(g.edges)
out = nx_subdue(g)

print(out)
