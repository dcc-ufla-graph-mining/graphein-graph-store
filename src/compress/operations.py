from PDBGraphStore import PDBGraphStore
import Builder
import concurrent.futures as futures

def split_graph_store(graph_store=PDBGraphStore, pdb_code_list=[]) -> tuple[PDBGraphStore, PDBGraphStore]:
    pass

def merge_graph_stores(graph_stores=[PDBGraphStore]) -> PDBGraphStore:
    pass

def extract_pdb_graphs_multiprocessing(pdb_store, pdb_codes=[], num_cpus=4) -> list:
    pass