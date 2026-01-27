from PDBGraphStore import PDBGraphStore
import Builder
from concurrent.futures import ProcessPoolExecutor, as_completed


#TODO
def split_graph_store(graph_store=PDBGraphStore, pdb_code_list=[]) -> tuple[PDBGraphStore, PDBGraphStore]:
    pass

#TODO
def merge_graph_stores(graph_stores=[PDBGraphStore]) -> PDBGraphStore:
    pass

def extract_pdb_graphs_multiprocessing(pdb_store, pdb_codes=[], num_cpus=4) -> list:
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(pdb_store.extract_pdb, pdb_code) for pdb_code in pdb_codes]

        extracted_graphs = []

        for future in as_completed(futures):
            extracted_graphs.append(future.result())
