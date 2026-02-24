from PDBGraphStore import PDBGraphStore
import Builder
from concurrent.futures import ProcessPoolExecutor, as_completed


def split_graph_store(pdb_store=PDBGraphStore, pdb_code_list=[]) -> tuple[PDBGraphStore, PDBGraphStore]:
    pdb_graphs = []
    for pdb_code in pdb_code_list:
        pdb_graphs.append(pdb_store.extract_pdb(pdb_code))
    
    pdb_store.remove_multi_pdb(pdb_code_list)

    pdb_store_2 = PDBGraphStore()

    for pdb_graph in pdb_graphs:
        pdb_store_2.insert_pdb(pdb_graph)

    return (pdb_store, pdb_store_2)

def merge_graph_stores(graph_stores=[PDBGraphStore]) -> PDBGraphStore:
    main_graph_store = graph_stores.pop()
    

    for graph_store in graph_stores:
        pdb_codes_to_insert = [pdb_code for pdb_code in graph_store.get_this_pdb_list() if pdb_code not in main_graph_store.get_this_pdb_list()]
        graphs_to_insert_list = [graph_store.extract(pdb_code) for pdb_code in pdb_codes_to_insert]
        [main_graph_store.insert_pdb(graph) for graph in graphs_to_insert_list]

        return main_graph_store

def extract_pdb_graphs_multiprocessing(pdb_store, pdb_codes=[], num_cpus=4) -> list:
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(pdb_store.extract_pdb, pdb_code) for pdb_code in pdb_codes]

        extracted_graphs = []

        for future in as_completed(futures):
            extracted_graphs.append(future.result())

        return extracted_graphs