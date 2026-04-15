from PDBGraphStore import PDBGraphStore
from concurrent.futures import ProcessPoolExecutor, as_completed

def remove_graph_from_graph_store(pdbs_to_remove: list, pdb_store: PDBGraphStore):
    pdbs_to_remove_set = set(pdbs_to_remove)

    all_pdbs = pdb_store.get_this_pdb_list()
    pdbs_to_keep = [x for x in all_pdbs if x not in pdbs_to_remove_set]

    graphs_to_insert = {}

    for pdb_code in pdbs_to_keep:
        graphs_to_insert[pdb_code] = [pdb_store.extract_pdb(pdb_code)]

    new_store = PDBGraphStore(None)
    
    new_store.insert_pdb(graphs_to_insert)

    return new_store

def split_graph_store(pdb_store: PDBGraphStore, pdb_code_list: list) -> tuple:
    pdb_graphs = []
    for pdb_code in pdb_code_list:
        pdb_graphs.append(pdb_store.extract_pdb(pdb_code))
    
    pdb_store.remove_multi_pdb(pdb_code_list)

    pdb_store_2 = PDBGraphStore()

    for pdb_graph in pdb_graphs:
        pdb_store_2.insert_pdb(pdb_graph)

    return (pdb_store, pdb_store_2)

def merge_graph_stores(graph_stores: list) -> PDBGraphStore:
    main_graph_store = graph_stores.pop()

    for graph_store in graph_stores:
        pdb_codes_to_insert = [pdb_code for pdb_code in graph_store.get_this_pdb_list() if pdb_code not in main_graph_store.get_this_pdb_list()]
        graphs_to_insert_list = [graph_store.extract(pdb_code) for pdb_code in pdb_codes_to_insert]
        [main_graph_store.insert_pdb(graph) for graph in graphs_to_insert_list]

        return main_graph_store

def extract_pdb_graphs_multiprocessing(pdb_store: PDBGraphStore, pdb_codes: list, num_cpus=4) -> list:
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(pdb_store.extract_pdb, pdb_code) for pdb_code in pdb_codes]

        extracted_graphs = []

        for future in as_completed(futures):
            extracted_graphs.append(future.result())

    return extracted_graphs