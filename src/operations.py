from PDBGraphStore import PDBGraphStore
from concurrent.futures import ProcessPoolExecutor, as_completed

def remove_graph_from_store(pdbs_to_remove: list, pdb_store: PDBGraphStore) -> PDBGraphStore:
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
    pdb_store_code_list = list(pdb_store.get_this_pdb_list())

    if not pdb_code_list:
        mid = len(pdb_store_code_list) // 2

        list_1 = pdb_store_code_list[:mid]
        list_2 = pdb_store_code_list[mid:]
    else:
        list_1 = [x for x in pdb_store_code_list if x in pdb_code_list]
        list_2 = [x for x in pdb_store_code_list if x not in pdb_code_list]

    store_1 = remove_graph_from_store(list_2, pdb_store)
    store_2 = remove_graph_from_store(list_1, pdb_store)

    return store_1, store_2

def merge_graph_stores(graph_stores: list) -> PDBGraphStore:
    main_graph_store = graph_stores.pop()
    main_pdbs = set(main_graph_store.get_this_pdb_list())

    for graph_store in graph_stores:
        for pdb_code in set(graph_store.get_this_pdb_list()):
            if pdb_code not in main_pdbs:
                graph = graph_store.extract_pdb(pdb_code)
                main_graph_store.insert_pdb({pdb_code: [graph]})
                main_pdbs.add(pdb_code)

    return main_graph_store

def extract_pdb_graphs_multiprocessing(pdb_store: PDBGraphStore, pdb_codes: list, num_cpus=4) -> list:
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(pdb_store.extract_pdb, pdb_code) for pdb_code in pdb_codes]

        extracted_graphs = []

        for future in as_completed(futures):
            extracted_graphs.append(future.result())

    return extracted_graphs