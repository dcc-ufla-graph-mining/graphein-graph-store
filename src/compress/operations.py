from PDBGraphStore import PDBGraphStore
import Builder
import concurrent.futures as futures

def split_graph_store(graph_store=PDBGraphStore, pdb_code_list=[]):
    for _ in range(len(pdb_code_list)):
        pdbs_to_insert = graph_store.get_pdb_list()
        pdbs_to_insert = list(filter(lambda x: x in pdb_code_list, pdbs_to_insert))
        edge_funcs = ["aromatic", "bb_carbonyl_carbonyl"]

        pdb_graph_list = graph_store.extract_pdb_graphs(pdbs_to_insert, edge_funcs)
        pdb_graphs_to_insert_list = {}

        for g in pdb_graph_list:
            pdb_graphs_to_insert_list.setdefault(g.graph["pdb_code"], [])
            pdb_graphs_to_insert_list[g.graph["pdb_code"]].append(g)

        graph_store2 = PDBGraphStore()

        print(f'pdbs to insert: {pdbs_to_insert}')
        graph_store2.insert_pdbs(pdb_graphs_to_insert_list)
        print(graph_store.remove_multiple_pdbs(pdbs_to_insert))

    return graph_store, graph_store2

def merge_graph_stores(graph_stores=[PDBGraphStore]):
    if len(graph_stores) < 2:
        return
    
    main_graph_store = graph_stores.pop(0)
    
    for _ in range(len(graph_stores)):
        aux_graph_store = graph_stores.pop(0)
        pdbs_to_insert = aux_graph_store.get_pdb_list()
        edge_funcs = ["aromatic", "bb_carbonyl_carbonyl"]

        pdb_graph_list = aux_graph_store.extract_pdb_graphs(pdbs_to_insert, edge_funcs)
        pdb_graphs_to_insert_list = {}
        for g in pdb_graph_list:
            pdb_graphs_to_insert_list.setdefault(g.graph["pdb_code"], [])
            pdb_graphs_to_insert_list[g.graph["pdb_code"]].append(g)
        
        main_graph_store.insert_pdbs(graphs=pdb_graphs_to_insert_list)

        del aux_graph_store


def extract_pdb_graphs_multiprocessing(pdb_store, pdb_codes=[], edge_constructions_functions=[], num_cpus=4):

    results = []

    with futures.ProcessPoolExecutor(max_workers=num_cpus) as executor:
        future_to_pdb_code = {
            executor.submit(pdb_store.extract_pdb_graphs, [pdb_code], edge_constructions_functions): pdb_code 
            for pdb_code in pdb_codes
        }
        
        for i, future in enumerate(futures.as_completed(future_to_pdb_code)):
            pdb_code = future_to_pdb_code[future]
            print(f"Processando {i+1}/{len(pdb_codes)}: {pdb_code}")
            try:
                data = future.result()
                results.append(data)
                print(f'✓ Sucesso: {pdb_code}')
            except Exception as e:
                print(f'✗ Erro em {pdb_code}: {type(e).__name__}: {e}')
                import traceback
                traceback.print_exc()
        
    return results