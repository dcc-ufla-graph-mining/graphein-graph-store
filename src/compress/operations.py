from PDBGraphStore import PDBGraphStore
import Builder

def split_graph_store(graph_store=PDBGraphStore, pdb_code_list=[]):
    for pdb_code in pdb_code_list:
        graph_store.remove_pdb(pdb_code)
    
    node_to_id,\
    edge_to_id,\
    pdb_to_nodes,\
    pdb_to_edges,\
    node_attrs,\
    edge_attrs,\
    node_attr_keys,\
    edge_attr_keys = Builder.compress_pdb_graphs()

    graph_store2 = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, node_attr_keys, edge_attr_keys)

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
