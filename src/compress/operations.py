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
    pass