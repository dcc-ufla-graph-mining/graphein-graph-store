import bisect
import pickle
import random
import time
import pandas as pd

import networkx as nx
import numpy as np
from bidict import bidict
from graphein.protein import add_atomic_edges
from pympler import asizeof
from pyroaring import BitMap, BitMap64


############################################################################################
#here happens the compression

# e = edge, v = vertex, g = graph
def compress_with_composition(protein_graphs):
    edge_to_pdbs = {}   
    node_to_pdbs = {}  
    edge_attrs = {}   
    node_attrs = {}    

    edge_attr_keys_list = list(list(next(iter(protein_graphs.values())).edges(data=True))[0][2].keys())
    node_attr_keys_list = list(list(next(iter(protein_graphs.values())).nodes(data=True))[0][1].keys())

    edge_attr_keys = {}
    node_attr_keys = {}  

    for key in edge_attr_keys_list: 
        edge_attr_keys[key] = []
    for key in node_attr_keys_list: 
        node_attr_keys[key] = []

    del edge_attr_keys_list
    del node_attr_keys_list
    
    pdb_to_view = {}

    for pdb_code, g in protein_graphs.items():
        pdb_to_view[pdb_code] = ([BitMap64()], [BitMap64()], [], [])  # nodes, edges, node attributes, edge attributes
        
        for node in g.nodes():
            if node not in node_to_pdbs:
                node_to_pdbs[node] = []
                node_attrs[node] = []
                
                for value in node_attr_keys:
                    attr_value = g.nodes[node][value]

                    if isinstance(attr_value, pd.Series):
                        attr_value = tuple([attr_value.tolist(), attr_value.name, attr_value.index])
                    elif isinstance(attr_value, np.ndarray):
                        attr_value = tuple([tuple(attr_value)])

                    if attr_value not in node_attr_keys[value]:
                        node_attr_keys[value].append(attr_value)

                    node_attrs[node].append(node_attr_keys[value].index(attr_value))
                
            if pdb_code not in node_to_pdbs[node]:
                node_to_pdbs[node].append(pdb_code)
        
        for u, v, data in g.edges(data=True):
            edge = (u, v)
            
            if edge not in edge_to_pdbs:
                edge_to_pdbs[edge] = []
                edge_attrs[edge] = []
                
                for value in edge_attr_keys:
                    attr_value = data[value]
                    if isinstance(attr_value, np.ndarray):
                        attr_value = tuple(attr_value)
                        
                    if attr_value not in edge_attr_keys[value]:
                        edge_attr_keys[value].append(attr_value)
                    
                    edge_attrs[edge].append(edge_attr_keys[value].index(attr_value))
            
            if pdb_code not in edge_to_pdbs[edge]:
                edge_to_pdbs[edge].append(pdb_code)

    edge_id = 0
    edge_to_id = {}
    for e in edge_to_pdbs:
        edge_to_id[e] = edge_id
        pdbs = edge_to_pdbs[e]
        for pdb_code in pdbs:
            pdb_to_view[pdb_code][1][0].add(edge_id)
            pdb_to_view[pdb_code][3].append(edge_attrs[e])
            
        edge_id += 1

    del edge_to_pdbs

    node_id = 0
    node_to_id = {}  
    for u in node_to_pdbs:
        node_to_id[u] = node_id
        pdbs = node_to_pdbs[u]
        for pdb_code in pdbs:
            pdb_to_view[pdb_code][0][0].add(node_id)
            pdb_to_view[pdb_code][2].append(node_attrs[u])
            
        node_id += 1

    del node_to_pdbs

    node_to_id = bidict(node_to_id)  
    edge_to_id = bidict(edge_to_id)

    for pdb_code, view in pdb_to_view.items():
        original_graph = protein_graphs[pdb_code]

        nodes = [node_to_id.inverse[node_id] for node_id in view[0][0]]  
        edges = [edge_to_id.inverse[edge_id] for edge_id in view[1][0]]
        extracted_graph = nx.Graph()
        extracted_graph.update(edges=edges, nodes=nodes)

        for node in nodes:
            if node in node_attrs:
                for i, key in enumerate(node_attr_keys):
                    index = node_attrs[node][i]
                    value = node_attr_keys[key][index]

                    if isinstance(value, tuple) and len(value) == 1:
                        value = np.array(value[0])
                    elif isinstance(value, tuple) and len(value) == 3:
                        value = pd.Series(value[0], name=value[1], index=value[2])
                        
                    extracted_graph.nodes[node][key] = value

        for u, v in edges:
            if (u, v) in edge_attrs:
                for i, key in enumerate(edge_attr_keys):
                    index = edge_attrs[(u, v)][i]
                    value = edge_attr_keys[key][index]
                    
                    if isinstance(value, tuple):
                        value = np.array(value)
                        
                    extracted_graph.edges[u, v][key] = value

        assert nx.utils.nodes_equal(extracted_graph.nodes, original_graph.nodes)
        assert nx.utils.edges_equal(extracted_graph.edges, original_graph.edges)

    print("nodeSizeCompressed", (asizeof.asizeof(node_to_id) + asizeof.asizeof(node_attrs) + asizeof.asizeof(node_attr_keys)) / 1024 / 1024, "MB")
    print("edgeSizeCompressed", (asizeof.asizeof(edge_to_id) + asizeof.asizeof(edge_attrs) + asizeof.asizeof(edge_attr_keys)) / 1024 / 1024, "MB")
    print("pdbToViewSizeCompressed", asizeof.asizeof(pdb_to_view) / 1024 / 1024, "MB")

    return PDBGraphStoreBitmap(node_to_id, edge_to_id, pdb_to_view, edge_attr_keys, node_attr_keys, edge_attrs, node_attrs)



############################################################################################

class PDBGraphStoreBitmap:

    def __init__(self, isolated_node_to_id, edge_to_id, pdb_to_view, edge_attr_keys, isolated_node_attr_keys, edge_attrs, isolated_node_attrs):
        self.isolated_node_to_id = isolated_node_to_id
        self.edge_to_id = edge_to_id
        self.pdb_to_view = pdb_to_view
        self.isolated_node_attr_keys = isolated_node_attr_keys
        self.edge_attr_keys = edge_attr_keys
        self.edge_attrs = edge_attrs
        self.isolated_node_attrs = isolated_node_attrs


    def extract_pdb_graph(self, pdb_code):
        view = self.pdb_to_view.get(pdb_code, None)
        
        def union_bitmaps(bms):
            union_bm = BitMap64()
            if bms is None: 
                return union_bm
            for bm in bms: 
                union_bm = union_bm | bm
            return union_bm
        
        isolated_nodes = [self.isolated_node_to_id.inverse[node_id] for node_id in union_bitmaps(view[0])]
        edges = [self.edge_to_id.inverse[edge_id] for edge_id in union_bitmaps(view[1])]
        
        extracted_graph = nx.Graph()
        extracted_graph.update(edges=edges, nodes=isolated_nodes)
        
        for node in isolated_nodes:
            if node in self.isolated_node_attrs:
                for i, key in enumerate(self.isolated_node_attr_keys):
                    index = self.isolated_node_attrs[node][i]
                    value = self.isolated_node_attr_keys[key][index]
                    
                    if isinstance(value, tuple) and len(value) == 1:
                        value = np.array(value[0])
                    elif isinstance(value, tuple) and len(value) == 3:
                        value = pd.Series(value[0], name=value[1], index=value[2])
                        
                    extracted_graph.nodes[node][key] = value
        
        for u, v in edges:
            if (u, v) in self.edge_attrs:
                for i, key in enumerate(self.edge_attr_keys):
                    index = self.edge_attrs[(u, v)][i]
                    value = self.edge_attr_keys[key][index]
                    
                    if isinstance(value, tuple):
                        value = np.array(value)  
                        
                    extracted_graph.edges[u, v][key] = value
        
        return extracted_graph

    def run_tree_compression(self):
        # build pairs jaccard
        def get_max_pair(pdb_to_bitmaps):
            best_pair = None
            best_jacc = None
            for pdb_code1, bms1 in pdb_to_bitmaps.items():
                bm_1 = bms1[0]
                for pdb_code2, bms2 in pdb_to_bitmaps.items():
                    if pdb_code1 >= pdb_code2: continue
                    bm_2 = bms2[0]
                    jacc = bm_1.intersection_cardinality(bm_2)
                    if best_jacc is None or jacc > best_jacc:
                        best_jacc = jacc
                        best_pair = (pdb_code1, pdb_code2)
            return best_pair

        def get_optimized_bitmaps(reconstructions):
            while len(reconstructions) >= 2:
                pdb_code1, pdb_code2 = get_max_pair(reconstructions)
                # print(pdb_code1, pdb_code2)
                left = reconstructions[pdb_code1][0]
                right = reconstructions[pdb_code2][0]
                intersection = left & right
                delta_1 = left - intersection
                delta_2 = right - intersection

                left_child = (delta_1, reconstructions[pdb_code1][1])
                right_child = (delta_2, reconstructions[pdb_code2][1])
                reconstructions[pdb_code1 + '_' + pdb_code2] = (
                    intersection, {pdb_code1: left_child, pdb_code2: right_child})

                del reconstructions[pdb_code1]
                del reconstructions[pdb_code2]

            def key_in_keys(key, keys):
                for k in keys:
                    if key in k: return k
                return None

            root_tree_node = list(reconstructions.items())[0]

            pdb_to_view_opt = {}

            for pdb_code in self.pdb_to_view:
                tree_node = root_tree_node
                bms = [tree_node[1][0]] if len(tree_node[1][0]) > 0 else []
                bm_height = 1
                key = key_in_keys(pdb_code, tree_node[1][1])
                while key:
                    bm_height += 1
                    tree_node = (key, tree_node[1][1][key])
                    if len(tree_node[1][0]) > 0:
                        tree_node[1][0].run_optimize()
                        bms.append(tree_node[1][0])
                    key = key_in_keys(pdb_code, tree_node[1][1])
                if len(bms) > 0:
                    pdb_to_view_opt[pdb_code] = bms

            return pdb_to_view_opt

        edge_bitmaps = get_optimized_bitmaps(
            {pdb_code: (view[1][0], {}) for pdb_code, view in self.pdb_to_view.items()})
        isolated_node_bitmaps = get_optimized_bitmaps(
            {pdb_code: (view[0][0], {}) for pdb_code, view in self.pdb_to_view.items()})

        # TODO: assertion code (remove) {
        for pdb_code, view in self.pdb_to_view.items():
            union_bm1 = BitMap64()
            for bm in isolated_node_bitmaps.get(pdb_code, BitMap64()):
                union_bm1 = union_bm1 | bm

            assert union_bm1.symmetric_difference_cardinality(view[0][0]) == 0

            union_bm1 = BitMap64()
            for bm in edge_bitmaps[pdb_code]:
                union_bm1 = union_bm1 | bm

            assert union_bm1.symmetric_difference_cardinality(view[1][0]) == 0
        # }

        # replace existing views with optimized ones
        self.pdb_to_view = {pdb_code: (isolated_node_bitmaps.get(pdb_code, None), edge_bitmaps.get(pdb_code, None)) for
                            pdb_code in self.pdb_to_view}
        self.pdb_to_view = {k: v for k, v in self.pdb_to_view.items() if (k, v) != (None, None)}


def print_graph(protein_graphs):
    return (graph.nodes(data=True))
    

if __name__ == "__main__":
    import os
    import metadata
    from graphein.protein.config import ProteinGraphConfig
    from graphein.protein.edges.distance import add_hydrogen_bond_interactions, add_peptide_bonds
    from graphein.protein.graphs import construct_graph
    from graphein.protein.utils import download_pdb
    import networkx as nx


    file_path = os.path.dirname(os.path.realpath(metadata.__file__))
    print(file_path)

    params_to_change = {"granularity": "atom", "edge_construction_functions": [add_atomic_edges]}

    config = ProteinGraphConfig(**params_to_change)
    # print(config.model_dump())

    pdb_codes = []
    with open(f'{file_path}/../../data/soybean_ppigremlin.txt', 'r') as f:
        for line in f:
            pdb_codes.append(line.strip())

    protein_graphs_with_data = {}

    i = 0
    for pdb_code in pdb_codes:
        print(i, pdb_code)
        i += 1

        if os.path.exists(f"{file_path}/../../data/pdb_files/{pdb_code}.pdb"):
            pdb_file = os.path.abspath(f"{file_path}/../../data/pdb_files/{pdb_code}.pdb")
        else:
            pdb_file = download_pdb(pdb_code, f"{file_path}/../../data/pdb_files/")
            if pdb_file is None:
                print(f"Failed to download {pdb_code}")
                continue
        
        graph = construct_graph(config=config, path=pdb_file)
        graph.graph.clear()
        print(graph)
        protein_graphs_with_data[pdb_code] = graph  # Store graph

    v_size = 0
    e_size = 0
    g_size = 0


    for k, v in protein_graphs_with_data.items():
        v.graph.clear()
        v_size += asizeof.asizeof(v._node) / 1024 / 1024
        e_size += asizeof.asizeof(v._adj) / 1024 / 1024
        g_size += asizeof.asizeof(v) / 1024 / 1024

    print("VertexSize", v_size, "MB")
    print("EdgeSize", e_size, "MB")
    print("UncompressedBytes", g_size, "MB")

    ob = compress_with_composition(protein_graphs_with_data)
    print("CompressedBytes", asizeof.asizeof(ob) / 1024 / 1024)

    an1 = ob.extract_pdb_graph("1AN1")

    assert nx.utils.nodes_equal(an1._node, protein_graphs_with_data["1AN1"]._node)
    assert nx.utils.edges_equal(an1._adj, protein_graphs_with_data["1AN1"]._adj)


    # print("CompressedBytesSer", len(pickle.dumps(pdb_graph_store_bitmap)) / 1024 / 1024)
    # pdb_graph_store_bitmap.run_tree_compression()
    # print("CompressedAndOptimizedBytes", asizeof.asizeof(pdb_graph_store_bitmap) / 1024 / 1024)
    # print("CompressedAndOptimizedBytesSer", len(pickle.dumps(pdb_graph_store_bitmap)) / 1024 / 1024)

    # # TODO: assertion code (remove) {

    # start_time = time.time()
    # for pdb_code, protein_graph in protein_graphs.items():
    #     continue
    # elapsed_time = time.time() - start_time
    # print("RuntimeIterationUncompressed", elapsed_time)

    # start_time = time.time()
    # for pdb_code, protein_graph in protein_graphs.items():
    #     computed_graph = pdb_graph_store_bitmap.extract_pdb_graph(pdb_code)
    #     assert nx.utils.graphs_equal(computed_graph, protein_graph)
    # elapsed_time = time.time() - start_time
    # print("RuntimeIterationCompressed", elapsed_time)

    # # }
