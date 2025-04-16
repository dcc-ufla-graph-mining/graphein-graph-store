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
#the ordered_set lib has an data type named OrderedSet
#which is a set that keeps the order of the elements
#because of this, we can search for the index or for the key, both with o(1)
#https://pypi.org/project/ordered-set/
from ordered_set import OrderedSet


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
        edge_attr_keys[key] = OrderedSet()
    for key in node_attr_keys_list: 
        node_attr_keys[key] = OrderedSet()

    del edge_attr_keys_list
    del node_attr_keys_list
    
    pdb_to_view = {}

    for pdb_code, g in protein_graphs.items():
        pdb_to_view[pdb_code] = ([BitMap64()], [BitMap64()], [], [])  # nodes, edges, node attributes, edge attributes
        
        for node in g.nodes():
            if node not in node_to_pdbs:
                node_to_pdbs[node] = []

                attr_indexes = []
                
                for value in node_attr_keys:
                    attr_value = g.nodes[node][value]

                    if isinstance(attr_value, pd.Series):
                        attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
                    elif isinstance(attr_value, np.ndarray):
                        attr_value = tuple([tuple(attr_value)])

                    attr_indexes.append(node_attr_keys[value].add(attr_value))

                    node_attrs[node] = attr_indexes
                
            if pdb_code not in node_to_pdbs[node]:
                node_to_pdbs[node].append(pdb_code)
        
        for u, v, data in g.edges(data=True):
            edge = (u, v)
            
            if edge not in edge_to_pdbs:
                edge_to_pdbs[edge] = []
                attr_indexes = []
                
                for value in edge_attr_keys:
                    attr_value = data[value]
                    if isinstance(attr_value, np.ndarray):
                        attr_value = tuple(attr_value)
                    elif isinstance(attr_value, (list, set)):
                        attr_value = tuple(attr_value)
                        
                    attr_indexes.append(edge_attr_keys[value].add(attr_value))
                    edge_attrs[edge] = attr_indexes
            
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
                        value = set(value)
                        
                    extracted_graph.edges[u, v][key] = value

        # assert nx.utils.nodes_equal(extracted_graph.nodes(data=True), original_graph.nodes(data=True))
        try:
            assert nx.utils.edges_equal(extracted_graph._adj, original_graph._adj)
        except AssertionError as e:
            print("error in edges")
            print(nx.difference(extracted_graph, original_graph).nodes(data=True))
            for (u, v, attr) in extracted_graph.edges.data():
                print(u, v, attr)
                for k, t in attr.items():
                    print(type(t))
                
                break

            for (u, v, attr) in original_graph.edges.data():  
                print(u, v, attr) 
                for k, t in attr.items():
                    print(type(t))
                break
            break

        try:    
            assert nx.utils.nodes_equal(extracted_graph._node, original_graph._node)
        except AssertionError as e:
            print("error in nodes")
            for node, attr in extracted_graph.nodes(data=True):
                print(node, attr)
                for k, t in attr.items():
                    print(type(t))
                break

            for node, attr in original_graph.nodes(data=True):
                print(node, attr)
                for k, t in attr.items():
                    print(type(t))
                break
            break
    
    return PDBGraphStoreBitmap(node_to_id, edge_to_id, pdb_to_view, edge_attr_keys, node_attr_keys)


class PDBGraphStoreBitmap:

    def __init__(self, node_to_id, edge_to_id, pdb_to_view, edge_attr_keys, node_attr_keys):
        self.node_to_id = node_to_id
        self.edge_to_id = edge_to_id
        self.pdb_to_view = pdb_to_view
        self.node_attr_keys = node_attr_keys
        self.edge_attr_keys = edge_attr_keys

    # return the value in MB
    def calculate_graph_complete_space_size(self):
        return (asizeof.asizeof(self.node_to_id) + asizeof.asizeof(self.edge_to_id) + \
               asizeof.asizeof(self.pdb_to_view) + asizeof.asizeof(self.edge_attr_keys) + \
               asizeof.asizeof(self.node_attr_keys)) / 1024 / 1024
    
    # return the value in MB
    def calculate_graph_structure_space_size(self):
        return (asizeof.asizeof(self.node_to_id) + asizeof.asizeof(self.edge_to_id) + \
               asizeof.asizeof(self.pdb_to_view)) / 1024 / 1024

    # return the value in MB
    def calculate_edges_complete_space_size(self):
        total_size = 0
        for pdb_code, view in self.pdb_to_view.items():
            total_size += asizeof.asizeof(view[1]) + asizeof.asizeof(view[3])
        
        return (total_size + asizeof.asizeof(self.edge_attr_keys) + asizeof.asizeof(self.edge_to_id)) / 1024 / 1024
    
    # return the value in MB
    def calculate_nodes_complete_space_size(self):
        total_size = 0
        for pdb_code, view in self.pdb_to_view.items():
            total_size += asizeof.asizeof(view[0]) + asizeof.asizeof(view[2])
        
        return (total_size + asizeof.asizeof(self.node_attr_keys) + asizeof.asizeof(self.node_to_id)) / 1024 / 1024

    def extract_pdb_graph(self, pdb_code):
        view = self.pdb_to_view.get(pdb_code, None)
        
        def union_bitmaps(bms):
            union_bm = BitMap64()
            if bms is None: 
                return union_bm
            for bm in bms: 
                union_bm = union_bm | bm
            return union_bm
        
        try:
            isolated_nodes = [self.node_to_id.inverse[node_id] for node_id in union_bitmaps(view[0])]
            edges = [self.edge_to_id.inverse[edge_id] for edge_id in union_bitmaps(view[1])]
            
            extracted_graph = nx.Graph()
            extracted_graph.update(edges=edges, nodes=isolated_nodes)

            node_attrs = self.pdb_to_view[pdb_code][2]
            edge_attrs = self.pdb_to_view[pdb_code][3]
            
            for node in isolated_nodes:
                if node in node_attrs:
                    for i, key in enumerate(self.node_attr_keys):
                        index = node_attrs[node][i]
                        value = self.node_attr_keys[key][index]
                        
                        if isinstance(value, tuple) and len(value) == 1:
                            value = np.array(value[0])
                        elif isinstance(value, tuple) and len(value) == 3:
                            value = pd.Series(value[0], name=value[1], index=value[2])
                            
                        extracted_graph.nodes[node][key] = value
            
            for u, v in edges:
                if (u, v) in edge_attrs:
                    for i, key in enumerate(self.edge_attr_keys):
                        index = edge_attrs[(u, v)][i]
                        value = self.edge_attr_keys[key][index]
                        
                        if isinstance(value, tuple):
                            value = np.array(value)  
                            
                        extracted_graph.edges[u, v][key] = value
        except Exception as e:
            print(f"Error extracting graph for {pdb_code}: {e}")
            extracted_graph = nx.Graph()
            extracted_graph.update(edges=[], nodes=[])
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
    import time


    file_path = os.path.dirname(os.path.realpath(metadata.__file__))
    print(file_path)

    params_to_change = {"granularity": "atom", "edge_construction_functions": [add_atomic_edges]}

    config = ProteinGraphConfig(**params_to_change)
    # print(config.model_dump())

    pdb_codes = []
    with open(f'{file_path}/../../data/bcl_ppigremlin.txt', 'r') as f:
        for line in f:
            pdb_codes.append(line.strip())

    protein_graphs_with_data = {}

    time_begin = time.time()

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

    time_end = time.time()

    print("Time to construct graphs:", time_end - time_begin)
    print("Number of graphs:", len(protein_graphs_with_data))

    v_size = 0
    e_size = 0

    for k, v in protein_graphs_with_data.items():
        v.graph.clear()
        v_size += asizeof.asizeof(v._node) / 1024 / 1024
        e_size += asizeof.asizeof(v._adj) / 1024 / 1024

    time_begin = time.time()

    global_graph_obj = compress_with_composition(protein_graphs_with_data)

    time_end = time.time()
    print("Time to compress:", time_end - time_begin)
    
    time_begin = time.time()
    an1 = global_graph_obj.extract_pdb_graph("1BXL")
    time_end = time.time()

    print("Time to extract:", time_end - time_begin)

    assert nx.utils.nodes_equal(an1._node, protein_graphs_with_data["1BXL"]._node)
    assert nx.utils.edges_equal(an1._adj, protein_graphs_with_data["1BXL"]._adj)

    print("uncompressed complete graph size", asizeof.asizeof(protein_graphs_with_data) / 1024 / 1024)
    print("uncompressed edge size", e_size)
    print("uncompressed node size", v_size)
    print("compressed graph complete space size", global_graph_obj.calculate_graph_complete_space_size())
    print("compressed graph structure space size", global_graph_obj.calculate_graph_structure_space_size())
    print("compressed edge size", global_graph_obj.calculate_edges_complete_space_size())
    print("compressed node size", global_graph_obj.calculate_nodes_complete_space_size())
    print("global graph object size", asizeof.asizeof(global_graph_obj) / 1024 / 1024)