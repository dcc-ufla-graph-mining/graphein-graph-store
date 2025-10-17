import networkx as nx
import numpy as np
from bidict import bidict
from pympler import asizeof
from pyroaring import BitMap, BitMap64
from sortedcontainers import SortedSet
import pandas as pd
import edge_functions_Model as edgeModel
import struct

edge_functions_func_list = []
edge_functions_str_list = []

def _initialize_structs():
    edge_to_pdbs = {}
    node_to_pdbs = {}
    edge_kinds = {}
    node_attrs = {}
    pdb_to_nodes = {}
    pdb_to_edges = {}
    pdb_codes_config = {}
    all_pdb_codes = SortedSet()

    edge_distances_temp = {}

    return edge_to_pdbs,node_to_pdbs,edge_kinds,node_attrs,pdb_to_nodes,pdb_to_edges,pdb_codes_config, all_pdb_codes,edge_distances_temp

def _extract_attribute_keys():
    edge_attr_keys_list = ["kind", "distance"]
    node_attr_keys_list = ["chain_id", "residue_name", "residue_number", "atom_type", "element_symbol", "coords", "b_factor", "meiler"]

    node_attr_keys = {}

    for key in node_attr_keys_list:
        node_attr_keys[key] = SortedSet()

    edge_kind_keys = SortedSet(set(edge_functions_str_list))

    return edge_kind_keys, node_attr_keys

def _process_node_attributes(node, graph, node_attr_keys):
    attr_indexes = []

    for value in node_attr_keys:
        attr_value = graph.nodes[node][value]

        if isinstance(attr_value, pd.Series):
            attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
        
        if isinstance(attr_value, np.ndarray):
            attr_value = tuple([tuple(attr_value)])

        node_attr_keys[value].add(attr_value)
        idx = node_attr_keys[value].index(attr_value)

        attr_indexes.append(idx)

    return attr_indexes

#Todo: implement all_pdb_code logic for edge distances
def _process_edge_kinds(edge_data, edge_kind_keys):
    attr_indexes = []

    attr_kind_value = list(edge_data["kind"]) 

    if attr_kind_value is None:
        raise ValueError("Edge kind attribute should not be None")
    
    kind_indexes = set()

    for kind in attr_kind_value:
        if kind not in edge_functions_func_list:
            kind_indexes.add(edge_kind_keys.index(kind))

    attr_indexes.append(kind_indexes)

    # distance_indexes = {}
    # distance_indexes[pdb_code] = process_edge_distances(edge_data, edge_kind_keys, pdb_code)

    # attr_indexes.append(distance_indexes)

    return attr_indexes

#todo: implement all_pdb_codes index logic for distance
def _check_if_edge_attr_matches(edge, data, edge_kinds, edge_kind_keys):
    kinds = set(filter(lambda x: x , data["kind"]))
    for k in kinds:
        edge_kind_keys.add(k)
    idx = set([edge_kind_keys.index(k) for k in kinds])
    if len(idx) > 0:
        for i in idx:
            edge_kinds[edge].append(i)

def _process_nodes(protein_graphs, node_to_pdbs, node_attrs, node_attr_keys, pdb_to_nodes, all_pdb_codes):
    for pdb_code, graphs in protein_graphs.items():
        all_pdb_codes.add(pdb_code)
        pdb_idx = all_pdb_codes.index(pdb_code)

        for g in graphs:
            pdb_to_nodes[pdb_idx] = BitMap64()
            for node in g.nodes():
                if node not in node_to_pdbs:
                    node_to_pdbs[node] = []
                    attr_indexes = _process_node_attributes(node, g, node_attr_keys)
                    node_attrs[node] = attr_indexes

                if pdb_code not in node_to_pdbs[node]:
                    node_to_pdbs[node].append(pdb_code)    

def _process_edges(protein_graphs, edge_to_pdbs, edge_kinds, edge_kind_keys, pdb_to_edges, all_pdb_codes, edge_distances_temp):
    for pdb_code, graphs in protein_graphs.items():

        pdb_idx = all_pdb_codes.index(pdb_code)

        for g in graphs:
            if pdb_idx not in pdb_to_edges:
                pdb_to_edges[pdb_idx] = BitMap64()
            
            for u, v, data in g.edges.data():
                edge = (u, v)

                if edge not in edge_to_pdbs:
                    edge_to_pdbs[edge] = []
                    try:
                        kind_indexes = _process_edge_kinds(data, edge_kind_keys)
                        # print(attr_indexes)
                    except ValueError as e:
                        print("error at process edges: ", e)

                    edge_kinds[edge] = kind_indexes
                else:
                    _check_if_edge_attr_matches(edge, data, edge_kinds, edge_kind_keys)

                if pdb_code not in edge_to_pdbs[edge]:
                    edge_to_pdbs[edge].append(pdb_code)
                    edge_distances_temp.setdefault(edge, {})[pdb_code] = data["distance"]
                
                #ao inves de append somente o pdb_code, apeend uma tupla (pdb_code, distance)

def _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes, all_pdb_codes, edge_distances_temp):
    edge_to_id = {}
    edge_distances = []

    for edge_id, e in enumerate(edge_to_pdbs):
        edge_to_id[e] = edge_id
        pdbs = edge_to_pdbs[e]
        for pdb_code in pdbs:
            pdb_idx = all_pdb_codes.index(pdb_code)

            edge_pair = tuple([edge_id, edge_distances_temp[e][pdb_code]])

            edge_pair = struct.pack("ld", *edge_pair)

            edge_distances.append(edge_pair)

            pdb_to_edges[pdb_idx].add(len(edge_distances)-1)

    node_to_id = {}

    for node_id, u in enumerate(node_to_pdbs):
        node_to_id[u] = node_id
        pdbs = node_to_pdbs[u]
        for pdb_code in pdbs:
            pdb_idx = all_pdb_codes.index(pdb_code)
            pdb_to_nodes[pdb_idx].add(node_id)

    node_to_id = bidict(node_to_id)
    edge_to_id = bidict(edge_to_id)

    return node_to_id, edge_to_id, edge_distances

def _process_pdb_codes_config(protein_graphs, pdb_codes_config):
    def add_edge_construction_function(fun):
        if fun not in pdb_codes_config[pdb_code]["edge_construction_functions"]:
            pdb_codes_config[pdb_code]["edge_construction_functions"].append(fun)

    for pdb_code, g in protein_graphs.items():
        if pdb_code not in pdb_codes_config:
            pdb_codes_config[pdb_code] = {}
        if "edge_construction_functions" not in pdb_codes_config[pdb_code]:
            pdb_codes_config[pdb_code]["edge_construction_functions"] = []

        for graph in g:
            for fun in graph.graph["config"].edge_construction_functions:
                if fun in edge_functions_func_list:
                    add_edge_construction_function(fun.__name__)    

def _reconstruct_node_attributes(extracted_graph, nodes, node_attrs, node_attr_keys):
    for node in nodes:
        if node in node_attrs:
            for i, key in enumerate(node_attr_keys):
                index = node_attrs[node][i]
                value = node_attr_keys[key][index]

                if isinstance(value, tuple) and len(value) == 1:
                    value = np.array(value[0])
                if isinstance(value, tuple) and len(value) == 3:
                    value = pd.Series(value[0], name=value[1], index=value[2])

                extracted_graph.nodes[node][key] = value

def _reconstruct_edge_attributes(extracted_graph, edges, edge_kinds, edge_kind_keys, edge_funcs, pdb_to_edges, pdb_code, edge_distances, edge_to_id, all_pdb_codes):
    pdb_idx = all_pdb_codes.index(pdb_code)
    for u, v in edges:
        edge = (u, v)
        
        if edge not in edge_kinds:
            continue
            
        kind_indexes_list = edge_kinds[edge]
        
        if len(kind_indexes_list) > 0:
            kind_indexes = kind_indexes_list[0]  
            kind_names = [edge_kind_keys[k] for k in kind_indexes]
            kind_names = list(filter(lambda x: edgeModel.edge_functions_dict[x] in edge_funcs, kind_names))
            
            if len(kind_names) > 0:
                if not extracted_graph.has_edge(*edge):
                    extracted_graph.add_edge(*edge)
                
                extracted_graph.edges[edge].setdefault("kind", set())
                
                for kind in kind_names:
                    extracted_graph.edges[edge]["kind"].add(kind)

    for edge_dist_idx in pdb_to_edges[pdb_idx]:
        edge_pair = edge_distances[edge_dist_idx]

        edge_pair = struct.unpack("ld", edge_pair)

        edge = edge_to_id.inverse[edge_pair[0]]
        extracted_graph.edges[edge]["distace"] = edge_pair[1]

def _reconstruct_and_validate_graphs(protein_graphs,
                                    node_to_id,
                                    edge_to_id,
                                    pdb_to_nodes,
                                    pdb_to_edges,
                                    node_attrs,
                                    edge_kinds,
                                    node_attr_keys,
                                    edge_kind_keys,
                                    edge_distances,
                                    all_pdb_codes):

    for pdb_code in protein_graphs:
        pdb_graphs = protein_graphs[pdb_code]

        pdb_idx = all_pdb_codes.index(pdb_code)

        for graph in pdb_graphs:
            original_graph = graph.copy()

            nodes = [node_to_id.inverse[node_id] for node_id in pdb_to_nodes[pdb_idx]]

            def get_edge_id(edge_distances_idx):
                edge_pair = edge_distances[edge_distances_idx]
                edge_pair = struct.unpack("ld", edge_pair)
                edge_id = edge_pair[0]
                return edge_id
            
            edges = [edge_to_id.inverse[get_edge_id(edge_distances_idx)] for edge_distances_idx in pdb_to_edges[pdb_idx]]

            extracted_graph = nx.Graph()
            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=nodes)

            _reconstruct_node_attributes(extracted_graph, nodes, node_attrs, node_attr_keys)
            _reconstruct_edge_attributes(extracted_graph, edges, edge_kinds, edge_kind_keys, original_graph.graph["config"].edge_construction_functions, pdb_to_edges, pdb_code, edge_distances, edge_to_id, all_pdb_codes)

            print(f"original {pdb_code} graph edge funcs: {original_graph.graph['config'].edge_construction_functions}")
            print(f"number of nodes in original graph: {len(original_graph.nodes())}")
            print(f"number of nodes in extracted graph: {len(extracted_graph.nodes())}")

            print(f"number of edges in original graph: {len(original_graph.edges())}")
            print(f"number of edges in extracted graph: {len(extracted_graph.edges())}")

            for e in original_graph.edges:
                print(f"original: {original_graph.edges[e]}")
                print(f"extracted: {extracted_graph.edges[e]}")

def compress_pdb_graphs(protein_graphs):
    print("reached compress")
    edge_to_pdbs, \
    node_to_pdbs, \
    edge_kinds,   \
    node_attrs,   \
    pdb_to_nodes, \
    pdb_to_edges, \
    pdb_codes_config, \
    all_pdb_codes, \
    edge_distances_temp = _initialize_structs()

    global edge_functions_func_list
    global edge_functions_str_list

    edge_functions_func_list = [v for _, v in edgeModel.edge_functions_dict.items()]
    edge_functions_str_list = [k for k, _ in edgeModel.edge_functions_dict.items()]

    edge_kind_keys, node_attr_keys = _extract_attribute_keys()

    _process_nodes(protein_graphs, node_to_pdbs, node_attrs, node_attr_keys, pdb_to_nodes, all_pdb_codes)
    _process_edges(protein_graphs, edge_to_pdbs, edge_kinds, edge_kind_keys, pdb_to_edges, all_pdb_codes, edge_distances_temp)

    print(all_pdb_codes)

    node_to_id, edge_to_id, edge_distances = _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes, all_pdb_codes, edge_distances_temp)

    _process_pdb_codes_config(protein_graphs, pdb_codes_config)

    del edge_to_pdbs
    del node_to_pdbs
    del edge_distances_temp

    _reconstruct_and_validate_graphs(protein_graphs,
                                    node_to_id,
                                    edge_to_id,
                                    pdb_to_nodes,
                                    pdb_to_edges,
                                    node_attrs,
                                    edge_kinds,
                                    node_attr_keys,
                                    edge_kind_keys,
                                    edge_distances,
                                    all_pdb_codes)

    return node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_kinds, node_attr_keys, edge_kind_keys, edge_distances,all_pdb_codes