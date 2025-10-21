import networkx as nx
import numpy as np
from bidict import bidict
from pyroaring import BitMap64
from sortedcontainers import SortedSet
import pandas as pd
import edge_functions_Model as edgeModel
import struct

from Indexed_set import IndexedSet

edge_functions_func_list = []
edge_functions_str_list = []

def _initialize_structs():
    edge_to_pdbs = {}
    node_to_pdbs = {}
    edge_kinds = {}
    node_attrs_global = {}
    pdb_to_nodes = {}
    pdb_to_edges = {}
    pdb_codes_config = {}
    node_attrs_unique_temp = {}
    all_pdb_codes = SortedSet()

    edge_distances_temp = {}

    return edge_to_pdbs,node_to_pdbs,edge_kinds,node_attrs_global,pdb_to_nodes,pdb_to_edges,pdb_codes_config, all_pdb_codes, node_attrs_unique_temp, edge_distances_temp

def _extract_attribute_keys():
    node_attr_keys_list = ["chain_id", "residue_name", "residue_number", "atom_type", "element_symbol", "coords", "b_factor", "meiler"]

    node_attr_keys = {}

    for key in node_attr_keys_list:
        node_attr_keys[key] = IndexedSet()

    edge_kind_keys = SortedSet(set(edge_functions_str_list))

    return edge_kind_keys, node_attr_keys

def _process_node_attributes(node, graph, node_attr_keys):
    attr_indexes_global = []
    attr_indexes_unique = []

    for i, value in enumerate(node_attr_keys):
        attr_value = graph.nodes[node][value]

        if i == 7: #attr: meiler; type: pd.Series
            attr_value = tuple([tuple(attr_value.tolist()), tuple(attr_value.name)])
        
        if i == 5: ##attr: coords; type: np.array
            attr_value = tuple([tuple(attr_value)])

        idx = node_attr_keys[value].add(attr_value)

        if i in [0, 1, 3, 4, 7]:
            attr_indexes_global.append(idx)
        else:
            attr_indexes_unique.append(idx)

    return attr_indexes_global, attr_indexes_unique


def _process_edge_kinds(edge_data, edge_kind_keys):
    attr_kind_value = list(edge_data["kind"]) 

    if attr_kind_value is None:
        raise ValueError("Edge kind attribute should not be None")
    
    kind_indexes = set()

    for kind in attr_kind_value:
        kind_indexes.add(edge_kind_keys.index(kind))

    return kind_indexes

def _check_if_edge_attr_matches(edge, data, edge_kinds, edge_kind_keys):
    kinds = set(filter(lambda x: x , data["kind"]))
    for k in kinds:
        edge_kind_keys.add(k)
    idx = [edge_kind_keys.index(k) for k in kinds]
    if idx:
        for i in idx:
            edge_kinds[edge].add(i)

def _process_nodes(protein_graphs, node_to_pdbs, node_attrs_global, node_attrs_unique_temp, node_attr_keys, pdb_to_nodes, all_pdb_codes):
    for pdb_code, graphs in protein_graphs.items():
        all_pdb_codes.add(pdb_code)
        pdb_idx = all_pdb_codes.index(pdb_code)
        
        if pdb_idx not in pdb_to_nodes:
            pdb_to_nodes[pdb_idx] = BitMap64()

        for g in graphs:
            for node in g.nodes():
                attr_indexes_global, attr_indexes_unique = _process_node_attributes(node, g, node_attr_keys)
                if node not in node_to_pdbs:
                    node_to_pdbs[node] = []
                    node_attrs_global[node] = attr_indexes_global
                    node_attrs_unique_temp.setdefault(node, {})

                node_attrs_unique_temp[node][pdb_code] = attr_indexes_unique

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


def _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes, all_pdb_codes, edge_distances_temp, node_attrs_unique_temp):
    edge_to_id = {}
    edge_distances = []
    node_attrs_unique = []

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
            # node_id, residue_number_idx, coords_idx, b_factor_idx 
            unique_attrs = node_attrs_unique_temp[u][pdb_code]
            node_pair = tuple([node_id, unique_attrs[0], unique_attrs[1], unique_attrs[2]])
            node_pair = struct.pack("llll", *node_pair)

            node_attrs_unique.append(node_pair)
            pdb_to_nodes[pdb_idx].add(len(node_attrs_unique)-1)

    node_to_id = bidict(node_to_id)
    edge_to_id = bidict(edge_to_id)

    return node_to_id, edge_to_id, node_attrs_unique, edge_distances

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

def _reconstruct_node_attributes(extracted_graph, node_attr_keys, node_attrs_global, node_attrs_unique, node_to_id, pdb_to_nodes, all_pdb_codes, pdb_code):
    pdb_idx = all_pdb_codes.index(pdb_code)
    
    for node_attr_idx in pdb_to_nodes[pdb_idx]:
        node_pair = struct.unpack("llll", node_attrs_unique[node_attr_idx])
        node_id = node_pair[0]
        node = node_to_id.inverse[node_id]
        
        chain_id_idx = node_attrs_global[node][0]
        residue_name_idx = node_attrs_global[node][1]
        atom_type_idx = node_attrs_global[node][2]
        element_symbol_idx = node_attrs_global[node][3]
        meiler_idx = node_attrs_global[node][4]

        residue_number_idx = node_pair[1]  # unique_attrs[0]
        coords_idx = node_pair[2]          # unique_attrs[1]
        b_factor_idx = node_pair[3]        # unique_attrs[2]
        
        extracted_graph.nodes[node]["chain_id"] = node_attr_keys["chain_id"].get(chain_id_idx)
        extracted_graph.nodes[node]["residue_name"] = node_attr_keys["residue_name"].get(residue_name_idx)
        extracted_graph.nodes[node]["residue_number"] = node_attr_keys["residue_number"].get(residue_number_idx)
        extracted_graph.nodes[node]["atom_type"] = node_attr_keys["atom_type"].get(atom_type_idx)
        extracted_graph.nodes[node]["element_symbol"] = node_attr_keys["element_symbol"].get(element_symbol_idx)
        extracted_graph.nodes[node]["coords"] = np.array(node_attr_keys["coords"].get(coords_idx))
        extracted_graph.nodes[node]["b_factor"] = node_attr_keys["b_factor"].get(b_factor_idx)
        
        meiler_value = node_attr_keys["meiler"].get(meiler_idx)
        extracted_graph.nodes[node]["meiler"] = pd.Series(
            meiler_value[0], 
            name=''.join(list(meiler_value[1])), 
            index=[f'dim_{x}' for x in [1,2,3,4,5,6,7]]
        )
        
        

def _reconstruct_edge_attributes(extracted_graph, edges, edge_kinds, edge_kind_keys, edge_funcs, pdb_to_edges, pdb_code, edge_distances, edge_to_id, all_pdb_codes):
    pdb_idx = all_pdb_codes.index(pdb_code)

    for u, v in edges:
        edge = (u, v)
        
        if edge not in edge_kinds:
            continue
            
        kind_indexes_list = edge_kinds[edge]
        
        if len(kind_indexes_list) > 0:
            kind_indexes = kind_indexes_list 
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
        distance = edge_pair[1]
        extracted_graph.edges[edge]["distance"] = edge_pair[1]

        # if "bb_carbonyl_carbonyl" in extracted_graph.edges[edge]["kind"] and distance >=3.20:
        #     extracted_graph.edges[edge]["kind"].remove("bb_carbonyl_carbonyl")


def _reconstruct_and_validate_graphs(protein_graphs,
                                    node_to_id,
                                    edge_to_id,
                                    pdb_to_nodes,
                                    pdb_to_edges,
                                    node_attrs_global,
                                    edge_kinds,
                                    node_attr_keys,
                                    edge_kind_keys,
                                    edge_distances,
                                    all_pdb_codes,
                                    node_attrs_unique):

    for pdb_code in protein_graphs:
        pdb_graphs = protein_graphs[pdb_code]

        pdb_idx = all_pdb_codes.index(pdb_code)

        for graph in pdb_graphs:
            original_graph = graph.copy()

            def get_node_id(node_attr_unique_idx):
                node_pair = node_attrs_unique[node_attr_unique_idx]
                node_pair = struct.unpack("llll", node_pair)
                node_id = node_pair[0]
                return node_id

            nodes = [node_to_id.inverse[get_node_id(node_id)] for node_id in pdb_to_nodes[pdb_idx]]

            def get_edge_id(edge_distances_idx):
                edge_pair = edge_distances[edge_distances_idx]
                edge_pair = struct.unpack("ld", edge_pair)
                edge_id = edge_pair[0]
                return edge_id
            
            edges = [edge_to_id.inverse[get_edge_id(edge_distances_idx)] for edge_distances_idx in pdb_to_edges[pdb_idx]]

            extracted_graph = nx.Graph()
            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=nodes)

            _reconstruct_node_attributes(extracted_graph, node_attr_keys,node_attrs_global,node_attrs_unique,node_to_id,pdb_to_nodes, all_pdb_codes, pdb_code)
            _reconstruct_edge_attributes(extracted_graph, edges, edge_kinds, edge_kind_keys, original_graph.graph["config"].edge_construction_functions, pdb_to_edges, pdb_code, edge_distances, edge_to_id, all_pdb_codes)

            print(f"original {pdb_code} graph edge funcs: {original_graph.graph['config'].edge_construction_functions}")
            print(f"number of nodes in original graph: {len(original_graph.nodes())}")
            print(f"number of nodes in extracted graph: {len(extracted_graph.nodes())}")

            print(f"number of edges in original graph: {len(original_graph.edges())}")
            print(f"number of edges in extracted graph: {len(extracted_graph.edges())}")

            # for n in original_graph.nodes:
            #     print(f"original: {original_graph.nodes[n]}")
            #     print(f"extracted: {extracted_graph.nodes[n]}")

            # # diff = set(original_graph.edges) - set(extracted_graph.edges)
            # g1 = original_graph
            # g2 = extracted_graph
            # print(original_graph.graph["config"].edge_construction_functions)

            # for e in set(g1.edges).intersection(g2.edges):
            #     attrs1 = g1.edges[e]
            #     attrs2 = g2.edges[e]

            #     for k in set(attrs1.keys()).union(attrs2.keys()):
            #         if attrs1.get(k) != attrs2.get(k):
            #             print(f"({pdb_code})Edge {e} differs in '{k}': {attrs1.get(k)} != {attrs2.get(k)}")
            #             print(original_graph.edges[e])
            #             print(extracted_graph.edges[e])

def compress_pdb_graphs(protein_graphs):
    print("reached compress")
    edge_to_pdbs, \
    node_to_pdbs, \
    edge_kinds,   \
    node_attrs_global,\
    pdb_to_nodes, \
    pdb_to_edges, \
    pdb_codes_config, \
    all_pdb_codes, \
    node_attrs_unique_temp, \
    edge_distances_temp = _initialize_structs()

    global edge_functions_func_list
    global edge_functions_str_list

    edge_functions_func_list = [v for _, v in edgeModel.edge_functions_dict.items()]
    edge_functions_str_list = [k for k, _ in edgeModel.edge_functions_dict.items()]

    edge_kind_keys, node_attr_keys = _extract_attribute_keys()

    _process_nodes(protein_graphs, node_to_pdbs, node_attrs_global, node_attrs_unique_temp, node_attr_keys, pdb_to_nodes, all_pdb_codes)
    _process_edges(protein_graphs, edge_to_pdbs, edge_kinds, edge_kind_keys, pdb_to_edges, all_pdb_codes, edge_distances_temp)

    print(all_pdb_codes)

    node_to_id, edge_to_id, node_attrs_unique, edge_distances = _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes, all_pdb_codes, edge_distances_temp, node_attrs_unique_temp)

    _process_pdb_codes_config(protein_graphs, pdb_codes_config)

    del edge_to_pdbs
    del node_to_pdbs
    del edge_distances_temp
    del node_attrs_unique_temp

    _reconstruct_and_validate_graphs(protein_graphs,
                                    node_to_id,
                                    edge_to_id,
                                    pdb_to_nodes,
                                    pdb_to_edges,
                                    node_attrs_global,
                                    edge_kinds,
                                    node_attr_keys,
                                    edge_kind_keys,
                                    edge_distances,
                                    all_pdb_codes,
                                    node_attrs_unique)
    
    body_parts = {}

    body_parts["node_to_id"] = node_to_id
    body_parts["edge_to_id"] = edge_to_id
    body_parts["pdb_to_nodes"] = pdb_to_nodes
    body_parts["pdb_to_edges"] = pdb_to_edges
    body_parts["node_attrs_global"] = node_attrs_global
    body_parts["node_attr_keys"] = node_attr_keys
    body_parts["edge_kinds"] = edge_kinds
    body_parts["edge_kind_keys"] = edge_kind_keys
    body_parts["edge_distances"] = edge_distances
    body_parts["all_pdb_codes"] = all_pdb_codes
    body_parts["node_attrs_unique"] = node_attrs_unique

    return body_parts