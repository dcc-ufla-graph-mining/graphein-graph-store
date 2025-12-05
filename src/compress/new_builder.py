import networkx as nx
from pyroaring import BitMap64
from bidict import bidict
import numpy as np
import pandas as pd
import pickle as pk
import time
import torch

def initialize_body_parts():
    body_parts = {
        "pdb_code_to_id": {},
        "node_label_to_node_id": bidict(),
        "edge_label_to_edge_id": bidict(),
        "pdb_id_to_nodes": {},
        "pdb_id_to_edges": {},
        "attr_keys": 
                [
                      "empty", "chain_id", "residue_name", "residue_number", 
                      "atom_type", "element_symbol", "coords", "b_factor", 
                      "meiler", "kind", "distance"
                ],
        "attr_values": bidict(),
        "node_global_attr_keyvalue_mapping": '',
        "node_local_attr_keyvalue_mapping": {},
        "edge_local_attr_keyvalue_mapping": {}
    }

    body_parts["attr_values"]["empty"] = 0

    return body_parts

def edge_label_undirected(edge_label: tuple) -> tuple:
    return tuple(sorted(edge_label))

def process_edge_distances(distance: float, body_parts: dict) -> list:
    distance_keyvalue_mapping_list = []

    attr_key = "distance"
    attr_value = distance

    if attr_value not in body_parts["attr_values"]:
        body_parts["attr_values"][attr_value] = len(body_parts["attr_values"])

    attr_key_id = body_parts["attr_keys"].index(attr_key)
    attr_value_id = body_parts["attr_values"][attr_value]  

    distance_keyvalue_mapping_list.append(attr_key_id)
    distance_keyvalue_mapping_list.append(attr_value_id)

    return distance_keyvalue_mapping_list

def process_edge_kinds(kinds: set, body_parts: dict) -> list:
    kind_keyvalue_mapping_list = []

    attr_key = "kind"
    attr_key_id = body_parts["attr_keys"].index(attr_key)

    kind_keyvalue_mapping_list.append(attr_key_id)

    for kind in kinds:
        attr_value = kind

        if attr_value not in body_parts["attr_values"]:
            body_parts["attr_values"][attr_value] = len(body_parts["attr_values"])
        
        attr_value_id = body_parts["attr_values"][attr_value]

        kind_keyvalue_mapping_list.append(attr_value_id)

    return kind_keyvalue_mapping_list

def process_edge_attrs(body_parts: dict, pdb_id: int, edge_id: int, edge: dict):
    local_attr_keyvalue_mapping = []

    local_attr_keyvalue_mapping.extend(process_edge_distances(edge["distance"],  body_parts))
    local_attr_keyvalue_mapping.extend(process_edge_kinds(edge["kind"], body_parts))
        
    body_parts["edge_local_attr_keyvalue_mapping"][(pdb_id, edge_id)] = local_attr_keyvalue_mapping


def process_edges(g: nx.Graph, body_parts: dict, pdb_id: int):
    for e in g.edges:
        edge_id = body_parts["edge_label_to_edge_id"][edge_label_undirected(e)]
        process_edge_attrs(body_parts, pdb_id, edge_id, g.edges[e])

def process_node_attrs(body_parts: dict, pdb_id: int, node_id: int, node: dict):
    global_idx = 0

    global_node_attributes= ["chain_id", "residue_name", "atom_type", "element_symbol", "meiler"]
    local_node_attributes = ["residue_number", "coords", "b_factor"]

    local_attr_list = []

    for attr_key, attr_value in node.items():

        if isinstance(attr_value, pd.Series):
            attr_value = tuple([tuple(attr_value.tolist()), tuple(attr_value.name)])
        
        if isinstance(attr_value, np.ndarray):
            attr_value = tuple(attr_value)

        if attr_value not in body_parts["attr_values"]:
            body_parts["attr_values"][attr_value] = len(body_parts["attr_values"])

        attr_key_id = body_parts["attr_keys"].index(attr_key)
        attr_value_id = body_parts["attr_values"][attr_value]
        
        if attr_key in local_node_attributes:
            local_attr_list.append(attr_key_id)
            local_attr_list.append(attr_value_id)

        if attr_key in global_node_attributes:
            body_parts["node_global_attr_keyvalue_mapping"][node_id][global_idx] = attr_key_id
            body_parts["node_global_attr_keyvalue_mapping"][node_id][global_idx+1] = attr_value_id
            global_idx+=2

    body_parts["node_local_attr_keyvalue_mapping"][(pdb_id, node_id)] = local_attr_list

def process_nodes(g: nx.Graph, body_parts: dict, pdb_id: int):
    for n in g.nodes:
        node_id = body_parts["node_label_to_node_id"][n]
        process_node_attrs(body_parts, pdb_id, node_id, g.nodes[n])

def construct_node_structure(g: nx.Graph, body_parts: dict, pdb_id: int):
    for node_label in g.nodes:
        if node_label not in body_parts["node_label_to_node_id"]:
            body_parts["node_label_to_node_id"][node_label] = len(body_parts["node_label_to_node_id"])

        node_id = body_parts["node_label_to_node_id"][node_label]
        body_parts["pdb_id_to_nodes"][pdb_id].add(node_id)

def construct_edge_structure(g: nx.Graph, body_parts: dict, pdb_id: int):
    for e in g.edges:
        edge_label = edge_label_undirected(e)

        if edge_label not in body_parts["edge_label_to_edge_id"]:
            body_parts["edge_label_to_edge_id"][edge_label] = len(body_parts["edge_label_to_edge_id"])
        
        edge_id = body_parts["edge_label_to_edge_id"][edge_label]
        body_parts["pdb_id_to_edges"][pdb_id].add(edge_id)

def construct_structure_attributes(g: nx.Graph, body_parts: dict, pdb_id: int):
    construct_node_structure(g, body_parts, pdb_id)
    construct_edge_structure(g, body_parts, pdb_id)

def initialize_structures(protein_graphs: dict[str, list[nx.Graph]], body_parts: dict):
    for pdb_code, pdb_graph_list in protein_graphs.items():
        if pdb_code not in body_parts["pdb_code_to_id"]:
            body_parts["pdb_code_to_id"][pdb_code] = len(body_parts["pdb_code_to_id"])

        pdb_id = body_parts["pdb_code_to_id"][pdb_code]

        if pdb_id not in body_parts["pdb_id_to_nodes"]:
            body_parts["pdb_id_to_nodes"][pdb_id] = BitMap64()

        if pdb_id not in body_parts["pdb_id_to_edges"]:
            body_parts["pdb_id_to_edges"][pdb_id] = BitMap64()

        for g in pdb_graph_list:
            construct_structure_attributes(g, body_parts, pdb_id)

    body_parts["node_global_attr_keyvalue_mapping"] = np.zeros((
                                                                len(body_parts["node_label_to_node_id"]), 
                                                                10
                                                                ), dtype=np.int32)

def reconstruct_nodes(body_parts: dict, g: nx.Graph, pdb_id: int):
    for node_label in g.nodes:
        node_id = body_parts["node_label_to_node_id"][node_label]

        global_attributes = body_parts["node_global_attr_keyvalue_mapping"][node_id]
        local_attributes = body_parts["node_local_attr_keyvalue_mapping"][(pdb_id, node_id)]

        for i in range(0, len(global_attributes), 2):
            attr_key_id = global_attributes[i]
            attr_value_id = global_attributes[i+1]

            attr_value = body_parts["attr_values"].inverse[attr_value_id]
            attr_key = body_parts["attr_keys"][attr_key_id]
            
            if isinstance(attr_value, tuple):
                attr_value = pd.Series(
                    attr_value[0],
                    name=''.join(list(attr_value[1])),
                    index=[f'dim_{x}' for x in [1,2,3,4,5,6,7]]
                )

            g.nodes[node_label][attr_key] = attr_value

        for i in range(0, len(local_attributes), 2):
            attr_key_id = local_attributes[i]
            attr_value_id = local_attributes[i+1]

            attr_value = body_parts["attr_values"].inverse[attr_value_id]
            attr_key = body_parts["attr_keys"][attr_key_id]

            if isinstance(attr_value, tuple):
                attr_value = np.array(attr_value)

            g.nodes[node_label][attr_key] = attr_value

def reconstruct_edge_distance(attributes: list, g: nx.Graph, body_parts: dict, edge_label):
    attr_key = "distance"
    attr_value_id = attributes[1]
    distance_value = body_parts["attr_values"].inverse[attr_value_id]

    g.edges[edge_label][attr_key] = distance_value

def reconstruct_edge_kinds(attributes: list, g: nx.Graph, body_parts: dict, edge_label):
    attr_key = "kind"

    kind_list = []

    for i in range(1, len(attributes)):
        attr_value_id = attributes[i]
        kind_value = body_parts["attr_values"].inverse[attr_value_id]

        kind_list.append(kind_value)
    
    kinds = set(kind_list)

    g.edges[edge_label][attr_key] = kinds

def reconstruct_edges(body_parts: dict, g: nx.Graph, pdb_id: int):
    for edge_label in g.edges:
        edge_id = body_parts["edge_label_to_edge_id"][edge_label_undirected(edge_label)]

        attributes = body_parts["edge_local_attr_keyvalue_mapping"][(pdb_id, edge_id)]
        
        reconstruct_edge_kinds(attributes[2:], g, body_parts, edge_label)
        reconstruct_edge_distance(attributes[:2], g, body_parts, edge_label)
        

def reconstruct_and_validate(protein_graphs: dict[str, list[nx.Graph]], body_parts: dict):
    for pdb_code, pdb_graph_list in protein_graphs.items():
        for original_graph in pdb_graph_list:
            extracted_graph = nx.Graph()
            pdb_id = body_parts["pdb_code_to_id"][pdb_code]

            nodes = [body_parts["node_label_to_node_id"].inverse[node_id] for node_id in body_parts["pdb_id_to_nodes"][pdb_id]]
            edges = [body_parts["edge_label_to_edge_id"].inverse[edge_id] for edge_id in body_parts["pdb_id_to_edges"][pdb_id]]

            extracted_graph.add_nodes_from(nodes)
            extracted_graph.add_edges_from(edges)

            reconstruct_nodes(body_parts, extracted_graph, pdb_id)
            reconstruct_edges(body_parts, extracted_graph, pdb_id)

            # print(f"\n Nodes in original graph: {len(original_graph.nodes)}")
            # print(f"\n Nodes in extracted graph: {len(extracted_graph.nodes)}")
            # print(f"\n Edges in original graph: {len(original_graph.edges)}")
            # print(f"\n Edges in extracted graph: {len(extracted_graph.edges)}")

            nodes_equal = nx.utils.nodes_equal(original_graph.nodes, extracted_graph.nodes)
            edges_equal = nx.utils.edges_equal(original_graph.edges, extracted_graph.edges)
            
            if not nodes_equal:
                print("nodes are not equal")
            
            if not edges_equal:
                print("edges are not equal")

            for u, v in set(original_graph.edges) & set(extracted_graph.edges):
                try:
                    e = (u, v)
                    original_edge = original_graph.edges[e]
                    print(original_edge)
                except:
                    e = (v, u)
                    original_edge = original_graph.edges[e]
                    print(original_edge)
                
                try:
                    e = (u, v)
                    extracted_edge = extracted_graph.edges[e]
                    print(extracted_edge)
                except:
                    e = (v, u)
                    extracted_edge = extracted_graph.edges[e]
                    print(extracted_edge)

                if original_edge != extracted_edge:
                    print(f"\ndifferent attributes in edges {e}: \n")
                    print(f"original: {original_edge}")
                    print(f"extracted: {extracted_edge}")

            for n in set(original_graph.nodes) & set(extracted_graph.nodes):
                original_node = original_graph.nodes[n]
                extracted_node = extracted_graph.nodes[n]

                for k, v in original_node.items():
                    w = extracted_node[k]
                    if isinstance(v, pd.Series):
                        original_node[k] = tuple([tuple(v.to_list()), tuple(v.name)])
                        extracted_node[k] = tuple([tuple(w.to_list()), tuple(w.name)])

                    elif isinstance(v, np.ndarray):
                        original_node[k] = tuple([tuple(v)])
                        extracted_node[k] = tuple([tuple(w)])

                if original_node != extracted_node:
                    print(f"\ndifferent attributes in node {n}: \n")
                    print(f"original: {original_node}")
                    print(f"extracted: {extracted_node}")

def process_graphs(protein_graphs: dict[str, list[nx.Graph]], body_parts: dict):
    for pdb_code, pdb_graph_list in protein_graphs.items():
        for g in pdb_graph_list:
            pdb_id = body_parts["pdb_code_to_id"][pdb_code]

            process_nodes(g, body_parts, pdb_id)
            process_edges(g, body_parts, pdb_id)

def compress_pdb_graphs(protein_graphs: dict[str, list[nx.Graph]]) -> dict:
    time_to_construct_start = time.time()

    body_parts = initialize_body_parts()
    initialize_structures(protein_graphs, body_parts)
    process_graphs(protein_graphs, body_parts) 

    time_to_construct = time.time() - time_to_construct_start

    reconstruct_and_validate(protein_graphs, body_parts)

    print(f'number of edges: {len(body_parts["edge_label_to_edge_id"])}')
    
    return body_parts, time_to_construct