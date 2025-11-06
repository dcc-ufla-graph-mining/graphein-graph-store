import networkx as nx
import numpy as np
from bidict import bidict
from pyroaring import BitMap64
import pandas as pd
import edge_functions_Model as edgeModel

from Indexed_set import IndexedSet

edge_functions_func_list = []
edge_functions_str_list = []

def _initialize_structs():
    node_attrs_global = {}
    pdb_to_nodes = {}
    pdb_to_edges = {}
    pdb_codes_config = {}
    all_pdb_codes = IndexedSet()
    
    node_to_id = bidict()
    edge_to_id = bidict()
    
    node_attrs_unique = []
    edge_attrs = []

    return (node_attrs_global, pdb_to_nodes, pdb_to_edges, pdb_codes_config, 
            all_pdb_codes, node_to_id, edge_to_id, node_attrs_unique, edge_attrs)

def _extract_attribute_keys():
    node_attr_keys_list = ["chain_id", "residue_name", "residue_number", 
                           "atom_type", "element_symbol", "coords", "b_factor", "meiler"]

    node_attr_keys = {}
    for key in node_attr_keys_list:
        node_attr_keys[key] = IndexedSet()

    edge_kind_keys = IndexedSet(set(edge_functions_str_list))

    return edge_kind_keys, node_attr_keys

def _process_node_attributes(node, graph, node_attr_keys):
    attr_indexes_global = []
    attr_indexes_unique = []

    for i, value in enumerate(node_attr_keys):
        attr_value = graph.nodes[node][value]

        if i == 7:  # attr: meiler; type: pd.Series
            attr_value = tuple([tuple(attr_value.tolist()), tuple(attr_value.name)])
        
        if i == 5:  # attr: coords; type: np.array
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
        if kind:
            kind_indexes.add(edge_kind_keys.index(kind))

    return kind_indexes

def _process_nodes(protein_graphs, node_attrs_global, node_attr_keys, 
                   pdb_to_nodes, all_pdb_codes, node_to_id, node_attrs_unique):
    for pdb_code, graphs in protein_graphs.items():
        all_pdb_codes.add(pdb_code)
        pdb_idx = all_pdb_codes.index(pdb_code)
        
        if pdb_idx not in pdb_to_nodes:
            pdb_to_nodes[pdb_idx] = BitMap64()

        for g in graphs:
            for node in g.nodes():
                attr_indexes_global, attr_indexes_unique = _process_node_attributes(node, g, node_attr_keys)
                
                # Map node to ID if not already mapped
                if node not in node_to_id:
                    node_id = len(node_to_id)
                    node_to_id[node] = node_id
                    node_attrs_global[node] = attr_indexes_global
                else:
                    node_id = node_to_id[node]
                
                # Create node attribute entry
                node_pair = tuple([node_id, attr_indexes_unique[0], 
                                 attr_indexes_unique[1], attr_indexes_unique[2]])
                node_attrs_unique.append(node_pair)
                pdb_to_nodes[pdb_idx].add(len(node_attrs_unique) - 1)

def _process_edges(protein_graphs, edge_kind_keys, pdb_to_edges, 
                   all_pdb_codes, edge_to_id, edge_attrs):
    for pdb_code, graphs in protein_graphs.items():
        pdb_idx = all_pdb_codes.index(pdb_code)
        
        if pdb_idx not in pdb_to_edges:
            pdb_to_edges[pdb_idx] = BitMap64()
            
        for g in graphs:
            for u, v, data in g.edges.data():
                edge = (u, v)
                
                # Map edge to ID if not already mapped
                if edge not in edge_to_id:
                    edge_id = len(edge_to_id)
                    edge_to_id[edge] = edge_id
                else:
                    edge_id = edge_to_id[edge]
                
                try:
                    kind_indexes = _process_edge_kinds(data, edge_kind_keys)
                except ValueError as e:
                    print("error at process edges: ", e)
                    continue
                
                # Store edge distance entry: (edge_id, distance, kind_indexes)
                edge_pair = tuple([edge_id, data["distance"], kind_indexes])
                edge_attrs.append(edge_pair)
                pdb_to_edges[pdb_idx].add(len(edge_attrs) - 1)

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

def _reconstruct_node_attributes(extracted_graph, node_attr_keys, node_attrs_global, 
                                node_attrs_unique, node_to_id, pdb_to_nodes, 
                                all_pdb_codes, pdb_code):
    pdb_idx = all_pdb_codes.index(pdb_code)
    
    for node_attr_idx in pdb_to_nodes[pdb_idx]:
        node_pair = node_attrs_unique[node_attr_idx]
        node_id = node_pair[0]
        node = node_to_id.inverse[node_id]
        
        chain_id_idx = node_attrs_global[node][0]
        residue_name_idx = node_attrs_global[node][1]
        atom_type_idx = node_attrs_global[node][2]
        element_symbol_idx = node_attrs_global[node][3]
        meiler_idx = node_attrs_global[node][4]

        residue_number_idx = node_pair[1]
        coords_idx = node_pair[2]
        b_factor_idx = node_pair[3]
        
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

def _reconstruct_edge_attributes(extracted_graph, edges, edge_kind_keys, 
                                pdb_to_edges, pdb_code, edge_attrs, 
                                edge_to_id, all_pdb_codes):
    pdb_idx = all_pdb_codes.index(pdb_code)

    for edge_attr_idx in pdb_to_edges[pdb_idx]:
        edge_pair = edge_attrs[edge_attr_idx]

        edge = edge_to_id.inverse[edge_pair[0]]
        if edge not in edges:
            continue

        if not extracted_graph.has_edge(*edge):
            extracted_graph.add_edge(*edge)

        distance = edge_pair[1]
        kinds = edge_pair[2]

        kinds = [edge_kind_keys.get(k) for k in kinds]

        extracted_graph.edges[edge].setdefault("kind", set())

        for k in kinds:
            extracted_graph.edges[edge]["kind"].add(k)
            
        extracted_graph.edges[edge]["distance"] = distance

def _reconstruct_and_validate_graphs(protein_graphs, node_to_id, edge_to_id,
                                    pdb_to_nodes, pdb_to_edges, node_attrs_global,
                                    node_attr_keys, edge_kind_keys, edge_attrs,
                                    all_pdb_codes, node_attrs_unique):

    for pdb_code in protein_graphs:
        pdb_graphs = protein_graphs[pdb_code]
        pdb_idx = all_pdb_codes.index(pdb_code)

        for graph in pdb_graphs:
            original_graph = graph.copy()

            def get_node_id(node_attr_unique_idx):
                node_pair = node_attrs_unique[node_attr_unique_idx]
                return node_pair[0]

            nodes = [node_to_id.inverse[get_node_id(node_id)] 
                    for node_id in pdb_to_nodes[pdb_idx]]

            def get_edge_id(edge_attrs_idx):
                edge_pair = edge_attrs[edge_attrs_idx]
                return edge_pair[0]
            
            edges = [edge_to_id.inverse[get_edge_id(edge_attrs_idx)] 
                    for edge_attrs_idx in pdb_to_edges[pdb_idx]]

            extracted_graph = nx.Graph()
            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=nodes)

            _reconstruct_node_attributes(extracted_graph, node_attr_keys, node_attrs_global,
                                       node_attrs_unique, node_to_id, pdb_to_nodes, 
                                       all_pdb_codes, pdb_code)
            _reconstruct_edge_attributes(extracted_graph, edges, edge_kind_keys,
                                       pdb_to_edges, pdb_code, edge_attrs, 
                                       edge_to_id, all_pdb_codes)

            print(f"original {pdb_code} graph edge funcs: {original_graph.graph['config'].edge_construction_functions}")
            print(f"number of nodes in original graph: {len(original_graph.nodes())}")
            print(f"number of nodes in extracted graph: {len(extracted_graph.nodes())}")
            print(f"number of edges in original graph: {len(original_graph.edges())}")
            print(f"number of edges in extracted graph: {len(extracted_graph.edges())}")

            # Validation
            def canonical_edges(G):
                return {(min(u, v), max(u, v)) for u, v in G.edges}
            
            msg = f'\n\n{original_graph.graph["config"].edge_construction_functions}\n\n'
            
            for e in set(canonical_edges(original_graph)) - set(canonical_edges(extracted_graph)):
                msg += f"\n{e}\n"
                msg += f"original: {original_graph.edges[e]}\n"
                msg += f"extracted: {extracted_graph.edges.get(e, 'NOT FOUND')}\n"
                
            for e in set(canonical_edges(extracted_graph)) - set(canonical_edges(original_graph)):
                msg += f"\n{e}\n"
                msg += f"original: {original_graph.edges.get(e, 'NOT FOUND')}\n"
                msg += f"extracted: {extracted_graph.edges[e]}\n"
                
            for e in set(original_graph.edges) & set(extracted_graph.edges):
                if original_graph.edges[e] != extracted_graph.edges[e]:
                    msg += f"\nDifferent attributes in edge {e}:\n"
                    msg += f"original:  {original_graph.edges[e]}\n"
                    msg += f"extracted: {extracted_graph.edges[e]}\n"

            print(msg)

def compress_pdb_graphs(protein_graphs):
    print("reached compress")
    
    (node_attrs_global, pdb_to_nodes, pdb_to_edges, pdb_codes_config,
     all_pdb_codes, node_to_id, edge_to_id, node_attrs_unique, 
     edge_attrs) = _initialize_structs()

    global edge_functions_func_list
    global edge_functions_str_list

    edge_functions_func_list = [v for _, v in edgeModel.edge_functions_dict.items()]
    edge_functions_str_list = [k for k, _ in edgeModel.edge_functions_dict.items()]

    edge_kind_keys, node_attr_keys = _extract_attribute_keys()

    _process_nodes(protein_graphs, node_attrs_global, node_attr_keys, 
                   pdb_to_nodes, all_pdb_codes, node_to_id, node_attrs_unique)
    _process_edges(protein_graphs, edge_kind_keys, pdb_to_edges, 
                   all_pdb_codes, edge_to_id, edge_attrs)

    print(all_pdb_codes)

    _process_pdb_codes_config(protein_graphs, pdb_codes_config)

    # _reconstruct_and_validate_graphs(protein_graphs, node_to_id, edge_to_id,
    #                                 pdb_to_nodes, pdb_to_edges, node_attrs_global,
    #                                 node_attr_keys, edge_kind_keys, edge_attrs,
    #                                 all_pdb_codes, node_attrs_unique)

    edge_distances_file = "edge_distances.txt"

    for edge_pair in edge_attrs:
        aux = edge_pair

        distance = edge_pair[1]

        with open(edge_distances_file, "a") as f:
            f.write(str(distance))
    
    body_parts = {
        "node_to_id": node_to_id,
        "edge_to_id": edge_to_id,
        "pdb_to_nodes": pdb_to_nodes,
        "pdb_to_edges": pdb_to_edges,
        "node_attrs_global": node_attrs_global,
        "node_attr_keys": node_attr_keys,
        "edge_kind_keys": edge_kind_keys,
        "edge_attrs": edge_attrs,
        "all_pdb_codes": all_pdb_codes,
        "node_attrs_unique": node_attrs_unique
    }

    return body_parts