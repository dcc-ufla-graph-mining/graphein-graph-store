import networkx as nx
from pyroaring import BitMap64, BitMap
from bidict import bidict
import numpy as np
import pandas as pd

class PDBGraphStore:
    def __init__(self, body_parts):
        if body_parts:
            self.__body_parts = body_parts
        else:
            self.__body_parts = {
                "pdb_code_to_id": {},
                "pdb_id_to_nodes": {},
                "pdb_id_to_edges": {},
                "node_label_to_node_id": bidict(),
                "edge_label_to_edge_id": bidict(),
                "edge_attr_keys": ["kind", "distance"],
                "node_global_attr_keys": ["chain_id", "residue_name", "atom_type", "element_symbol", "meiler"],
                "node_local_attr_keys": ["residue_number", "coords", "b_factor"],
                "node_attr_values": bidict(),
                "edge_attr_values": bidict(),
                "node_global_attr_keyvalue_mapping": '',
                "node_local_attr_keyvalue_mapping": {},
                "edge_local_attr_keyvalue_mapping": {}
            }
        print(self.__body_parts.keys())

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_code_list())} pdbs'
    
    def get_body_parts(self):
        return self.__body_parts

    def get_this_pdb_list(self):
        return self.__body_parts["pdb_code_to_id"].keys()
    
    def __edge_label_undirected(self, edge_label: tuple) -> tuple:
        return tuple(sorted(edge_label))

    # TODO considerar transformar esse metodo em uma classe a parte, eliminando assim a necessidade do Builder.py tambem
    def insert_pdb(self, pdb_to_insert: dict):
        def __process_edge_distances(distance: float)-> list:
            distance_keyvalue_mapping_list = []

            attr_key = "distance"
            attr_value = distance

            if attr_value not in self.__body_parts["edge_attr_values"]:
                self.__body_parts["edge_attr_values"][attr_value] = len(self.__body_parts["edge_attr_values"])

            attr_key_id = self.__body_parts["edge_attr_keys"].index(attr_key)
            attr_value_id = self.__body_parts["edge_attr_values"][attr_value]

            distance_keyvalue_mapping_list.append(attr_key_id)
            distance_keyvalue_mapping_list.append(attr_value_id)

            return distance_keyvalue_mapping_list

        def __process_edge_kinds(kinds: set) -> list:
            kind_keyvalue_mapping_list = []

            attr_key = "kind"
            attr_key_id = self.__body_parts["edge_attr_keys"].index(attr_key)

            kind_keyvalue_mapping_list.append(attr_key_id)

            for kind in kinds:
                attr_value = kind

                if attr_value not in self.__body_parts["edge_attr_values"]:
                    self.__body_parts["edge_attr_values"][attr_value] = len(self.__body_parts["edge_attr_values"])

                attr_value_id = self.__body_parts["edge_attr_values"][attr_value]

                kind_keyvalue_mapping_list.append(attr_value_id)

            return kind_keyvalue_mapping_list

        def __process_edge_attrs(pdb_id: int, edge_id: int, edge: dict):
            local_attr_keyvalue_mapping = []

            local_attr_keyvalue_mapping.extend(__process_edge_distances(edge["distance"]))
            local_attr_keyvalue_mapping.extend(__process_edge_kinds(edge["kind"]))

            self.__body_parts["edge_local_attr_keyvalue_mapping"][(pdb_id, edge_id)] = local_attr_keyvalue_mapping

        def __process_edges(g: nx.Graph, pdb_id: int):
            for e in g.edges:
                edge_id = self.__body_parts["edge_label_to_edge_id"][self.__edge_label_undirected(e)]
                __process_edge_attrs(pdb_id, edge_id, g.edges[e])

        def __process_global_node_attrs(node: dict, node_id: int):
            for attr_idx, node_global_attr in enumerate(self.__body_parts["node_global_attr_keys"]):
                attr_value = node[node_global_attr]

                if isinstance(attr_value, pd.Series):
                    attr_value = tuple([tuple(attr_value.tolist()), tuple(attr_value.name)])
                
                if attr_value not in self.__body_parts["node_attr_values"]:
                    self.__body_parts["node_attr_values"][attr_value] = len(self.__body_parts["node_attr_values"])

                attr_value_id = self.__body_parts["node_attr_values"][attr_value]
                self.__body_parts["node_global_attr_keyvalue_mapping"][node_id][attr_idx] = attr_value_id

        def __process_local_node_attrs(node: dict) -> list:
            local_attr_list = []
            for node_local_attr in self.__body_parts["node_local_attr_keys"]:
                attr_value = node[node_local_attr]

                if isinstance(attr_value, np.ndarray):
                    attr_value = tuple(attr_value)
                
                if attr_value not in self.__body_parts["node_attr_values"]:
                    self.__body_parts["node_attr_values"][attr_value] = len(self.__body_parts["node_attr_values"])

                attr_value_id = self.__body_parts["node_attr_values"][attr_value]
                local_attr_list.append(attr_value_id)

            return local_attr_list

        def __process_node_attrs(pdb_id: int, node_id: int, node: dict):
            __process_global_node_attrs(node, node_id)
            local_attr_list = __process_local_node_attrs(node)

            self.__body_parts["node_local_attr_keyvalue_mapping"][(pdb_id, node_id)] = local_attr_list

        def __process_nodes(g: nx.Graph, pdb_id: int):
            for n in g.nodes:
                node_id = self.__body_parts["node_label_to_node_id"][n]
                __process_node_attrs(pdb_id, node_id, g.nodes[n])

        def __construct_node_structure(g: nx.graph, pdb_id: int):
            for node_label in g.nodes:
                if node_label not in self.__body_parts["node_label_to_node_id"]:
                    self.__body_parts["node_label_to_node_id"][node_label] = len(self.__body_parts["node_label_to_node_id"])

                node_id = self.__body_parts["node_label_to_node_id"][node_label]
                self.__body_parts["pdb_id_to_nodes"][pdb_id].add(node_id)

        def __construct_edge_structure(g: nx.Graph, pdb_id: int): 
            for e in g.edges:
                edge_label = self.__edge_label_undirected(e)

                if edge_label not in self.__body_parts["edge_label_to_edge_id"]:
                    self.__body_parts["edge_label_to_edge_id"][edge_label] = len(self.__body_parts["edge_label_to_edge_id"])
                
                edge_id = self.__body_parts["edge_label_to_edge_id"][edge_label]
                self.__body_parts["pdb_id_to_edges"][pdb_id].add(edge_id)

        def __construct_structure_attributes(g: nx.Graph, pdb_id: int):
            __construct_node_structure(g, pdb_id)
            __construct_edge_structure(g, pdb_id)

        def insert():
            for pdb_code, pdb_graph in pdb_to_insert.items():
                if pdb_code not in self.__body_parts["pdb_code_to_id"]:
                    self.__body_parts["pdb_code_to_id"][pdb_code] = len(self.__body_parts["pdb_code_to_id"])
                
                pdb_id = self.__body_parts["pdb_code_to_id"][pdb_code]
                if pdb_id not in self.__body_parts["pdb_id_to_edges"]:
                    self.__body_parts["pdb_id_to_edges"][pdb_id] = BitMap64()
                if pdb_id not in self.__body_parts["pdb_id_to_nodes"]:
                    self.__body_parts["pdb_id_to_nodes"][pdb_id] = BitMap64()

                __construct_structure_attributes(pdb_graph[0], pdb_id)

                old = self.__body_parts["node_global_attr_keyvalue_mapping"]

                mshape = 5
                mtype = np.int32

                new_size = len(self.__body_parts["node_label_to_node_id"])
                new = np.zeros((new_size, mshape), dtype=mtype)

                if (type(old) != str):
                    new[:old.shape[0]] = old

                self.__body_parts["node_global_attr_keyvalue_mapping"] = new

                __process_nodes(pdb_graph[0], pdb_id)
                __process_edges(pdb_graph[0], pdb_id)
        
        insert()

    # TODO: considerar transformar esse método em uma classe a parte
    def extract_pdb(self, pdb_to_extract: str) -> nx.Graph:        
        def __reconstruct_node_global_attrs(node_id: int, extracted_graph: nx.Graph):
            global_attributes = self.__body_parts["node_global_attr_keyvalue_mapping"][node_id]
            global_attribute_keys = self.__body_parts["node_global_attr_keys"]

            node_label = self.__body_parts["node_label_to_node_id"].inverse[node_id]

            for attr_idx, attr_key in enumerate(global_attribute_keys):
                attr_value_id = global_attributes[attr_idx]
                attr_value = self.__body_parts["node_attr_values"].inverse[attr_value_id]

                if isinstance(attr_value, tuple):
                    attr_value = pd.Series(attr_value[0],
                                        name=''.join(list(attr_value[1])),
                                        index=[f'dim_{x}' for x in [1,2,3,4,5,6,7]]                                           
                                    )
                    
                extracted_graph.nodes[node_label][attr_key] = attr_value

        def __reconstruct_node_local_attrs(node_id: int, pdb_id: int, extracted_graph):
            local_attributes = self.__body_parts["node_local_attr_keyvalue_mapping"][(pdb_id, node_id)]
            local_attributes_keys = self.__body_parts["node_local_attr_keys"]

            node_label = self.__body_parts["node_label_to_node_id"].inverse[node_id]

            for attr_idx, attr_key in enumerate(local_attributes_keys):
                attr_value_id = local_attributes[attr_idx]
                attr_value = self.__body_parts["node_attr_values"].inverse[attr_value_id]

                if isinstance(attr_value, tuple):
                    attr_value = np.array(attr_value)

                extracted_graph.nodes[node_label][attr_key] = attr_value
        
        def __reconstruct_edge_kinds(attributes: list, g: nx.Graph, edge_label: str):
            attr_key = "kind"
            kind_list = []

            for i in range(1, len(attributes)):
                attr_value_id = attributes[i]
                kind_value = self.__body_parts["edge_attr_values"].inverse[attr_value_id]
                kind_list.append(kind_value)

            kinds = set(kind_list)

            g.edges[edge_label][attr_key] = kinds

        def __reconstruct_edge_distance(attributes: list, g: nx.Graph, edge_label: str):
            attr_key = "distance"
            attr_value_id = attributes[1]

            distance_value = self.__body_parts["edge_attr_values"].inverse[attr_value_id]

            g.edges[edge_label][attr_key] = distance_value

        def __reconstruct_nodes(extracted_graph: nx.Graph, pdb_id: int):
            for node_label in extracted_graph.nodes:
                node_id = self.__body_parts["node_label_to_node_id"][node_label]

                __reconstruct_node_global_attrs(node_id, extracted_graph)
                __reconstruct_node_local_attrs(node_id, pdb_id, extracted_graph)

        def __reconstruct_edges(extracted_graph: nx.Graph, pdb_id: int):
            for edge_label in extracted_graph.edges:
                edge_id = self.__body_parts["edge_label_to_edge_id"][self.__edge_label_undirected(edge_label)]
                attributes = self.__body_parts["edge_local_attr_keyvalue_mapping"][(pdb_id, edge_id)]

                __reconstruct_edge_kinds(attributes[2:], extracted_graph, edge_label)
                __reconstruct_edge_distance(attributes[:2], extracted_graph, edge_label)

        extracted_graph = nx.Graph()
        pdb_id = self.__body_parts["pdb_code_to_id"][pdb_to_extract]

        nodes = [self.__body_parts["node_label_to_node_id"].inverse[node_id] for node_id in self.__body_parts["pdb_id_to_nodes"][pdb_id]]
        edges = [self.__body_parts["edge_label_to_edge_id"].inverse[edge_id] for edge_id in self.__body_parts["pdb_id_to_edges"][pdb_id]]

        extracted_graph.add_nodes_from(nodes)
        extracted_graph.add_edges_from(edges)

        __reconstruct_nodes(extracted_graph, pdb_id)
        __reconstruct_edges(extracted_graph, pdb_id)

        return extracted_graph
