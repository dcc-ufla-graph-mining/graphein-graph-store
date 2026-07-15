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
                "node_local_attr_keys": ["coords", "b_factor"],
                "node_attr_values": bidict(),
                "edge_attr_values": bidict(),
                "node_local_attr_keyvalue_mapping": {},
                "edge_local_attr_keyvalue_mapping": {}
            }

        self.granularity = ""

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_this_pdb_list())} pdbs'
    
    def __set_granularity(self, granularity: str):
        self.granularity = granularity

    def __get_granularity(self):
        return self.granularity
    
    def print_attr(self):
        for k, v in self.__body_parts.items():
            print(f"{k}: {v}")
            print("\n")
    
    def get_body_parts(self):
        return self.__body_parts

    def get_this_pdb_list(self):
        return self.__body_parts["pdb_code_to_id"].keys()
    
    def __edge_label_undirected(self, edge_label: tuple) -> tuple:
        return tuple(sorted(edge_label))

    def __get_residue_name_to_meiler_dict(self):
        '''
        retorna um dicionario mapeando um `residue_name` para seu respectivo `meiler`
        '''
        MEILER = {
            "ALA": (1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23),
            "ARG": (2.34, 0.29, 6.13, -1.01, 10.74, 0.36, 0.25),
            "ASN": (1.60, 0.13, 2.95, -0.60, 6.52, 0.21, 0.22),
            "ASP": (1.60, 0.11, 2.78, -0.77, 2.95, 0.25, 0.20),
            "CYS": (1.77, 0.13, 2.43, 1.54, 6.35, 0.17, 0.41),
            "GLN": (1.56, 0.18, 3.95, -0.22, 5.65, 0.35, 0.25),
            "GLU": (1.56, 0.15, 3.78, -0.64, 3.09, 0.42, 0.21),
            "GLY": (0.00, 0.00, 0.00, 0.00, 6.07, 0.13, 0.15),
            "HIS": (2.99, 0.23, 4.66, 0.13, 7.69, 0.27, 0.30),
            "ILE": (4.19, 0.19, 4.00, 1.80, 6.04, 0.30, 0.45),
            "LEU": (2.59, 0.19, 4.00, 1.70, 6.04, 0.39, 0.31),
            "LYS": (1.89, 0.22, 4.77, -0.99, 9.99, 0.32, 0.27),
            "MET": (2.35, 0.22, 4.43, 1.23, 5.71, 0.38, 0.32),
            "PHE": (2.94, 0.29, 5.89, 1.79, 5.67, 0.30, 0.38),
            "PRO": (2.67, 0.00, 2.72, 0.72, 6.80, 0.13, 0.34),
            "SER": (1.31, 0.06, 1.60, -0.04, 5.70, 0.20, 0.28),
            "THR": (3.03, 0.11, 2.60, 0.26, 5.60, 0.21, 0.36),
            "TRP": (3.21, 0.41, 8.08, 2.25, 5.94, 0.32, 0.42),
            "TYR": (2.94, 0.30, 6.47, 0.96, 5.66, 0.25, 0.41),
            "VAL": (3.67, 0.14, 3.00, 1.22, 6.02, 0.27, 0.49),
            "UNK": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        }

        return MEILER
    
    def __get_meiler_by_residue(self, residue_name: str):
        return self.__get_residue_name_to_meiler_dict()[residue_name]

    def __get_atom_type_to_element_symbol_dict(self):
        '''
        retorna um dicionario que mapeia todos os `atom_type` existentes para seus respectivos `element_symbol` 
        '''
        atom_types = [
            "N", "CA", "C", "O", "CB", "OG", "CG", "CD1", "CD2", "CE1", "CE2",
            "CZ", "OD1", "ND2", "CG1", "CG2", "CD", "CE", "NZ", "OD2", "OE1",
            "NE2", "OE2", "OH", "NE", "NH1", "NH2", "OG1", "SD", "ND1", "SG",
            "NE1", "CE3", "CZ2", "CZ3", "CH2", "OXT"
        ]

        return {atom: atom[0] for atom in atom_types}
    
    def __get_element_symbol_by_atom_type(self, atom_type: str):
        return self.__get_atom_type_to_element_symbol_dict()[atom_type]

    def insert(self, pdb_to_insert: dict):
        '''
        input: dict[str: nx.Graph]
        '''
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
            if not self.granularity:
                k = list(pdb_to_insert.keys())
                gr = pdb_to_insert[k[0]][0].graph['config'].dict()['granularity']
                self.__set_granularity(gr)

            for pdb_code, pdb_graph in pdb_to_insert.items():
                if pdb_code not in self.__body_parts["pdb_code_to_id"]:
                    self.__body_parts["pdb_code_to_id"][pdb_code] = len(self.__body_parts["pdb_code_to_id"])
                
                pdb_id = self.__body_parts["pdb_code_to_id"][pdb_code]
                if pdb_id not in self.__body_parts["pdb_id_to_edges"]:
                    self.__body_parts["pdb_id_to_edges"][pdb_id] = BitMap64()
                if pdb_id not in self.__body_parts["pdb_id_to_nodes"]:
                    self.__body_parts["pdb_id_to_nodes"][pdb_id] = BitMap64()

                __construct_structure_attributes(pdb_graph, pdb_id)

                __process_nodes(pdb_graph, pdb_id)
                __process_edges(pdb_graph, pdb_id)
        
        insert()

    def extract(self, pdb_to_extract: str) -> nx.Graph:        
        def __reconstruct_node_global_attrs(node_label: str, extracted_graph: nx.Graph):
            # global_attribute_keys = ["chain_id", "residue_name","residue_number", "atom_type", "element_symbol", "meiler"]
            # node_label = {chain_id} : {residue_name} : {residue_number} : {atom_type}

            attr_aux = node_label.split(":")
            chain_id = attr_aux[0]
            residue_name = attr_aux[1]
            residue_number = attr_aux[2]

            granularity = self.__get_granularity()

            if granularity == 'atom':
                atom_type = attr_aux[3]
            else:
                atom_type = granularity

            element_symbol = self.__get_element_symbol_by_atom_type(atom_type)
            meiler_tuple = self.__get_meiler_by_residue(residue_name)

            meiler = pd.Series(meiler_tuple,
                            name=residue_name,
                            index=[f'dim_{x}' for x in [1,2,3,4,5,6,7]])
                
            extracted_graph.nodes[node_label]['chain_id'] = chain_id
            extracted_graph.nodes[node_label]['residue_name'] = residue_name
            extracted_graph.nodes[node_label]['residue_number'] = residue_number
            extracted_graph.nodes[node_label]['atom_type'] = atom_type
            extracted_graph.nodes[node_label]['element_symbol'] = element_symbol
            extracted_graph.nodes[node_label]['meiler'] = meiler


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

                __reconstruct_node_global_attrs(node_label, extracted_graph)
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
