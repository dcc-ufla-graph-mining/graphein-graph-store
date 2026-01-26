import networkx as nx
from pyroaring import BitMap64
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

    def insert_pdb(self, pdb_to_insert: str):
        def __edge_label_undirected(edge_label: tuple) -> tuple:
            pass

        def __process_edge_distances(distance: float)-> list:
            pass

        def __process_edge_kinds(kinds: set) -> list:
            pass

        def __process_edge_attrs(pdb_id: int, edge_id: int, edge: dict):
            pass

        def __process_edges(g: nx.Graph, pdb_id: int):
            pass

        def __process_global_node_attrs(node: dict, node_id: int):
            pass

        def __process_local_node_attrs(node: dict) -> list:
            pass

        def __process_node_attrs(pdb_id: int, node_id: int, node: dict):
            pass

        def __process_nodes(g: nx.Graph, pdb_id: int):
            pass

        def __construct_node_structure(g: nx.graph, pdb_id: int):
            pass

        def __construct_edge_structure(g: nx.Graph, pdb_id: int): 
            pass

        def __construct_structure_attributes(g: nx.Graph, pdb_id):
            pass
    
    # return pdb
    def extract_pdb(self, pdb_to_extract: str) -> nx.Graph:
        pass

    # return str(pdb) se o pdb foi removido com sucesso 
    def remove_pdb(self, pdb_to_remove: str) -> str:
        pass

    # return list(str(pdb)) se o pdb foi removido com sucejsso
    def remove_multi_pdb(self, pdb_to_remove_list: list) -> list:
        pass
