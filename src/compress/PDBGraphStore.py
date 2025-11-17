import networkx as nx
import numpy as np
from pyroaring import BitMap64
from bidict import bidict
import pandas as pd
import pickle as pk

class PDBGraphStore:
    def __init__(self, body_parts):
        if body_parts:
            self.__body_parts = body_parts
            print(self.__body_parts.keys())
        else:
            raise ValueError("'body_parts' should not be None")

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_code_list())} pdbs'


    def node_to_id_size(self):
        return asizeof.asizeof(self.__body_parts["node_to_id"]) / 1024 / 1024

    def edge_to_id_size(self):
        return asizeof.asizeof(self.__body_parts["edge_to_id"]) / 1024 / 1024

    def pdb_to_nodes_size(self):
        return asizeof.asizeof(self.__body_parts["pdb_to_nodes"]) / 1024 / 1024

    def pdb_to_edges_size(self):
        return asizeof.asizeof(self.__body_parts["pdb_to_edges"]) / 1024 / 1024

    def node_attrs_size(self):
        return asizeof.asizeof(self.__body_parts["node_attrs_global"]) / 1024 / 1024 + \
                asizeof.asizeof(self.__body_parts["node_attrs_unique"]) / 1024 / 1024

    def edge_attrs_size(self):
        return asizeof.asizeof(self.__body_parts["edge_attrs"]) / 1024 / 1024

    def node_attr_keys_size(self):
        return asizeof.asizeof(self.__body_parts["node_attr_keys"]) / 1024 / 1024

    def edge_kind_keys_size(self):
        return asizeof.asizeof(self.__body_parts["edge_kind_keys"]) / 1024 / 1024
    
    #TODO consertar essa mediçao
    def compressible_edge_parts_size(self):
        return (asizeof.asizeof(self.__body_parts["edge_kind_keys"])
        + asizeof.asizeof(self.__body_parts["edge_attrs"]) 
        )/ 1024 / 1024

    def compressible_node_parts_size(self):
        return ((asizeof.asizeof(self.__body_parts["node_attrs_global"]) 
                +asizeof.asizeof(self.__body_parts["node_attr_keys"]["chain_id"])
                +asizeof.asizeof(self.__body_parts["node_attr_keys"]["residue_name"])
                +asizeof.asizeof(self.__body_parts["node_attr_keys"]["atom_type"])
                +asizeof.asizeof(self.__body_parts["node_attr_keys"]["element_symbol"])
                )/1024 / 1024)    
    
    def incompressible_node_parts_size(self):
        return ((asizeof.asizeof(self.__body_parts["node_attrs_unique"]) 
                 +asizeof.asizeof(self.__body_parts["node_attr_keys"]["residue_name"])
                 +asizeof.asizeof(self.__body_parts["node_attr_keys"]["coords"])
                 +asizeof.asizeof(self.__body_parts["node_attr_keys"]["b_factor"])
                 +asizeof.asizeof(self.__body_parts["node_attr_keys"]["meiler"])
                 )/ 1024 / 1024)

    def calculate_graph_complete_space_size(self):
        return (
            self.node_to_id_size()
            + self.edge_to_id_size()
            + self.pdb_to_nodes_size()
            + self.pdb_to_edges_size()
            + self.node_attrs_size()
            + self.edge_attrs_size()
            + self.node_attr_keys_size()
            + self.edge_kind_keys_size()
        )

    def calculate_total_nodes_size(self):
        return (
            self.node_to_id_size()
            + self.pdb_to_nodes_size()
            + self.node_attrs_size()
            + self.node_attr_keys_size()
        )

    def calculate_total_edges_size(self):
        return (
            self.edge_to_id_size()
            + self.pdb_to_edges_size()
            + self.edge_attrs_size()
            + self.edge_kind_keys_size()
        )

