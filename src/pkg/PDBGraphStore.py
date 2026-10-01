from __future__ import annotations
from typing import Optional, Union, Dict, Any, Tuple
import networkx as nx
from pyroaring import BitMap64
from bidict import bidict
import numpy as np
import pandas as pd
from graphein.protein.config import ProteinGraphConfig
from pkg.Insert import Insert
from pkg.Extract import Extract

class PDBGraphStore:
    def __init__(self, body_parts=None):
        self.__body_parts: dict | None = None
        self.granularity: str | None = None
        self.config: ProteinGraphConfig | None = None

        self.__set_body_parts(body_parts)
        self.inserter = Insert(self)
        self.extractor = Extract(self)

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_list())} pdbs'

    def close(self):
        self.set_config(None)
        self.set_granularity(None)
        self.__set_body_parts(None)
    
    def __set_body_parts(self, body_parts: dict=None):
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
                "edge_attr_values": bidict(),
                "node_attr_keyvalue_mapping": {},
                "edge_attr_keyvalue_mapping": {},
                "node_coords_values": bidict(), #array de float32
                "node_b_factor_values": bidict() #array de float64
            }

    def set_granularity(self, granularity: str | None):
        self.granularity = granularity
    
    __set_granularity = set_granularity

    def set_config(self, config: dict | ProteinGraphConfig | None):
        self.config = config
        if config is not None:
            if hasattr(config, 'dict'):
                granularity = config.dict()['granularity']
            elif isinstance(config, dict):
                granularity = config.get('granularity')
            else:
                granularity = getattr(config, 'granularity', None)
            self.set_granularity(granularity)
        else:
            self.set_granularity(None)

    __set_config = set_config

    def get_granularity(self):
        return self.granularity

    __get_granularity = get_granularity

    def get_config(self):
        return self.config
    
    def print_attr(self):
        for k, v in self.__body_parts.items():
            print(f"{k}: {v}")
            print("\n")
    
    def get_body_parts(self):
        return self.__body_parts

    def get_pdb_list(self):
        return self.__body_parts["pdb_code_to_id"].keys()
    
    @staticmethod
    def edge_label_undirected(edge_label: tuple) -> tuple:
        return tuple(sorted(edge_label))

    __edge_label_undirected = edge_label_undirected

    def get_meiler_by_residue(self, residue_name: str) -> pd.Series:
        '''
        retorna o Meiler de um residuo especifico
        '''
        return self.extractor.get_meiler_by_residue(residue_name)

    __get_meiler_by_residue = get_meiler_by_residue

    def get_one_hot_by_residue(self, residue_name: str) -> np.ndarray:
        '''
        retorna o one_hot array de um residuo especifico
        '''
        return self.extractor.get_one_hot_by_residue(residue_name)

    __get_one_hot_by_residue = get_one_hot_by_residue

    def insert(self, pdb_to_insert: dict):
        '''
        input: dict[str: nx.Graph]
        '''
        return self.inserter.insert(pdb_to_insert)

    def extract(self, pdb_to_extract: str) -> nx.Graph:
        return self.extractor.extract(pdb_to_extract)
