import networkx as nx
from pyroaring import BitMap64
from bidict import bidict
import pandas as pd
import pickle as pk
from pympler import asizeof

class PDBGraphStore:
    def __init__(self, body_parts):
        if body_parts:
            self.__body_parts = body_parts
            print(self.__body_parts.keys())
        else:
            raise ValueError("'body_parts' should not be None")

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_code_list())} pdbs'

    
