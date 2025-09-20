import networkx as nx
import numpy as np
from bidict import bidict
from pympler import asizeof
from pyroaring import BitMap64
from ordered_set import OrderedSet
import pandas as pd
import edge_functions_Model as edgeModel

class PDBGraphStore:
    def __init__(
            self, 
            node_to_id={},
            edge_to_id={},
            pdb_to_nodes={},
            pdb_to_edges={},
            node_attrs={},
            edge_attrs={},
            node_attr_keys={},
            edge_attr_keys={}
            ):
        self.node_to_id = node_to_id #mapeamento de de node para id, e vice versa, global
        self.edge_to_id = edge_to_id #mapeamento de de edge para id, e vice versa, global
        self.pdb_to_nodes = pdb_to_nodes #bitmap indicando quais pdbs cada node pertence
        self.pdb_to_edges = pdb_to_edges #bitmap indicando quais pdbs cada edge pertence
        self.node_attrs = node_attrs #lista de indices para o attr de cada node
        self.edge_attrs = edge_attrs #lista de indices para o attr de cada edge
        self.node_attr_keys = node_attr_keys #dicionario de atributos indexados para cada node
        self.edge_attr_keys = edge_attr_keys #dicionario de atributos indexados para cada edge

    def _reconstruct_node_attributes(self, extracted_graph, nodes):
        for node in nodes:
            if node in self.node_attrs:
                for i, key in enumerate(self.node_attr_keys):
                    index = self.node_attrs[node][i]
                    value = self.node_attr_keys[key][index]

                    if isinstance(value, tuple) and len(value) == 1:
                        value = np.array(value[0])

                    if isinstance(value, tuple) and len(value) == 3:
                        value = pd.Series(value[0], name=value[1], index=value[2])

                    extracted_graph.nodes[node][key] = value

    def _reconstruct_edge_attributes(self, extracted_graph, edges, edge_funcs):
        for u, v in edges:
            edge = u, v

            kinds = self.edge_attrs[edge][0]
            if len(self.edge_attrs[edge]) > 1:
                distance = self.edge_attrs[edge][1]

            kind_names = [self.edge_attr_keys["kind"][k] for k in kinds]
            kind_names = list(filter(lambda x: edgeModel.edge_functions_dict[x] in edge_funcs, kind_names))

            if len(kind_names) > 0:
                if not extracted_graph.has_edge(*edge):
                    extracted_graph.add_edge(*edge)

                extracted_graph.edges[edge].setdefault("kind", set())

                for kind in kind_names:
                    extracted_graph.edges[edge]["kind"].add(kind)

            try:
                extracted_graph.edges[edge]["distance"] = self.edge_attr_keys["distance"][self.edge_attrs[edge][1]]
            except Exception as e:
                print(e)

    def extract_pdb_graphs(self, pdb_codes=[], edge_construction_functions=[]):
        for pdb_code in pdb_codes:
            extracted_graph = nx.Graph()

            extracted_nodes = [self.node_to_id.inverse[self.node_id] for node_id in self.pdb_to_nodes[pdb_code][0]]
            extracted_edges = [self.edge_to_id.inverse[self.edge_id] for edge_id in self.pdb_to_edges[pdb_code][0]]

            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=extracted_nodes)

            self._reconstruct_node_attributes(extracted_graph, extracted_nodes)
            self._reconstruct_edge_attributes(extracted_graph, extracted_edges, edge_construction_functions)

            return extracted_graph

    def insert_pdb(self, pdb_code, graph):
        pass

    def remove_pdb(self, pdb_code):
        pass

    def node_to_id_size(self):
        return asizeof.asizeof(self.node_to_id) / 1024 / 1024
    
    def edge_to_id_size(self):
        return asizeof.asizeof(self.edge_to_id) / 1024 / 1024
    
    def pdb_to_nodes_size(self):
        return asizeof.asizeof(self.pdb_to_nodes) / 1024 / 1024
    
    def pdb_to_edges_size(self):
        return asizeof.asizeof(self.pdb_to_edges) / 1024 / 1024
    
    def node_attrs_size(self):
        return asizeof.asizeof(self.node_attrs) / 1024 / 1024
    
    def edge_attrs_size(self):
        return asizeof.asizeof(self.edge_attrs) / 1024 / 1024
    
    def node_attr_keys_size(self):
        return asizeof.asizeof(self.node_attr_keys) / 1024 / 1024
    
    def edge_attr_keys_size(self):
        return asizeof.asizeof(self.edge_attr_keys) / 1024 / 1024
    
    def calculate_graph_complete_space_size(self):
        return self.node_to_id_size() + \
            self.edge_to_id_size() + \
            self.pdb_to_nodes_size() + \
            self.pdb_to_edges_size() + \
            self.node_attrs_size() + \
            self.edge_attrs_size() + \
            self.node_attr_keys_size() + \
            self.edge_attr_keys_size()
    
    def calculate_total_nodes_size(self):
        return self.node_to_id_size() + \
            self.pdb_to_nodes_size() + \
            self.node_attrs_size() + \
            self.node_attr_keys_size()
    
    def calculate_total_edges_size(self):
        return self.edge_to_id_size() + \
            self.pdb_to_edges_size() + \
            self.edge_attrs_size() + \
            self.edge_attr_keys_size()
