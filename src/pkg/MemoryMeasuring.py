from pkg.PDBGraphStore import PDBGraphStore
import pickle as pk
from pympler import asizeof
import numpy as np 

class MemoryMeasuring:
    def __init__(self, pdbObject: PDBGraphStore):
        self.__body_parts = pdbObject.get_body_parts()

    def pdb_code_to_id_memory(self):
        return asizeof.asizeof(self.__body_parts["pdb_code_to_id"])/1024/1024
    
    def pdb_code_to_id_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["pdb_code_to_id"]))/1024/1024

    def node_label_to_node_id_memory(self):
        return asizeof.asizeof(self.__body_parts["node_label_to_node_id"])/1024/1024
    
    def node_label_to_node_id_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["node_label_to_node_id"]))/1024/1024

    def edge_label_to_edge_id_memory(self):
        return asizeof.asizeof(self.__body_parts["edge_label_to_edge_id"])/1024/1024
    
    def edge_label_to_edge_id_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["edge_label_to_edge_id"]))/1024/1024

    def pdb_id_to_nodes_memory(self):
        return asizeof.asizeof(self.__body_parts["pdb_id_to_nodes"])/1024/1024
    
    def pdb_id_to_nodes_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["pdb_id_to_nodes"]))/1024/1024
    
    def pdb_id_to_edges_memory(self):
        return asizeof.asizeof(self.__body_parts["pdb_id_to_edges"])/1024/1024
    
    def pdb_id_to_edges_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["pdb_id_to_edges"]))/1024/1024

    def node_attr_keys_memory(self):
        keys = self.__body_parts.get("node_attr_keys", self.__body_parts.get("node_local_attr_keys"))
        if keys is not None:
            return asizeof.asizeof(keys)/1024/1024
        return 0.0

    def node_local_attr_keys_memory(self):
        return self.node_attr_keys_memory()

    def edge_attr_keys_memory(self):
        return asizeof.asizeof(self.__body_parts["edge_attr_keys"])/1024/1024

    def attr_keys_memory(self):
        return (
            self.node_attr_keys_memory() +
            self.edge_attr_keys_memory()
        )
    
    def attr_keys_serialized_memory(self):
        node_keys = self.__body_parts.get("node_attr_keys", self.__body_parts.get("node_local_attr_keys"))
        node_size = len(pk.dumps(node_keys))/1024/1024 if node_keys is not None else 0.0
        edge_keys = self.__body_parts.get("edge_attr_keys")
        edge_size = len(pk.dumps(edge_keys))/1024/1024 if edge_keys is not None else 0.0
        return node_size + edge_size

    def edge_attr_values_memory(self):
        return asizeof.asizeof(self.__body_parts["edge_attr_values"])/1024/1024
    
    def edge_attr_values_serialized_memory(self):
        return len(pk.dumps(self.__body_parts["edge_attr_values"]))/1024/1024

    def node_coords_values_memory(self):
        if "node_coords_values" in self.__body_parts:
            return asizeof.asizeof(self.__body_parts["node_coords_values"])/1024/1024
        return 0.0

    def node_coords_values_serialized_memory(self):
        if "node_coords_values" in self.__body_parts:
            return len(pk.dumps(self.__body_parts["node_coords_values"]))/1024/1024
        return 0.0

    def node_b_factor_values_memory(self):
        if "node_b_factor_values" in self.__body_parts:
            return asizeof.asizeof(self.__body_parts["node_b_factor_values"])/1024/1024
        return 0.0

    def node_b_factor_values_serialized_memory(self):
        if "node_b_factor_values" in self.__body_parts:
            return len(pk.dumps(self.__body_parts["node_b_factor_values"]))/1024/1024
        return 0.0
    
    def node_attr_values_memory(self):
        if "node_attr_values" in self.__body_parts:
            return asizeof.asizeof(self.__body_parts["node_attr_values"])/1024/1024
        return self.node_coords_values_memory() + self.node_b_factor_values_memory()
    
    def node_attr_values_serialized_memory(self):
        if "node_attr_values" in self.__body_parts:
            return len(pk.dumps(self.__body_parts["node_attr_values"]))/1024/1024
        return self.node_coords_values_serialized_memory() + self.node_b_factor_values_serialized_memory()

    def node_attr_keyvalue_mapping_memory(self):
        vetor = self.__body_parts.get("node_attr_keyvalue_mapping", self.__body_parts.get("node_local_attr_keyvalue_mapping"))
        return asizeof.asizeof(vetor)/1024/1024
    
    def node_attr_keyvalue_mapping_serialized_memory(self):
        vetor = self.__body_parts.get("node_attr_keyvalue_mapping", self.__body_parts.get("node_local_attr_keyvalue_mapping"))
        return len(pk.dumps(vetor))/1024/1024

    def node_local_attr_keyvalue_mapping_memory(self):
        return self.node_attr_keyvalue_mapping_memory()
    
    def node_local_attr_keyvalue_mapping_serialized_memory(self):
        return self.node_attr_keyvalue_mapping_serialized_memory()

    def edge_attr_keyvalue_mapping_memory(self):
        vetor = self.__body_parts.get("edge_attr_keyvalue_mapping", self.__body_parts.get("edge_local_attr_keyvalue_mapping"))
        return asizeof.asizeof(vetor)/1024/1024
    
    def edge_attr_keyvalue_mapping_serialized_memory(self):
        vetor = self.__body_parts.get("edge_attr_keyvalue_mapping", self.__body_parts.get("edge_local_attr_keyvalue_mapping"))
        return len(pk.dumps(vetor))/1024/1024

    def edge_local_attr_keyvalue_mapping_memory(self):
        return self.edge_attr_keyvalue_mapping_memory()
    
    def edge_local_attr_keyvalue_mapping_serialized_memory(self):
        return self.edge_local_attr_keyvalue_mapping_serialized_memory()

    def graph_structure_memory(self):
        return (
        self.pdb_code_to_id_memory() +
        self.pdb_id_to_nodes_memory() +
        self.pdb_id_to_edges_memory() +
        self.node_label_to_node_id_memory() +
        self.edge_label_to_edge_id_memory() 
        )

    def dict_attributes_memory(self):
        return (
        self.attr_keys_memory() +
        self.edge_attr_values_memory() +
        self.node_attr_values_memory()
        )

    def node_attributes_memory(self):
        return (
        self.node_attr_keyvalue_mapping_memory()
        )

    def edge_attributes_memory(self):
        return self.edge_attr_keyvalue_mapping_memory()
    
    def total_body_parts_memory(self):
        return asizeof.asizeof(self.__body_parts)/1024/1024

    def total_memory(self):
        return self.total_body_parts_memory()

    def total_body_parts_serialized_memory(self):
        return len(pk.dumps(self.__body_parts))/1024/1024

    def total_serialized_memory(self):
        return self.total_body_parts_serialized_memory()