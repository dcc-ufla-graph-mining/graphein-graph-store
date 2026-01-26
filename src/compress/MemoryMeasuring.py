import PDBGraphStore
import pickle as pk
from pympler import asizeof 

class MemoryMeasuring:
    def __init__(self, pdbObject: PDBGraphStore):
        self.pdbObject = pdbObject

    def pdb_code_to_id_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["pdb_code_to_id"])/1024/1024
    
    def pdb_code_to_id_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["pdb_code_to_id"]))/1024/1024

    def node_label_to_node_id_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["node_label_to_node_id"])/1024/1024
    
    def node_label_to_node_id_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["node_label_to_node_id"]))/1024/1024

    def edge_label_to_edge_id_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["edge_label_to_edge_id"])/1024/1024
    
    def edge_label_to_edge_id_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["edge_label_to_edge_id"]))/1024/1024

    def pdb_id_to_nodes_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["pdb_id_to_nodes"])/1024/1024
    
    def pdb_id_to_nodes_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["pdb_id_to_nodes"]))/1024/1024
    
    def pdb_id_to_edges_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["pdb_id_to_edges"])/1024/1024
    
    def pdb_id_to_edges_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["pdb_id_to_edges"]))/1024/1024

    def node_global_attr_keys_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["node_global_attr_keys"])/1024/1024

    def node_local_attr_keys_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["node_local_attr_keys"])/1024/1024

    def edge_attr_keys_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["edge_attr_keys"])/1024/1024

    def attr_keys_memory(self):
        return (
            self.node_global_attr_keys_memory() +
            self.node_local_attr_keys_memory() +
            self.edge_attr_keys_memory()
        )
    
    def attr_keys_serialized_memory(self):
        return (
            len(pk.dumps(self.pdbObject.__body_parts["node_global_attr_keys"]))/1024/1024 +
            len(pk.dumps(self.pdbObject.__body_parts["node_local_attr_keys"]))/1024/1024 +
            len(pk.dumps(self.pdbObject.__body_parts["edge_attr_keys"]))/1024/1024
            )

    def edge_attr_values_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["edge_attr_values"])/1024/1024
    
    def node_attr_values_memory(self):
        return asizeof.asizeof(self.pdbObject.__body_parts["node_attr_values"])/1024/1024
    
    def edge_attr_values_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["edge_attr_values"]))/1024/1024
    
    def node_attr_values_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["node_attr_values"]))/1024/1024

    def node_global_attr_keyvalue_mapping_memory(self):
        vetor = self.pdbObject.__body_parts["node_global_attr_keyvalue_mapping"]

        return asizeof.asizeof(vetor)/1024/1024
    
    def node_global_attr_keyvalue_mapping_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["node_global_attr_keyvalue_mapping"]))/1024/1024

    def node_local_attr_keyvalue_mapping_memory(self):
        vetor = self.pdbObject.__body_parts["node_local_attr_keyvalue_mapping"]

        return asizeof.asizeof(vetor)/1024/1024
    
    def node_local_attr_keyvalue_mapping_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["node_local_attr_keyvalue_mapping"]))/1024/1024

    def edge_local_attr_keyvalue_mapping_memory(self):
        vetor = self.pdbObject.__body_parts["edge_local_attr_keyvalue_mapping"]

        return asizeof.asizeof(vetor)/1024/1024
    
    def edge_local_attr_keyvalue_mapping_serialized_memory(self):
        return len(pk.dumps(self.pdbObject.__body_parts["edge_local_attr_keyvalue_mapping"]))/1024/1024

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
        self.node_global_attr_keyvalue_mapping_memory() +
        self.node_local_attr_keyvalue_mapping_memory()
        )

    def edge_attributes_memory(self):
        return self.edge_local_attr_keyvalue_mapping_memory()
    
    def total_memory(self):
        return (
            self.graph_structure_memory() +
            self.dict_attributes_memory() +
            self.node_attributes_memory() +
            self.edge_attributes_memory() 
        )

    def total_serialized_memory(self):
        return (
            self.attr_keys_serialized_memory() +
            self.edge_attr_values_serialized_memory() +
            self.node_attr_values_serialized_memory() +
            self.edge_label_to_edge_id_serialized_memory() +
            self.edge_local_attr_keyvalue_mapping_serialized_memory() +
            self.node_global_attr_keyvalue_mapping_serialized_memory() +
            self.node_local_attr_keyvalue_mapping_serialized_memory() +
            self.pdb_id_to_edges_serialized_memory() +
            self.pdb_code_to_id_serialized_memory() +
            self.pdb_id_to_nodes_serialized_memory()
        )