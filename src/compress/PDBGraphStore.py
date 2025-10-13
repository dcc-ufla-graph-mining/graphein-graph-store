import networkx as nx
import numpy as np
from bidict import bidict
from pympler import asizeof
from pyroaring import BitMap64
from sortedcontainers import SortedSet
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
            edge_kinds={},
            node_attr_keys={},
            edge_kind_keys={},
            edge_distances={},
            all_pdb_codes={},
            ):
        self.node_to_id = node_to_id 
        self.edge_to_id = edge_to_id 
        self.pdb_to_nodes = pdb_to_nodes 
        self.pdb_to_edges = pdb_to_edges 
        self.node_attrs = node_attrs
        self.edge_kinds = edge_kinds 
        self.node_attr_keys = node_attr_keys 
        self.edge_kind_keys = edge_kind_keys 
        self.edge_distances = edge_distances
        self.all_pdb_codes = all_pdb_codes

        if len(node_to_id) == 0:
            self.init_attribute_keys()

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_list())} pdbs'

    def init_attribute_keys(self):
        node_attr_keys_list = ["chain_id", "residue_name", "residue_number", "atom_type", "element_symbol", "coords", "b_factor", "meiler"]

        for k in node_attr_keys_list:
            self.node_attr_keys[k] = SortedSet()
        
        self.edge_kind_keys = SortedSet(set([k for k, _ in edgeModel.edge_functions_dict.items()]))

    def get_len_edges(self):
        return len(self.edge_to_id.keys())
    
    def get_len_nodes(self):
        return len(self.node_to_id.keys())

    def get_pdb_list(self):
        return [pdb_code for pdb_code in self.all_pdb_codes]

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

    def _reconstruct_edge_attributes(self, extracted_graph, edges, edge_kinds, edge_kind_keys, 
                                 edge_funcs, pdb_to_edges, pdb_code, edge_distances, 
                                 edge_to_id, all_pdb_codes):
        
        pdb_idx = all_pdb_codes.index(pdb_code)
        
        for u, v in edges:
            edge = (u, v)
            
            if edge not in edge_kinds:
                continue
            
            kind_indexes_list = edge_kinds[edge]
            
            if len(kind_indexes_list) > 0:
                kind_indexes = kind_indexes_list[0]
                kind_names = [edge_kind_keys[k] for k in kind_indexes]
                kind_names = list(filter(lambda x: edgeModel.edge_functions_dict[x] in edge_funcs, kind_names))
                
                if len(kind_names) > 0:
                    if not extracted_graph.has_edge(*edge):
                        extracted_graph.add_edge(*edge)
                    
                    extracted_graph.edges[edge].setdefault("kind", set())
                    
                    for kind in kind_names:
                        extracted_graph.edges[edge]["kind"].add(kind)
        
        for edge_dist_idx in pdb_to_edges[pdb_idx]:
            edge_pair = edge_distances[edge_dist_idx]
            edge = edge_to_id.inverse[edge_pair[0]]
            extracted_graph.edges[edge]["distance"] = edge_pair[1]

    def extract_pdb_graphs(self, pdb_codes=[], edge_construction_functions=[]):
        extracted_graphs = []
        edge_funcs = [edgeModel.edge_functions_dict[x] for x in edge_construction_functions]
        print(f"dentro de extract_pdb. pdb_code={pdb_codes}, {edge_funcs}")
        
        for pdb_code in pdb_codes:
            pdb_idx = self.all_pdb_codes.index(pdb_code)
            extracted_graph = nx.Graph()
            
            # Extract nodes
            extracted_nodes = [self.node_to_id.inverse[node_id] for node_id in self.pdb_to_nodes[pdb_idx]]
            
            # Extract edges 
            def get_edge_id(edge_distances_idx):
                edge_pair = self.edge_distances[edge_distances_idx]
                edge_id = edge_pair[0]
                return edge_id
            
            extracted_edges = [self.edge_to_id.inverse[get_edge_id(edge_distances_idx)] 
                              for edge_distances_idx in self.pdb_to_edges[pdb_idx]]
            
            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=extracted_nodes)
            
            self._reconstruct_node_attributes(extracted_graph, extracted_nodes)
            
            self._reconstruct_edge_attributes(extracted_graph, extracted_edges, 
                                             self.edge_kinds, self.edge_kind_keys,
                                             edge_funcs,
                                             self.pdb_to_edges, pdb_code,
                                             self.edge_distances, self.edge_to_id,
                                             self.all_pdb_codes)
            
            extracted_graphs.append(extracted_graph)
        
        return extracted_graphs

    def insert_pdbs(self, graphs={}):
        if not graphs:
            return
        
        for pdb_code, graph_list in graphs.items():
            
            
            if pdb_code not in self.all_pdb_codes:
                self.all_pdb_codes.add(pdb_code)

            pdb_idx = self.all_pdb_codes.index(pdb_code)
            
            if pdb_idx not in self.pdb_to_nodes:
                self.pdb_to_nodes[pdb_idx] = BitMap64()
            if pdb_idx not in self.pdb_to_edges:
                self.pdb_to_edges[pdb_idx] = BitMap64()
            
            for graph in graph_list:
                for node in graph.nodes():
                    if node not in self.node_to_id:
                        new_node_id = len(self.node_to_id)
                        self.node_to_id[node] = new_node_id
                        
                        attr_indexes = []
                        for key in self.node_attr_keys:
                            attr_value = graph.nodes[node][key]
                            
                            if isinstance(attr_value, pd.Series):
                                attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
                            elif isinstance(attr_value, np.ndarray):
                                attr_value = tuple([tuple(attr_value)])
                            
                            attr_index = self.node_attr_keys[key].add(attr_value)
                            attr_indexes.append(attr_index)
                        
                        self.node_attrs[node] = attr_indexes

                    node_id = self.node_to_id[node]
                    self.pdb_to_nodes[pdb_idx].add(node_id)
                
                for u, v, data in graph.edges.data():
                    edge = (u, v)
                    
                    if edge not in self.edge_to_id:
                        new_edge_id = len(self.edge_to_id)
                        self.edge_to_id[edge] = new_edge_id
                        
                        attr_kind_value = list(data["kind"]) if "kind" in data else []
                        kind_indexes = set()
                        
                        for kind in attr_kind_value:
                            if kind in [k for k, _ in edgeModel.edge_functions_dict.items()]:
                                kind_index = self.edge_kind_keys.add(kind)
                                kind_indexes.add(kind_index)
                        
                        self.edge_kinds[edge] = [kind_indexes]
                        
                        if "distance" in data and data["distance"] is not None:
                            if new_edge_id not in self.edge_distances:
                                self.edge_distances[new_edge_id] = {}
                            self.edge_distances[new_edge_id][pdb_idx] = data["distance"]
                    
                    else:
                        if "kind" in data:
                            kinds = set(filter(lambda x: x, data["kind"]))
                            for kind in kinds:
                                if kind in [k for k, _ in edgeModel.edge_functions_dict.items()]:
                                    kind_index = self.edge_kind_keys.add(kind)
                                    self.edge_kinds[edge][0].add(kind_index)

                        if "distance" in data and data["distance"] is not None:
                            edge_id = self.edge_to_id[edge]
                            if edge_id not in self.edge_distances:
                                self.edge_distances[edge_id] = {}
                            self.edge_distances[edge_id][pdb_idx] = data["distance"]
                    
                    edge_id = self.edge_to_id[edge]
                    self.pdb_to_edges[pdb_idx].add(edge_id)
                    print(f'edge {edge} added')
        
        self.pdb_list = [k for k in self.pdb_to_nodes.keys()]
        
        print(f"Inseridos {len(graphs)} novos PDBs no supergrafo.")
        print(f"Total de nós únicos: {len(self.node_to_id)}")
        print(f"Total de arestas únicas: {len(self.edge_to_id)}")
        print(f"Total de PDBs: {len(self.pdb_list)}")

#daqui pra baixo nao mudou ainda
#TODO aplicar a mudanca das keys de pdb_to_edge e node que foi de str para int
    def remove_pdb(self, pdb_code):        
        if pdb_code not in self.pdb_to_nodes:
            print(f"PDB {pdb_code} não encontrado no supergrafo.")
            return {"removed": False, "reason": "PDB not found"}
        
        stats = {
            "removed": True,
            "nodes_removed": 0,
            "edges_removed": 0,
            "nodes_kept": 0,
            "edges_kept": 0,
            "pdb_code": pdb_code
        }
        
        nodes_to_check = list(self.pdb_to_nodes[pdb_code])
        edges_to_check = list(self.pdb_to_edges[pdb_code])
        
        nodes_to_remove = []
        edges_to_remove = []
        
        for node_id in nodes_to_check:
            node = self.node_to_id.inverse[node_id]
            is_exclusive = True
            
            for other_pdb, bitmap in self.pdb_to_nodes.items():
                if other_pdb != pdb_code and node_id in bitmap:
                    is_exclusive = False
                    break
            
            if is_exclusive:
                nodes_to_remove.append((node_id, node))
            else:
                stats["nodes_kept"] += 1
        
        pdb_code_index = self.all_pdb_codes.index(pdb_code)
        
        for edge_id in edges_to_check:
            edge = self.edge_to_id.inverse[edge_id]
            is_exclusive = True
            
            for other_pdb, bitmap in self.pdb_to_edges.items():
                if other_pdb != pdb_code and edge_id in bitmap:
                    is_exclusive = False
                    break
            
            if is_exclusive:
                edges_to_remove.append((edge_id, edge))
            else:
                stats["edges_kept"] += 1
                # Remove a distância específica deste pdb_code
                if edge_id in self.edge_distances:
                    if pdb_code_index in self.edge_distances[edge_id]:
                        del self.edge_distances[edge_id][pdb_code_index]
        
        for node_id, node in nodes_to_remove:
            del self.node_to_id[node]
            
            if node in self.node_attrs:
                del self.node_attrs[node]
            
            stats["nodes_removed"] += 1
        
        for edge_id, edge in edges_to_remove:
            del self.edge_to_id[edge]
            
            if edge in self.edge_kinds:
                del self.edge_kinds[edge]
            
            if edge_id in self.edge_distances:
                del self.edge_distances[edge_id]
            
            stats["edges_removed"] += 1
        
        del self.pdb_to_nodes[pdb_code]
        del self.pdb_to_edges[pdb_code]
        
        self.pdb_list = [k for k in self.pdb_to_nodes.keys()]
        
        if nodes_to_remove or edges_to_remove:
            self._recompact_ids(nodes_to_remove, edges_to_remove)
        
        print(f"PDB {pdb_code} removido com sucesso:")
        print(f"  - Nós removidos: {stats['nodes_removed']}")
        print(f"  - Arestas removidas: {stats['edges_removed']}")
        print(f"  - Nós mantidos (compartilhados): {stats['nodes_kept']}")
        print(f"  - Arestas mantidas (compartilhadas): {stats['edges_kept']}")
        print(f"  - Total de PDBs restantes: {len(self.pdb_list)}")
        
        return stats

    def _recompact_ids(self, removed_nodes, removed_edges):        
        if removed_nodes:
            removed_node_ids = sorted([node_id for node_id, _ in removed_nodes], reverse=True)
            
            for pdb_code, bitmap in self.pdb_to_nodes.items():
                for removed_id in removed_node_ids:
                    bitmap.discard(removed_id)
                    
                    new_bitmap = BitMap64()
                    for node_id in bitmap:
                        if node_id > removed_id:
                            new_bitmap.add(node_id - 1)
                        else:
                            new_bitmap.add(node_id)
                    self.pdb_to_nodes[pdb_code] = new_bitmap
            
            old_mapping = dict(self.node_to_id)
            self.node_to_id.clear()
            
            id_offset = 0
            for old_id in sorted(old_mapping.values()):
                while id_offset in [node_id for node_id, _ in removed_nodes]:
                    id_offset += 1
                
                node = None
                for n, oid in old_mapping.items():
                    if oid == old_id:
                        node = n
                        break
                
                if node is not None:
                    self.node_to_id[node] = old_id - sum(1 for rid, _ in removed_nodes if rid < old_id)
        
        if removed_edges:
            removed_edge_ids = sorted([edge_id for edge_id, _ in removed_edges], reverse=True)
            
            for pdb_code, bitmap in self.pdb_to_edges.items():
                for removed_id in removed_edge_ids:
                    bitmap.discard(removed_id)
                    
                    new_bitmap = BitMap64()
                    for edge_id in bitmap:
                        if edge_id > removed_id:
                            new_bitmap.add(edge_id - 1)
                        else:
                            new_bitmap.add(edge_id)
                    self.pdb_to_edges[pdb_code] = new_bitmap
            
            old_mapping = dict(self.edge_to_id)
            self.edge_to_id.clear()
            
            # Criar novo mapeamento de edge_distances
            new_edge_distances = {}
            
            for old_id in sorted(old_mapping.values()):
                edge = None
                for e, oid in old_mapping.items():
                    if oid == old_id:
                        edge = e
                        break
                
                if edge is not None and (old_id, edge) not in removed_edges:
                    new_id = old_id - sum(1 for rid, _ in removed_edges if rid < old_id)
                    self.edge_to_id[edge] = new_id
                    
                    # Atualiza edge_distances com o novo ID
                    if old_id in self.edge_distances:
                        new_edge_distances[new_id] = self.edge_distances[old_id]
            
            self.edge_distances = new_edge_distances

    def remove_multiple_pdbs(self, pdb_codes):
        if not pdb_codes:
            return {"removed": False, "reason": "No PDB codes provided"}
        
        combined_stats = {
            "removed": True,
            "nodes_removed": 0,
            "edges_removed": 0,
            "nodes_kept": 0,
            "edges_kept": 0,
            "pdbs_removed": [],
            "pdbs_not_found": []
        }
        
        valid_pdbs = []
        for pdb_code in pdb_codes:
            if pdb_code in self.pdb_to_nodes:
                valid_pdbs.append(pdb_code)
            else:
                combined_stats["pdbs_not_found"].append(pdb_code)
        
        if not valid_pdbs:
            combined_stats["removed"] = False
            combined_stats["reason"] = "No valid PDBs found"
            return combined_stats
        
        for pdb_code in valid_pdbs:
            stats = self.remove_pdb(pdb_code)
            if stats["removed"]:
                combined_stats["nodes_removed"] += stats["nodes_removed"]
                combined_stats["edges_removed"] += stats["edges_removed"]
                combined_stats["nodes_kept"] += stats["nodes_kept"]
                combined_stats["edges_kept"] += stats["edges_kept"]
                combined_stats["pdbs_removed"].append(pdb_code)
        
        return combined_stats

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
        return (asizeof.asizeof(self.edge_distances) / 1024 / 1024) + (asizeof.asizeof(self.edge_kinds) /1024/1024)
    
    def node_attr_keys_size(self):
        return asizeof.asizeof(self.node_attr_keys) / 1024 / 1024
    
    def edge_kind_keys_size(self):
        return asizeof.asizeof(self.edge_kind_keys) / 1024 / 1024
    
    def calculate_graph_complete_space_size(self):
        return self.node_to_id_size() + \
            self.edge_to_id_size() + \
            self.pdb_to_nodes_size() + \
            self.pdb_to_edges_size() + \
            self.node_attrs_size() + \
            self.edge_attrs_size() + \
            self.node_attr_keys_size() + \
            self.edge_kind_keys_size()
    
    def calculate_total_nodes_size(self):
        return self.node_to_id_size() + \
            self.pdb_to_nodes_size() + \
            self.node_attrs_size() + \
            self.node_attr_keys_size()
    
    def calculate_total_edges_size(self):
        return self.edge_to_id_size() + \
            self.pdb_to_edges_size() + \
            self.edge_attrs_size() + \
            self.edge_kind_keys_size()
