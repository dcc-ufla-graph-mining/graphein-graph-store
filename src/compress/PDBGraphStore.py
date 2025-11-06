import networkx as nx
import numpy as np
from bidict import bidict
from pympler import asizeof
from pyroaring import BitMap64
from sortedcontainers import SortedSet
import pandas as pd
import edge_functions_Model as edgeModel
import struct
from Indexed_set import IndexedSet

class PDBGraphStore:
    def __init__(self, body_parts):
        if body_parts:
            self.__body_parts = body_parts
            print(self.__body_parts.keys())
        else:
            raise ValueError("'body_parts' should not be None")

    def __str__(self):
        return f'PDBGraphStore with {len(self.get_pdb_code_list())} pdbs'

    def get_pdb_code_list(self):
        return list(self.__body_parts["all_pdb_codes"])

    def __reconstruct_node_attributes(self, extracted_graph, pdb_code):
        pdb_idx = self.__body_parts["all_pdb_codes"].index(pdb_code)
        
        for node_attr_idx in self.__body_parts["pdb_to_nodes"][pdb_idx]:
            node_pair = self.__body_parts["node_attrs_unique"][node_attr_idx]
            node_id = node_pair[0]
            node = self.__body_parts["node_to_id"].inverse[node_id]
            
            chain_id_idx = self.__body_parts["node_attrs_global"][node][0]
            residue_name_idx = self.__body_parts["node_attrs_global"][node][1]
            atom_type_idx = self.__body_parts["node_attrs_global"][node][2]
            element_symbol_idx = self.__body_parts["node_attrs_global"][node][3]
            meiler_idx = self.__body_parts["node_attrs_global"][node][4]

            residue_number_idx = node_pair[1]
            coords_idx = node_pair[2]
            b_factor_idx = node_pair[3]
            
            extracted_graph.nodes[node]["chain_id"] = self.__body_parts["node_attr_keys"]["chain_id"].get(chain_id_idx)
            extracted_graph.nodes[node]["residue_name"] = self.__body_parts["node_attr_keys"]["residue_name"].get(residue_name_idx)
            extracted_graph.nodes[node]["residue_number"] = self.__body_parts["node_attr_keys"]["residue_number"].get(residue_number_idx)
            extracted_graph.nodes[node]["atom_type"] = self.__body_parts["node_attr_keys"]["atom_type"].get(atom_type_idx)
            extracted_graph.nodes[node]["element_symbol"] = self.__body_parts["node_attr_keys"]["element_symbol"].get(element_symbol_idx)
            extracted_graph.nodes[node]["coords"] = np.array(self.__body_parts["node_attr_keys"]["coords"].get(coords_idx))
            extracted_graph.nodes[node]["b_factor"] = self.__body_parts["node_attr_keys"]["b_factor"].get(b_factor_idx)
            
            meiler_value = self.__body_parts["node_attr_keys"]["meiler"].get(meiler_idx)
            extracted_graph.nodes[node]["meiler"] = pd.Series(
                meiler_value[0], 
                name=''.join(list(meiler_value[1])), 
                index=[f'dim_{x}' for x in [1, 2, 3, 4, 5, 6, 7]]
            )

    def __reconstruct_edge_attributes(self, extracted_graph, edges, pdb_idx, edge_funcs):
        for edge_dist_idx in self.__body_parts["pdb_to_edges"][pdb_idx]:
            edge_pair = self.__body_parts["edge_attrs"][edge_dist_idx]
            
            edge_id = edge_pair[0]
            distance = edge_pair[1]
            kind_indexes = edge_pair[2]
            
            edge = self.__body_parts["edge_to_id"].inverse[edge_id]
            
            if edge not in edges:
                continue

            kinds = [self.__body_parts["edge_kind_keys"].get(k) for k in kind_indexes]
            
            kinds = [k for k in kinds if k in edge_funcs]
            
            if not kinds:
                continue
            
            if not extracted_graph.has_edge(*edge):
                extracted_graph.add_edge(*edge)

            extracted_graph.edges[edge].setdefault("kind", set())
            for k in kinds:
                extracted_graph.edges[edge]["kind"].add(k)
                
            extracted_graph.edges[edge]["distance"] = distance

    def __extract_nodes(self, pdb_idx):
        def get_node_id(node_attr_unique_idx):
            node_pair = self.__body_parts["node_attrs_unique"][node_attr_unique_idx]
            return node_pair[0]
        
        return [self.__body_parts["node_to_id"].inverse[get_node_id(node_id)] 
                for node_id in self.__body_parts["pdb_to_nodes"][pdb_idx]]

    def __extract_edges(self, pdb_idx):
        def get_edge_id(edge_attrs_idx):
            edge_pair = self.__body_parts["edge_attrs"][edge_attrs_idx]
            return edge_pair[0]
        
        return [self.__body_parts["edge_to_id"].inverse[get_edge_id(edge_attrs_idx)] 
                for edge_attrs_idx in self.__body_parts["pdb_to_edges"][pdb_idx]]
    
    def extract_pdb_graphs(self, pdb_codes=[], edge_construction_functions=[]):
        if not pdb_codes:
            raise ValueError("pdb_codes should not be empty")

        if not edge_construction_functions:
            raise ValueError("edge_construction_functions should not be empty")

        extracted_graphs = []

        for pdb_code in pdb_codes:
            pdb_idx = self.__body_parts["all_pdb_codes"].index(pdb_code)
            
            extracted_graph = nx.Graph()
            extracted_graph.graph["pdb_code"] = pdb_code

            # Extrai nós e arestas
            extracted_nodes = self.__extract_nodes(pdb_idx)
            extracted_edges = self.__extract_edges(pdb_idx)
            
            # Adiciona os nós ao grafo
            extracted_graph.update(nodes=extracted_nodes)
            
            # Reconstrói os atributos dos nós e arestas
            self.__reconstruct_node_attributes(extracted_graph, pdb_code)
            self.__reconstruct_edge_attributes(extracted_graph, extracted_edges, pdb_idx, edge_construction_functions)
            
            extracted_graphs.append(extracted_graph)
        
        return extracted_graphs

    def __insert_node_attrs(self, node_to_insert):
        attr_indexes = []

        for value in self.__body_parts["node_attr_keys"]:
            attr_value = node_to_insert[value]
                        
            if isinstance(attr_value, pd.Series):
                attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
            elif isinstance(attr_value, np.ndarray):
                attr_value = tuple([tuple(attr_value)])
                            
            self.__body_parts["node_attr_keys"][value].add(attr_value)
            idx = self.__body_parts["node_attr_keys"][value].index(attr_value)
            attr_indexes.append(idx)

        return attr_indexes

    def __insert_nodes(self, pdb_idx, graph):
        if pdb_idx not in self.__body_parts["pdb_to_nodes"]:
            self.__body_parts["pdb_to_nodes"][pdb_idx] = BitMap64()
        
        for node in graph.nodes():
            if node not in self.__body_parts["node_to_id"]:
                new_node_id = len(self.__body_parts["node_to_id"])
                self.__body_parts["node_to_id"][node] = new_node_id

                self.__body_parts["node_attrs"][node] = self.__insert_node_attrs(node)

            node_id = self.__body_parts["node_to_id"][node]
            self.__body_parts["pdb_to_nodes"][pdb_idx].add(node_id)

    def __insert_edge_kinds(self, edge_to_insert, kinds_to_insert):
        kind_indexes = set()
        attr_indexes = []

        for kind in kinds_to_insert:
            if kind not in list(edgeModel.edge_functions_dict.values()):
                raise ValueError(f"edge_kind should be in the kind list, based on edge_construction_funcs. Error at \
                                     {edge_to_insert}: {kinds_to_insert}. kind {kind} not in edge_construction funcs")

            kind_indexes.add(self.__body_parts["edge_kind_keys"].index(kind))

        attr_indexes.append(kind_indexes)
                
        self.__body_parts["edge_kinds"][edge_to_insert] = attr_indexes

    def __handle_with_existent_edge_insertion(self, edge_to_insert, kinds_to_insert):
        idx = set([self.__body_parts["edge_kind_keys"].index(k) for k in kinds_to_insert])
        if idx:
            for i in idx:
                self.__body_parts["edge_kinds"][edge_to_insert].add(i)

    def __insert_edge_distances(self, edge_to_insert, distance, pdb_idx):
        edge_id = self.__body_parts["edge_to_id"][edge_to_insert]
        edge_pair = tuple([edge_id, distance])

        self.__body_parts["edge_distances"].append(edge_pair)
        self.__body_parts["pdb_to_edges"][pdb_idx].add(len(self.__body_parts["edge_distances"])-1)

    def __insert_edges(self, pdb_idx, graph):
        if pdb_idx not in self.__body_parts["pdb_to_edges"]:
            self.__body_parts["pdb_to_edges"][pdb_idx] = BitMap64()

        for u, v, data in graph.edges.data():
            edge = (u, v)
            if edge not in self.__body_parts["edge_to_id"]:
                new_edge_id = len(self.__body_parts["edge_to_id"])
                self.__body_parts["edge_to_id"][edge] = new_edge_id

                self.__insert_edge_kinds(edge, list(data["kind"]))

            else:
                self.__handle_with_existent_edge_insertion(edge, set(data["kind"]))

            self.__insert_edge_distances(edge, list(data["distance"]), pdb_idx)

    def insert_pdbs(self, graphs={}):
        if not graphs:
            raise ValueError("'graph' should not be empty")

        for pdb_code, graph_list in graphs.items():
            if pdb_code not in self.__body_parts["all_pdb_codes"]:
                self.__body_parts["all_pdb_codes"].add(pdb_code)

            pdb_idx = self.__body_parts["all_pdb_codes"].index(pdb_code)

            for graph in graph_list:
                self.__insert_nodes(pdb_idx, graph)
                self.__insert_edges(pdb_idx, graph)
        
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

