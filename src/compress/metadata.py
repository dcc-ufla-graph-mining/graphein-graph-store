import pickle
import random
import time
import pandas as pd
import os
import metadata

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.graphs import compute_edges
from graphein.protein.utils import download_pdb

from graphein.protein.edges.atomic import (
    add_atomic_edges, #kind = covalent; bond_length = ?
    add_bond_order, #kind = SINGLE, DOUBLE, TRIPLE
    add_ring_status #kind = RING
)

from graphein.protein.edges.distance import (
    add_aromatic_interactions, #kind = aromatic
    add_aromatic_sulphur_interactions, #kind = aromatic_sulphur
    add_backbone_carbonyl_carbonyl_interactions, #kind = bb_carbonyl_carbonyl
    add_cation_pi_interactions, #kind = cation_pi
    add_distance_to_edges, #distance = ?
    add_distance_window, #kind = f"distance_window_{min}_{max}"
    add_delaunay_triangulation, #kind = delaunay
    # add_distance_threshold, #kind distance_threshold  // nao funciona eu ainda nao investiguei porque
    add_disulfide_interactions, #kind = disulfide
    add_fully_connected_edges, #kind = fully_connected
    add_hydrogen_bond_interactions, #kind = hbond
    add_hydrophobic_interactions, #kind = hydrophobic
    add_ionic_interactions, #kind = ionic
    add_k_nn_edges, #kind = knn  //obs: o nome é escolha do usuario e pode ser diferente do padrao knn
    add_peptide_bonds, #kind = peptide_bond
    add_pi_stacking_interactions, #kind = pi_stacking,
    add_t_stacking, #kind = t_stacking
    add_salt_bridges, #kind = salt_bridge
    add_vdw_interactions, #kind = vdw // obs: o nome é escolha do usuario e pode ser diferente do padrao vdw
    add_vdw_clashes, #kind = vdw_clash 
)


all_edge_funcs = [
    add_atomic_edges, 
    add_bond_order, 
    add_ring_status, 
    add_aromatic_interactions, 
    add_aromatic_sulphur_interactions, 
    add_backbone_carbonyl_carbonyl_interactions, 
    add_cation_pi_interactions, 
    add_distance_to_edges, 
    #add_distance_window, 
    add_delaunay_triangulation, 
    add_disulfide_interactions, 
    add_fully_connected_edges, 
    add_hydrogen_bond_interactions, 
    add_hydrophobic_interactions, 
    add_ionic_interactions, 
    add_k_nn_edges, 
    add_peptide_bonds, 
    add_pi_stacking_interactions, 
    add_t_stacking
]

edge_mutable_func_attributes = {
    "covalent": add_atomic_edges,
    "SINGLE": add_bond_order,
    "DOUBLE": add_bond_order,
    "TRIPLE": add_bond_order,
    "RING": add_ring_status,
}

edge_imutable_func_attributes = {
    "aromatic": add_aromatic_interactions,
    "aromatic_sulphur": add_aromatic_sulphur_interactions,
    "bb_carbonyl_carbonyl": add_backbone_carbonyl_carbonyl_interactions,
    "cation_pi": add_cation_pi_interactions,
    "distance": add_distance_to_edges,
    "distance_window": add_distance_window,
    "delaunay": add_delaunay_triangulation,
    "disulfide": add_disulfide_interactions,
    "fully_connected": add_fully_connected_edges,
    "hbond": add_hydrogen_bond_interactions,
    "hydrophobic": add_hydrophobic_interactions,
    "ionic": add_ionic_interactions,
    "knn": add_k_nn_edges,
    "peptide_bond": add_peptide_bonds,
    "pi_stacking": add_pi_stacking_interactions,
    "t_stacking": add_t_stacking,
    "salt_bridge": add_salt_bridges, 
}

import networkx as nx
import numpy as np
from bidict import bidict
from pympler import asizeof
from pyroaring import BitMap, BitMap64
from ordered_set import OrderedSet

############################################################################################

kind_attr = OrderedSet(set(edge_imutable_func_attributes.keys()))

print(kind_attr)

dataset = os.environ.get("DATASET")
file_path = os.path.dirname(os.path.realpath(metadata.__file__))
print(file_path)

if os.environ.get("DATA_DIR") is not None:
    data = os.environ.get("DATA_DIR")
else:
    data = os.path.abspath(f"{file_path}/../../data/")

errors_path = os.path.abspath(f"{file_path}/../../errors/")
results_path = os.path.abspath(f"{file_path}/../../results/")
    
pdb_dir = os.path.abspath(f"{data}/pdb_files/")

if not os.path.exists(pdb_dir):
    os.makedirs(pdb_dir)
if not os.path.exists(errors_path):
    os.makedirs(errors_path)
if not os.path.exists(results_path):
    os.makedirs(results_path)

dataset_name = dataset.split(".")[0]


#contar o tempo so ate a parte de buildar o grafo grandao, excluindo o extract e o assert
def _initialize_data_structures():
    """Inicializa as estruturas necessárias"""
    edge_to_pdbs = {}   
    node_to_pdbs = {}  
    edge_attrs = {}  
    node_attrs = {}
    pdb_to_nodes = {}
    pdb_to_edges = {}
    pdb_codes_config = {}
    
    return edge_to_pdbs, node_to_pdbs, edge_attrs, node_attrs, pdb_to_nodes, pdb_to_edges, pdb_codes_config


def _extract_attribute_keys():
    """Extrai as chaves dos atributos de nodes e edges dos grafos"""        
    
    edge_attr_keys_list = ["kind", "distance"]
    node_attr_keys_list = ['chain_id', 'residue_name', 'residue_number', 'atom_type', 'element_symbol', 'coords', 'b_factor', 'meiler']

    edge_attr_keys = {}
    node_attr_keys = {}

    edge_attr_keys["kind"] = kind_attr
    edge_attr_keys["distance"] = OrderedSet()

    # for key in node_attr_keys_list: 
    #     node_attr_keys[key] = OrderedSet()

    node_attr_keys["chain_id"] = OrderedSet()

    return edge_attr_keys, node_attr_keys


def _process_node_attributes(node, graph, node_attr_keys):
    """Processa os atributos de um nó específico"""
    attr_indexes = []
    
    for value in node_attr_keys:
        attr_value = graph.nodes[node][value]

        if isinstance(attr_value, pd.Series):
            attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
        elif isinstance(attr_value, np.ndarray):
            attr_value = tuple([tuple(attr_value)])

        attr_indexes.append(node_attr_keys[value].add(attr_value))

    return attr_indexes


def _process_edge_attributes(edge_data, edge_attr_keys):
    """Processa os atributos de uma aresta específica"""
    attr_indexes = []

    attr_kind_value = list(edge_data["kind"])
    attr_distance_value = edge_data["distance"]

    if attr_kind_value is None:
        raise ValueError("Edge kind attribute should not be None")
    
    if attr_distance_value is None:
        raise ValueError("Edge distance attribute should not be None")
    
    kind_indexes = set()

    for kind in attr_kind_value:
        if kind not in edge_mutable_func_attributes:
            kind_indexes.add(edge_attr_keys["kind"].index(kind))
        
    distance_indexes = edge_attr_keys["distance"].add(attr_distance_value)

    attr_indexes.append(kind_indexes)
    attr_indexes.append(distance_indexes)

    return attr_indexes



def _check_if_edge_attr_matches(edge, data, edge_attrs, edge_attr_keys):
    """Verifica e atualiza os atributos de uma aresta específica"""
    
    real_data = set(filter(lambda x: x not in edge_mutable_func_attributes, data["kind"]))

    diff = set([edge_attr_keys["kind"].index(k) for k in real_data]) - set(edge_attrs[edge][0])
    if len(diff) != 0:
        for i in diff:
            edge_attrs[edge][0].add(i)


def _process_nodes(protein_graphs, node_to_pdbs, node_attrs, node_attr_keys, pdb_to_nodes):
    """Processa todos os nós dos grafos"""
    for pdb_code, graphs in protein_graphs.items():
        for g in graphs:
            pdb_to_nodes[pdb_code] = [BitMap64()]
            
            for node in g.nodes():
                if node not in node_to_pdbs:
                    node_to_pdbs[node] = []
                    attr_indexes = _process_node_attributes(node, g, node_attr_keys)
                    node_attrs[node] = attr_indexes
                    
                if pdb_code not in node_to_pdbs[node]:
                    node_to_pdbs[node].append(pdb_code)


def _process_edges(protein_graphs, edge_to_pdbs, edge_attrs, edge_attr_keys, pdb_to_edges):
    """Processa todas as arestas dos grafos"""
    for pdb_code, graphs in protein_graphs.items():
        for g in graphs:
            if pdb_code not in pdb_to_edges:
                pdb_to_edges[pdb_code] = [BitMap64()]
            
            for u, v, data in g.edges(data=True):
                edge = (u, v)
                
                if edge not in edge_to_pdbs:
                    edge_to_pdbs[edge] = []
                    attr_indexes = _process_edge_attributes(data, edge_attr_keys)
                    edge_attrs[edge] = attr_indexes

                else:
                    _check_if_edge_attr_matches(edge, data, edge_attrs, edge_attr_keys)
                
                if pdb_code not in edge_to_pdbs[edge]:
                    edge_to_pdbs[edge].append(pdb_code)


def _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes):
    """Cria mapeamentos entre arestas/nós e seus IDs únicos"""

    edge_id = 0
    edge_to_id = {}
    for e in edge_to_pdbs:
        edge_to_id[e] = edge_id
        pdbs = edge_to_pdbs[e]
        for pdb_code in pdbs:
            pdb_to_edges[pdb_code][0].add(edge_id)
        edge_id += 1

    node_id = 0
    node_to_id = {}  
    for u in node_to_pdbs:
        node_to_id[u] = node_id
        pdbs = node_to_pdbs[u]
        for pdb_code in pdbs:
            pdb_to_nodes[pdb_code][0].add(node_id)
        node_id += 1

    node_to_id = bidict(node_to_id)  
    edge_to_id = bidict(edge_to_id)
    
    return node_to_id, edge_to_id


def _reconstruct_node_attributes(extracted_graph, nodes, node_attrs, node_attr_keys):
    """Reconstrói os atributos dos nós no grafo extraído."""
    for node in nodes:
        if node in node_attrs:
            for i, key in enumerate(node_attr_keys):
                index = node_attrs[node][i]
                value = node_attr_keys[key][index]

                if isinstance(value, tuple) and len(value) == 1:
                    value = np.array(value[0])
                elif isinstance(value, tuple) and len(value) == 3:
                    value = pd.Series(value[0], name=value[1], index=value[2])
                    
                extracted_graph.nodes[node][key] = value


def _reconstruct_edge_attributes(extracted_graph, edges, edge_attrs, edge_attr_keys, edge_funcs):
    """Reconstrói os atributos das arestas no grafo extraído."""
    for u, v in edges:
        edge = u, v
        kinds = edge_attrs[edge][0]
        distance = edge_attrs[edge][1]
        kind_names = [edge_attr_keys["kind"][k] for k in kinds]
        kind_names = list(filter(lambda x: edge_imutable_func_attributes[x] in edge_funcs, kind_names))

        if len(kind_names) > 0:
            try:
                extracted_graph.edges[edge]
            except KeyError:
                extracted_graph.add_edge(*edge)
            
            try:
                extracted_graph.edges[edge]["kind"]
            except KeyError:
                extracted_graph.edges[edge]["kind"] = set()

            for kind in kind_names:
                extracted_graph.edges[edge]["kind"].add(kind)


def _validate_graph_reconstruction(extracted_graph, original_graph, pdb_code):
    """Valida se o grafo foi reconstruído corretamente."""
    try:
        assert nx.utils.edges_equal(extracted_graph._adj, original_graph._adj)
    except Exception as e:
        print(e)
        print("error in edges")
        for (u, v, attr) in extracted_graph.edges.data():
            print(u, v, attr)
            for k, t in attr.items():
                print(type(t))
            break
        for (u, v, attr) in original_graph.edges.data():  
            print(u, v, attr) 
            for k, t in attr.items():
                print(type(t))
            break
        return False

    try:    
        assert nx.utils.nodes_equal(extracted_graph._node, original_graph._node)
    except Exception as e:
        print("error in nodes")
        with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
            f.write(f"Error in graph extraction for {pdb_code}: {e}\n")
        for node, attr in extracted_graph.nodes(data=True):
            print(node, attr)
            for k, t in attr.items():
                print(type(t))
            break

        for node, attr in original_graph.nodes(data=True):
            print(node, attr)
            for k, t in attr.items():
                print(type(t))
            break
        return False
    
    return True


def _reconstruct_and_validate_graphs(protein_graphs, node_to_id, edge_to_id, 
                                   pdb_to_nodes, pdb_to_edges, 
                                   node_attrs, edge_attrs, 
                                   node_attr_keys, edge_attr_keys, funcs):
    """Reconstrói e valida todos os grafos"""

    is_different = False

    for pdb_code in protein_graphs:
        pdb_graphs = protein_graphs[pdb_code]
        for graph in pdb_graphs:
            original_graph = graph.copy()
            print(original_graph.graph["config"].model_dump())

            nodes = [node_to_id.inverse[node_id] for node_id in pdb_to_nodes[pdb_code][0]]  
            edges = [edge_to_id.inverse[edge_id] for edge_id in pdb_to_edges[pdb_code][0]]

                
            extracted_graph = nx.Graph()
            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.update(nodes=nodes)

            _reconstruct_node_attributes(extracted_graph, nodes, node_attrs, node_attr_keys)
            _reconstruct_edge_attributes(extracted_graph, edges, edge_attrs, edge_attr_keys, original_graph.graph["config"].edge_construction_functions)

            print("original graph edge funcs: ", original_graph.graph["config"].edge_construction_functions)
            
            print("number of nodes in original graph: ", len(original_graph.nodes()))
            
            print("number of nodes in extracted graph: ", len(extracted_graph.nodes()))

            # print("original graph", original_graph.edges(data=True))
            # print("extracted graph", extracted_graph.edges(data=True))

            if (len(extracted_graph.edges()) != len(original_graph.edges())):
                is_different = True
                print("number of edges in extracted graph: ", len(extracted_graph.edges()))
                print("number of edges in original graph: ", len(original_graph.edges()))

                with open(f"extracted.txt", "w") as f:
                    
                    f.write("\n\n")
                    f.write("sample\n")
                    sample = random.sample(list(original_graph.edges), 10)

                    for e in sample:
                        f.write(str(extracted_graph.edges[e]))
                        f.write("\n")
                        f.write(str(original_graph.edges[e]))
                        f.write("\n")


                    f.write("\n\n")
                    f.write("differences between original and extracted graph:\n")
                    original_to_extracted = set(original_graph.edges) - set(extracted_graph.edges)
                    f.write(str(original_to_extracted))
                    f.write("\n\n")
                    for e in original_to_extracted:
                        f.write(f"edge {e} in original graph: ")
                        f.write(str(original_graph.edges[e]))
                        f.write("\n")
                    f.write("\n\n")
                    f.write("differences between extracted and original graph:\n")

                    extracted_to_original = set(extracted_graph.edges) - set(original_graph.edges)
                    for e in extracted_to_original:
                        f.write(str(e))
                        f.write(str(extracted_graph.edges[e]))
                        f.write("\n")
                    f.write("\n\n")
                    for e in extracted_to_original:
                        f.write(f"edge {e} in extracted graph: ")
                        f.write(str(extracted_graph.edges[e]))
                        f.write("\n")
                        for g in pdb_graphs:
                            
                            try:
                                text = f"edge {e} in graph {str(g.graph['config'].model_dump())}: \n {g.edges[e]}"
                                f.write(text)
                            except KeyError:
                                continue
                            f.write("\n")

                    f.write(str(original_graph.graph["config"].model_dump()))

                    f.write("\n\n")
                    for e in extracted_to_original:
                        f.write(str(edge_attrs[e][0]))
                        f.write("\n")
                    
                    for u, v in extracted_to_original:
                        f.write(f"node {u} in extracted graph: ")
                        f.write(str(extracted_graph.nodes[u]))
                        f.write("\n")
                        f.write(f"node {u} in original graph: ")
                        f.write(str(original_graph.nodes[u]))
                        f.write("\n\n")
                        f.write("all edges from this node in extracted graph:\n")
                        for e in extracted_graph.edges(u, data=True):
                            f.write(str(e))
                            f.write("\n")
                        f.write("\n\n")
                        f.write("all edges from this node in original graph:\n")
                        for e in original_graph.edges(u, data=True):
                            f.write(str(e))
                            f.write("\n")

                        f.write("\n\n")

                        f.write(f"node {v} in extracted graph: ")
                        f.write(str(extracted_graph.nodes[v]))
                        f.write("\n\n")
                        f.write("all edges from this node in extracted graph:\n")
                        for e in extracted_graph.edges(v, data=True):
                            f.write(str(e))
                            f.write("\n")
                        f.write("\n\n")
                        f.write(f"node {v} in original graph: ")
                        f.write(str(original_graph.nodes[v]))
                        f.write("\n\n")
                        f.write("all edges from this node in original graph:\n")
                        for e in extracted_graph.edges(v, data=True):
                            f.write(str(e))
                            f.write("\n")
                        f.write("\n\n")




            # if not _validate_graph_reconstruction(extracted_graph, original_graph, pdb_code):
            #     break

def _process_pdb_codes_config(protein_graphs, pdb_codes_config):
    """Processa as configurações dos códigos PDB"""

    def add_edge_construction_function(fun):
        """Adiciona uma função de construção de arestas à configuração do PDB"""

        if fun not in pdb_codes_config[pdb_code]["edge_construction_functions"]:
            pdb_codes_config[pdb_code]["edge_construction_functions"].append(fun)

    for pdb_code, g in protein_graphs.items():
        if pdb_code not in pdb_codes_config:
            pdb_codes_config[pdb_code] = {}
        if "edge_construction_functions" not in pdb_codes_config[pdb_code]:
            pdb_codes_config[pdb_code]["edge_construction_functions"] = []
        for graph in g:
            for fun in graph.graph["config"].edge_construction_functions:
                if fun in all_edge_funcs:
                    add_edge_construction_function(fun.__name__)


def compress_with_composition(protein_graphs):
    """
    Comprime grafos de proteínas usando composição de bitmaps.
    
    Args:
        protein_graphs: Dicionário de grafos de proteínas indexados por código PDB
        
    Returns:
        PDBGraphStoreBitmap: Estrutura comprimida contendo os grafos
    """

    edge_to_pdbs, node_to_pdbs, edge_attrs, node_attrs, pdb_to_nodes, pdb_to_edges, pdb_codes_config = _initialize_data_structures()
    
    edge_attr_keys, node_attr_keys = _extract_attribute_keys()
    
    _process_nodes(protein_graphs, node_to_pdbs, node_attrs, node_attr_keys, pdb_to_nodes)
    _process_edges(protein_graphs, edge_to_pdbs, edge_attrs, edge_attr_keys, pdb_to_edges)
    
    # print(edge_attr_keys.items())
    # print(edge_attrs.items())

    node_to_id, edge_to_id = _create_id_mappings(edge_to_pdbs, node_to_pdbs, pdb_to_edges, pdb_to_nodes)

    # print(node_to_id.items())

    _process_pdb_codes_config(protein_graphs, pdb_codes_config)

    # print(pdb_codes_config["1F7Z"].items())

    for k, v in pdb_codes_config.items():
        print(k)
        print(v)
    
    del edge_to_pdbs
    del node_to_pdbs

    # print("node_to_id: ", node_to_id)
    # print("edge_to_id: ", edge_to_id)

    # print("pdb_to_nodes: ", pdb_to_nodes)
    # print("pdb_to_edges: ", pdb_to_edges)

    # print("node_attrs: ", node_attrs)
    # print("edge_attrs: ", edge_attrs)

    # print("node_attr_keys: ", node_attr_keys)
    # print("edge_attr_keys: ", edge_attr_keys)

    funcs = [add_atomic_edges, add_bond_order, add_ring_status]
    
    _reconstruct_and_validate_graphs(protein_graphs, node_to_id, edge_to_id, 
                                   pdb_to_nodes, pdb_to_edges, 
                                   node_attrs, edge_attrs, 
                                   node_attr_keys, edge_attr_keys, funcs)
    

    return PDBGraphStoreBitmap(node_to_id, edge_to_id, 
                            pdb_to_nodes, pdb_to_edges, 
                            node_attrs, edge_attrs, 
                            edge_attr_keys, node_attr_keys)



def split_PDB_store(pdb_store, pdb_codes: list):
    #TODO
    pass

def merge_PDB_stores(pdb_store1, pdb_store2):
    #TODO
    pass

class PDBGraphStoreBitmap:
    def __init__(
            self, node_to_id={}, 
            edge_to_id={},
            pdb_to_nodes={}, 
            pdb_to_edges={}, 
            node_attrs={}, 
            edge_attrs={}, 
            edge_attr_keys={}, 
            node_attr_keys={}
            ):
        
        self.node_to_id = node_to_id #mapeamento de de node para id, e vice versa, global
        self.edge_to_id = edge_to_id #mapeamento de de edge para id, e vice versa, global
        self.pdb_to_nodes = pdb_to_nodes #bitmap indicando quais pdbs cada node pertence
        self.pdb_to_edges = pdb_to_edges #bitmap indicando quais pdbs cada edge pertence
        self.node_attrs = node_attrs #lista de indices para o attr de cada node
        self.edge_attrs = edge_attrs #lista de indices para o attr de cada edge
        self.node_attr_keys = node_attr_keys #dicionario de atributos indexados para cada node
        self.edge_attr_keys = edge_attr_keys #dicionario de atributos indexados para cada edge

    
    def _reconstruct_edge_attributes(self, extracted_graph, edges, edge_funcs):
        for u, v in edges:

            if (u, v) in self.edge_attrs:
                kinds = self.edge_attrs[(u, v)][0]
                distance = self.edge_attrs[(u, v)][1]

                kind_names = [self.edge_attr_keys["kind"][k] for k in kinds]
                
                for kind_name in kind_names:
                    if kind_name in edge_imutable_func_attributes:
                        if edge_imutable_func_attributes[kind_name] in edge_funcs:
                            if not "kind" in extracted_graph.edges[(u, v)].keys():
                                extracted_graph.edges[(u, v)]["kind"] = set()

                            extracted_graph.edges[(u, v)]["kind"].add(kind_name)

                if "kind" not in extracted_graph.edges[(u, v)].keys():
                    # print(f"Removing edge ({u}, {v}) from extracted graph because it has no kind attribute")
                    extracted_graph.remove_edge(u, v)
                else:
                    extracted_graph.edges[(u, v)]["distance"] = self.edge_attr_keys["distance"][distance]
            else:
                # print(f"Removing edge ({u}, {v}) from extracted graph because it has no attributes")
                extracted_graph.remove_edge(u, v)

        return extracted_graph

    def extract_pdb_graph(self, pdb_code, edge_construction_functions=[add_atomic_edges, add_peptide_bonds]):
        nodes_view = self.pdb_to_nodes.get(pdb_code, None)
        edges_view = self.pdb_to_edges.get(pdb_code, None)
        
        def union_bitmaps(bms):
            union_bm = BitMap64()
            if bms is None: 
                return union_bm
            for bm in bms: 
                union_bm = union_bm | bm
            return union_bm
        
        try:
            nodes = [self.node_to_id.inverse[node_id] for node_id in union_bitmaps(nodes_view)]
            edges = [self.edge_to_id.inverse[edge_id] for edge_id in union_bitmaps(edges_view)]
            
            extracted_graph = nx.Graph()
            extracted_graph.update(edges=edges, nodes=nodes)

            for node in nodes:
                if node in self.node_attrs:
                    for i, key in enumerate(self.node_attr_keys):
                        index = self.node_attrs[node][i]
                        value = self.node_attr_keys[key][index]
                        
                        if isinstance(value, tuple) and len(value) == 1:
                            value = np.array(value[0])
                        elif isinstance(value, tuple) and len(value) == 3:
                            value = pd.Series(value[0], name=value[1], index=value[2])
                            
                        extracted_graph.nodes[node][key] = value

            extracted_graph = self._reconstruct_edge_attributes(extracted_graph, edges, edge_construction_functions)

            extracted_graph.graph["pdb_code"] = pdb_code
            extracted_graph.graph["config"] = ProteinGraphConfig(
                granularity="N",
                edge_construction_functions=edge_construction_functions
            )
            
        except Exception as e:
            # print(f"Error extracting graph for {pdb_code}: {e}")
            extracted_graph = nx.Graph()
            extracted_graph.update(edges=[], nodes=[])
            with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
                f.write(f"Error extracting graph for {pdb_code}: {e}\n")

        return extracted_graph

    def insert_pdb(self, pdb_code, graph):
        #TODO
        pass

    def remove_pdb(self, pdb_code):
        #TODO
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
        return self.node_to_id_size() + self.edge_to_id_size() + self.pdb_to_nodes_size() + self.pdb_to_edges_size() + self.node_attrs_size() + self.edge_attrs_size() + self.node_attr_keys_size() + self.edge_attr_keys_size()
    def calculate_total_nodes_size(self):
        return self.node_to_id_size() + self.pdb_to_nodes_size() + self.node_attrs_size() + self.node_attr_keys_size()
    def calculate_total_edges_size(self):
        return self.edge_to_id_size() + self.pdb_to_edges_size() + self.edge_attrs_size() + self.edge_attr_keys_size()


def main():
    with open(f"{errors_path}/{dataset_name}_errors.log", "w") as f:
        f.write(f"Errors log for {dataset_name} \n")

    with open(f"{results_path}/{dataset_name}_results.txt", "w") as f:
        f.write(f"Results log for {dataset_name}\n")

    print("main")

    pdb_codes = []
    
    params_to_change_list = []
    config_list = []


    with open(f'{data}/{dataset}', 'r') as f:
        for line in f:
            pdb_codes.append(line.strip().upper())

    for i in range(len(pdb_codes)):
        size = random.randint(1, len(all_edge_funcs))
        edge_construction_functions = random.sample(all_edge_funcs, size)
        # edge_construction_functions = [add_atomic_edges, add_t_stacking, add_bond_order, add_delaunay_triangulation, add_hydrogen_bond_interactions]

        print(edge_construction_functions)

        params_to_change_list.append({"granularity": "N", 
                                    "edge_construction_functions": edge_construction_functions})
        
        config_list.append(ProteinGraphConfig(**params_to_change_list[i]))

    # for i in range(len(pdb_codes)):
    #     match i:
    #         case 0:
    #             edge_construction_functions = [add_k_nn_edges, add_hydrogen_bond_interactions]
    #         case 1:
    #             edge_construction_functions = [add_t_stacking, add_ionic_interactions, add_atomic_edges,
    #                                             add_hydrogen_bond_interactions, add_pi_stacking_interactions,
    #                                             add_backbone_carbonyl_carbonyl_interactions, 
    #                                             add_hydrophobic_interactions, add_distance_to_edges]
    #         case 2:
    #             edge_construction_functions = [add_aromatic_interactions, add_k_nn_edges, add_t_stacking,
    #                                             add_aromatic_sulphur_interactions, add_cation_pi_interactions,
    #                                             add_peptide_bonds, add_bond_order, add_hydrogen_bond_interactions,
    #                                             add_ring_status, add_atomic_edges, add_fully_connected_edges, add_hydrophobic_interactions,
    #                                             add_pi_stacking_interactions, add_disulfide_interactions, add_backbone_carbonyl_carbonyl_interactions,
    #                                             add_ionic_interactions, add_delaunay_triangulation, add_distance_to_edges]
    #         case 3:
    #             edge_construction_functions = [add_aromatic_interactions, add_aromatic_sulphur_interactions]
    #         case 4:
    #             edge_construction_functions = [add_peptide_bonds, add_hydrogen_bond_interactions, add_ring_status, add_cation_pi_interactions, 
    #                                             add_aromatic_interactions, add_hydrophobic_interactions, add_aromatic_sulphur_interactions, 
    #                                             add_delaunay_triangulation, add_t_stacking, add_pi_stacking_interactions, add_bond_order, add_ionic_interactions]
    #         case 5:
    #             edge_construction_functions = [add_cation_pi_interactions, add_hydrogen_bond_interactions, add_ionic_interactions, add_ring_status,
    #                                             add_disulfide_interactions, add_atomic_edges, add_t_stacking, add_delaunay_triangulation,
    #                                             add_backbone_carbonyl_carbonyl_interactions, add_aromatic_interactions, add_distance_to_edges, add_k_nn_edges, 
    #                                             add_bond_order, add_fully_connected_edges, add_aromatic_sulphur_interactions]
        
    #     params_to_change_list.append({"granularity": "N", 
    #                                   "edge_construction_functions": edge_construction_functions})
        
    #     config_list.append(ProteinGraphConfig(**params_to_change_list[i]))

    protein_graphs_with_data = {}
    protein_graphs_without_data = {}

    time_begin = time.time()

    i = 0

    pdb_codes_copy = pdb_codes.copy()

    number_of_nodes = []
    number_of_edges = []

    config = iter(config_list)

    for pdb_code in pdb_codes_copy:
        i += 1
        
        if os.path.exists(f"{pdb_dir}/{pdb_code}.pdb"):
            print(f"Reading {pdb_code} from local directory")
            try:
                pdb_file = os.path.abspath(f"{pdb_dir}/{pdb_code}.pdb")
            except Exception as e:
                print(f"Error reading {pdb_code}: {e}")
                pdb_codes.remove(pdb_code)
                

                with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
                    f.write(f"Error reading {pdb_code}: {e}\n")
                continue
        else:
            print(f"Downloading {pdb_code} from PDB")
            try:
                pdb_file = download_pdb(pdb_code, f"{pdb_dir}/")
                if pdb_file is None:
                    print(f"Failed to download {pdb_code}")
                    pdb_codes.remove(pdb_code)
                    with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
                        f.write(f"Failed to download {pdb_code}\n")
                    continue
            except Exception as e:
                print(f"Error downloading {pdb_code}: {e}")
                pdb_codes.remove(pdb_code)

                with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
                    f.write(f"Error downloading {pdb_code}: {e}\n")
                continue

        graph = construct_graph(config=next(config), path=pdb_file)
        print(len(graph.edges))
        print(graph.graph["pdb_code"])

        aux = graph.graph.copy()
        aux.clear()
        aux["config"] = graph.graph["config"]
        aux["pdb_code"] = graph.graph["pdb_code"]
        graph.graph = aux
        print(graph)

        try:
            protein_graphs_with_data[pdb_code].append(graph.copy())  
        except KeyError:
            protein_graphs_with_data[pdb_code] = []
            protein_graphs_with_data[pdb_code].append(graph.copy())  

        for node in graph.nodes():
            graph.nodes[node].clear()
        for u, v in graph.edges():
            graph.edges[u, v].clear()

        try:
            protein_graphs_without_data[pdb_code].append(graph.copy())
        except KeyError:
            protein_graphs_without_data[pdb_code] = []
            protein_graphs_without_data[pdb_code].append(graph.copy())

        number_of_nodes.append(len(graph.nodes()))
        number_of_edges.append(len(graph.edges()))

        del graph

    with open(f"{results_path}/{dataset_name}_results.txt", "w") as f:
        f.write(f"Average number of nodes: {np.mean(number_of_nodes)}\n")
        f.write(f"Average number of edges: {np.mean(number_of_edges)}\n")
    
    del number_of_nodes
    del number_of_edges

    del pdb_codes_copy


    time_end = time.time()

    construct_time = time_end - time_begin
    print("Time to construct graphs:", construct_time)
    print("Number of graphs:", len(protein_graphs_with_data))

    v_size = 0
    e_size = 0

    v_serialized = 0
    e_serialized = 0

    for k, v in protein_graphs_with_data.items():
        for g in v:
            v_size += asizeof.asizeof(g._node) / 1024 / 1024
            e_size += asizeof.asizeof(g._adj) / 1024 / 1024

            v_serialized += asizeof.asizeof(pickle.dumps(g._node)) / 1024 / 1024
            e_serialized += asizeof.asizeof(pickle.dumps(g._adj)) / 1024 / 1024

    time_begin = time.time()

    global_graph_obj = compress_with_composition(protein_graphs_with_data)

    time_end = time.time()
    compress_time = time_end - time_begin
    print("Time to compress:", compress_time)

    random.shuffle(pdb_codes)

    times_to_extract = []

    pdb_codes = set(pdb_codes)  
    for pdb_code in pdb_codes:

        for graph in protein_graphs_with_data[pdb_code]:
            extracted_graph = global_graph_obj.extract_pdb_graph(pdb_code, graph.graph["config"].edge_construction_functions)
            print("\n\n")
            print("Number of edges in original graph: ", len(graph.edges()))
            print("Number of nodes in original graph: ", len(graph.nodes()))
            print("Number of edges in extracted graph: ", len(extracted_graph.edges()))
            print("Number of nodes in extracted graph: ", len(extracted_graph.nodes()))
            print("\n\n")
            print("original graph config:", graph.graph["config"].model_dump())
            print("extracted graph config:", extracted_graph.graph["config"].model_dump())
            print("\n\n")
            print(pdb_codes, protein_graphs_with_data[pdb_code])

            del extracted_graph


        time_begin = time.time()
        g = global_graph_obj.extract_pdb_graph(pdb_code, protein_graphs_with_data[pdb_code][0].graph["config"].edge_construction_functions)
        time_end = time.time()
        times_to_extract.append(time_end - time_begin)

        try:
            assert nx.utils.nodes_equal(g._node, protein_graphs_with_data[pdb_code][0]._node) 
            assert nx.utils.edges_equal(g._adj, protein_graphs_with_data[pdb_code][0]._adj)
        except AssertionError as e:
            print(f"Error in graph extraction for {pdb_code}: {e}")
            print("Extracted graph nodes:", list(g.nodes(data=True))[0:10])
            print("Original graph nodes:", list(protein_graphs_with_data[pdb_code][0].nodes(data=True))[0:10])
            print("Extracted graph edges:", list(g.edges(data=True))[0:10])
            print("Original graph edges:", list(protein_graphs_with_data[pdb_code][0].edges(data=True))[0:10])

            with open(f"{errors_path}/{dataset_name}_errors.log", "a") as f:
                f.write(f"Error in graph extraction for {pdb_code}: {e}\n")

            continue

    extract_time = sum(times_to_extract)/len(times_to_extract)
    print("Time to extract:", extract_time)
    

    print("\n\n")
    print("Number of edges in original graph: ", len(protein_graphs_with_data[pdb_code][0].edges()))
    print("Number of nodes in original graph: ", len(protein_graphs_with_data[pdb_code][0].nodes()))
    print("Number of edges in extracted graph: ", len(g.edges()))
    print("Number of nodes in extracted graph: ", len(g.nodes()))
    print("\n\n")
    print("original graph config:", protein_graphs_with_data[pdb_code][0].graph["config"].model_dump())
    print("extracted graph config:", g.graph["config"].model_dump())
    print("\n\n")
    print("uncompressed complete graph size", asizeof.asizeof(protein_graphs_with_data) / 1024 / 1024)
    print("uncompressed structure graph size", asizeof.asizeof(protein_graphs_without_data) / 1024 / 1024)
    print("uncompressed edge size", e_size)
    print("uncompressed node size", v_size)
    print("uncompressed complete graph serialized", asizeof.asizeof(pickle.dumps(protein_graphs_with_data)) / 1024 / 1024)
    print("uncompressed edge serialized", e_serialized)
    print("uncompressed node serialized", v_serialized)
    print("\n\n")
    print("compressed graph complete size", global_graph_obj.calculate_graph_complete_space_size())
    print("compressed complete node size", global_graph_obj.calculate_total_nodes_size())
    print("compressed complete edge size", global_graph_obj.calculate_total_edges_size())
    print("compressed node attributes size", global_graph_obj.node_attrs_size())
    print("compressed edge attributes size", global_graph_obj.edge_attrs_size())
    print("compressed node attributes keys size", global_graph_obj.node_attr_keys_size())
    print("compressed edge attributes keys size", global_graph_obj.edge_attr_keys_size())
    print("compressed pdb to nodes size", global_graph_obj.pdb_to_nodes_size())
    print("compressed pdb to edges size", global_graph_obj.pdb_to_edges_size())
    print("compressed node to id size", global_graph_obj.node_to_id_size())
    print("compressed edge to id size", global_graph_obj.edge_to_id_size())
    print("\n\n")
    print("compressed graph object size", asizeof.asizeof(global_graph_obj) / 1024 / 1024)
    print("compressed graph complete size serialized", asizeof.asizeof(pickle.dumps(global_graph_obj)) / 1024 / 1024)

    with open(f"{results_path}/{dataset_name}_results.txt", "a") as f:
        f.write(f"Time to construct graphs: {construct_time}\n")
        f.write(f"Time to compress: {compress_time}\n")
        f.write(f"Time to extract: {extract_time}\n")
        f.write(f"uncompressed complete graph size: {asizeof.asizeof(protein_graphs_with_data) / 1024 / 1024}\n")
        f.write(f"uncompressed structure graph size: {asizeof.asizeof(protein_graphs_without_data) / 1024 / 1024}\n")
        f.write(f"uncompressed edge size: {e_size}\n")
        f.write(f"uncompressed node size: {v_size}\n")
        f.write(f"uncompressed complete graph serialized: {asizeof.asizeof(pickle.dumps(protein_graphs_with_data)) / 1024 / 1024}\n")
        f.write(f"uncompressed edge serialized: {e_serialized}\n")
        f.write(f"uncompressed node serialized: {v_serialized}\n")
        f.write("\n\n")
        f.write(f"compressed graph complete size: {global_graph_obj.calculate_graph_complete_space_size()}\n")
        f.write(f"compressed complete node size: {global_graph_obj.calculate_total_nodes_size()}\n")
        f.write(f"compressed complete edge size: {global_graph_obj.calculate_total_edges_size()}\n")
        f.write(f"compressed node attributes size: {global_graph_obj.node_attrs_size()}\n")
        f.write(f"compressed edge attributes size: {global_graph_obj.edge_attrs_size()}\n")
        f.write(f"compressed node attributes keys size: {global_graph_obj.node_attr_keys_size()}\n")
        f.write(f"compressed edge attributes keys size: {global_graph_obj.edge_attr_keys_size()}\n")
        f.write(f"compressed pdb to nodes size: {global_graph_obj.pdb_to_nodes_size()}\n")
        f.write(f"compressed pdb to edges size: {global_graph_obj.pdb_to_edges_size()}\n")
        f.write(f"compressed node to id size: {global_graph_obj.node_to_id_size()}\n")
        f.write(f"compressed edge to id size: {global_graph_obj.edge_to_id_size()}\n")
        f.write("\n\n")
        f.write(f"compressed graph object size: {asizeof.asizeof(global_graph_obj) / 1024 / 1024}\n")
        f.write(f"compressed graph complete size serialized: {asizeof.asizeof(pickle.dumps(global_graph_obj)) / 1024/ 1024}")

def toy_example():
    pass




def test_compress():
    g1 = nx.Graph()
    g2 = nx.Graph()
    g3 = nx.Graph()
    g4 = nx.Graph()
    g1.add_node(1, chain_id="value1")
    g1.add_node(2, chain_id="value2")
    g1.add_edge(1, 2, kind=set(["covalent", "RING"]), distance=1.0)
    g2.add_node(3, chain_id="value2")
    g2.add_node(4, chain_id="value2")
    g2.add_edge(3, 4, kind=set(["TRIPLE", "aromatic"]), distance=2.0)
    g3.add_node(5, chain_id="value2")
    g3.add_node(6, chain_id="value1")
    g3.add_edge(5, 6, kind=set(["aromatic", "RING"]), distance=3.0)
    g4.add_node(7, chain_id="value1")
    g4.add_node(8, chain_id="value1")
    g4.add_edge(7, 8, kind=set(["TRIPLE", "hbond"]), distance=4.0)

    config1 = ProteinGraphConfig(
        granularity="N",
        edge_construction_functions=[add_atomic_edges, add_ring_status]
    )

    config2 = ProteinGraphConfig(
        granularity="N",
        edge_construction_functions=[add_aromatic_interactions, add_bond_order]
    )

    config3 = ProteinGraphConfig(
        granularity="N",
        edge_construction_functions=[add_ring_status, add_aromatic_interactions]
    )

    config4 = ProteinGraphConfig(
        granularity="N",
        edge_construction_functions=[add_bond_order, add_hydrogen_bond_interactions]
    )

    g1.graph["config"] = config1
    g2.graph["config"] = config2
    g3.graph["config"] = config3
    g4.graph["config"] = config4

    protein_graphs = {
        "1ABC": [g1, g2],
        "2DEF": [g3, g4]
    }

    print(g1.edges(data=True))

    compress_with_composition(protein_graphs)


if __name__ == "__main__":
    import os
    if os.environ["PROGRAM_NAME"] == "main":
        main()
    elif os.environ["PROGRAM_NAME"] == "toy":
        toy_example()









# all granularity opt:

# [
# ‘N’, ‘CA’, ‘C’, ‘O’, ‘CB’, ‘OG’, ‘CG’, ‘CD1’, ‘CD2’, ‘CE1’, ‘CE2’, ‘CZ’, ‘OD1’, 
# ‘ND2’, ‘CG1’, ‘CG2’, ‘CD’, ‘CE’, ‘NZ’, ‘OD2’, ‘OE1’, ‘NE2’, ‘OE2’, ‘OH’, ‘NE’, 
# ‘NH1’, ‘NH2’, ‘OG1’, ‘SD’, ‘ND1’, ‘SG’, ‘NE1’, ‘CE3’, ‘CZ2’, ‘CZ3’, ‘CH2’, ‘OXT’
# 'atom', 'centroids'
# ]

# testar apenas pra ca e atom
# medir o espaco de cada coisa serializada

