import pickle
import random
import time
import pandas as pd

import networkx as nx
import numpy as np
from bidict import bidict
from graphein.protein import add_atomic_edges
from pympler import asizeof
from pyroaring import BitMap, BitMap64
#the ordered_set lib has an data type named OrderedSet
#which is a set that keeps the order of the elements
#because of this, we can search for the index or for the key, both with o(1)
#https://pypi.org/project/ordered-set/
from ordered_set import OrderedSet
import metadata
import os

############################################################################################

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
#here happens the compression

# e = edge, v = vertex, g = graph

#contar o tempo so ate a parte de buildar o grafo grandao, excluindo o extract e o assert
def compress_with_composition(protein_graphs):
    edge_to_pdbs = {}   
    node_to_pdbs = {}  
    edge_attrs = {}   
    node_attrs = {}    

    edge_attr_keys_list = list(list(next(iter(protein_graphs.values())).edges(data=True))[0][2].keys())
    node_attr_keys_list = list(list(next(iter(protein_graphs.values())).nodes(data=True))[0][1].keys())

    edge_attr_keys = {}
    node_attr_keys = {}  

    for key in edge_attr_keys_list: 
        edge_attr_keys[key] = OrderedSet()
    for key in node_attr_keys_list: 
        node_attr_keys[key] = OrderedSet()

    del edge_attr_keys_list
    del node_attr_keys_list
    
    # Split pdb_to_view into four separate dictionaries
    pdb_to_nodes = {}
    pdb_to_edges = {}

    for pdb_code, g in protein_graphs.items():
        # Initialize separate dictionaries for each component
        pdb_to_nodes[pdb_code] = [BitMap64()]
        pdb_to_edges[pdb_code] = [BitMap64()]
        
        for node in g.nodes():
            if node not in node_to_pdbs:
                node_to_pdbs[node] = []

                attr_indexes = []
                
                for value in node_attr_keys:
                    attr_value = g.nodes[node][value]

                    if isinstance(attr_value, pd.Series):
                        attr_value = tuple([tuple(attr_value.tolist()), attr_value.name, tuple(attr_value.index)])
                    elif isinstance(attr_value, np.ndarray):
                        attr_value = tuple([tuple(attr_value)])

                    attr_indexes.append(node_attr_keys[value].add(attr_value))

                    node_attrs[node] = attr_indexes
                
            if pdb_code not in node_to_pdbs[node]:
                node_to_pdbs[node].append(pdb_code)
        
        for u, v, data in g.edges(data=True):
            edge = (u, v)
            
            if edge not in edge_to_pdbs:
                edge_to_pdbs[edge] = []
                attr_indexes = []
                
                for value in edge_attr_keys:
                    attr_value = data[value]
                    if isinstance(attr_value, np.ndarray):
                        attr_value = tuple(attr_value)
                    elif isinstance(attr_value, (list, set)):
                        attr_value = tuple(attr_value)
                        
                    attr_indexes.append(edge_attr_keys[value].add(attr_value))
                    edge_attrs[edge] = attr_indexes
            
            if pdb_code not in edge_to_pdbs[edge]:
                edge_to_pdbs[edge].append(pdb_code)

    edge_id = 0
    edge_to_id = {}
    for e in edge_to_pdbs:
        edge_to_id[e] = edge_id
        pdbs = edge_to_pdbs[e]
        for pdb_code in pdbs:
            pdb_to_edges[pdb_code][0].add(edge_id)
            
        edge_id += 1

    del edge_to_pdbs

    node_id = 0
    node_to_id = {}  
    for u in node_to_pdbs:
        node_to_id[u] = node_id
        pdbs = node_to_pdbs[u]
        for pdb_code in pdbs:
            pdb_to_nodes[pdb_code][0].add(node_id)
            
        node_id += 1

    del node_to_pdbs

    node_to_id = bidict(node_to_id)  
    edge_to_id = bidict(edge_to_id)

    for pdb_code in protein_graphs:
        original_graph = protein_graphs[pdb_code]

        nodes = [node_to_id.inverse[node_id] for node_id in pdb_to_nodes[pdb_code][0]]  
        edges = [edge_to_id.inverse[edge_id] for edge_id in pdb_to_edges[pdb_code][0]]
        extracted_graph = nx.Graph()
        extracted_graph.update(edges=edges, nodes=nodes)

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

        for u, v in edges:
            if (u, v) in edge_attrs:
                for i, key in enumerate(edge_attr_keys):
                    index = edge_attrs[(u, v)][i]
                    value = edge_attr_keys[key][index]
                    
                    if isinstance(value, tuple):
                        value = set(value)
                        
                    extracted_graph.edges[u, v][key] = value

        # assert nx.utils.nodes_equal(extracted_graph.nodes(data=True), original_graph.nodes(data=True))
        try:
            assert nx.utils.edges_equal(extracted_graph._adj, original_graph._adj)
        except Exception as e:
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
            break

        try:    
            assert nx.utils.nodes_equal(extracted_graph._node, original_graph._node)
        except Exception as e:
            print("error in nodes")
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
            break

            with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                f.write(f"Error in graph extraction for {pdb_code}: {e}\n")
    
    # Pass the four separate dictionaries to the constructor
    return PDBGraphStoreBitmap(node_to_id, edge_to_id, 
                            pdb_to_nodes, pdb_to_edges, 
                            node_attrs, edge_attrs, 
                            edge_attr_keys, node_attr_keys)


class PDBGraphStoreBitmap:

    #fazer o calculo de memoria discriminado para cada um dos attr abaixo, incluindo acada uma das posicoes da view separado
    def __init__(self, node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, edge_attr_keys, node_attr_keys):
        self.node_to_id = node_to_id #mapeamento de de node para id, e vice versa, global
        self.edge_to_id = edge_to_id #mapeamento de de edge para id, e vice versa, global
        self.pdb_to_nodes = pdb_to_nodes #bitmap indicando quais pdbs cada node pertence
        self.pdb_to_edges = pdb_to_edges #bitmap indicando quais pdbs cada edge pertence
        self.node_attrs = node_attrs #lista de indices para o attr de cada node
        self.edge_attrs = edge_attrs #lista de indices para o attr de cada edge
        self.node_attr_keys = node_attr_keys #dicionario de atributos indexados para cada node
        self.edge_attr_keys = edge_attr_keys #dicionario de atributos indexados para cada edge

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
    

    def extract_pdb_graph(self, pdb_code):
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
            
            for u, v in edges:
                if (u, v) in self.edge_attrs:
                    for i, key in enumerate(self.edge_attr_keys):
                        index = self.edge_attrs[(u, v)][i]
                        value = self.edge_attr_keys[key][index]
                        
                        if isinstance(value, tuple):
                            value = np.array(value)  
                            
                        extracted_graph.edges[u, v][key] = value
        except Exception as e:
            print(f"Error extracting graph for {pdb_code}: {e}")
            extracted_graph = nx.Graph()
            extracted_graph.update(edges=[], nodes=[])
            with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                f.write(f"Error extracting graph for {pdb_code}: {e}\n")
        return extracted_graph

    # def run_tree_compression(self):
    #     # build pairs jaccard
    #     def get_max_pair(pdb_to_bitmaps):
    #         best_pair = None
    #         best_jacc = None
    #         for pdb_code1, bms1 in pdb_to_bitmaps.items():
    #             bm_1 = bms1[0]
    #             for pdb_code2, bms2 in pdb_to_bitmaps.items():
    #                 if pdb_code1 >= pdb_code2: continue
    #                 bm_2 = bms2[0]
    #                 jacc = bm_1.intersection_cardinality(bm_2)
    #                 if best_jacc is None or jacc > best_jacc:
    #                     best_jacc = jacc
    #                     best_pair = (pdb_code1, pdb_code2)
    #         return best_pair

    #     def get_optimized_bitmaps(reconstructions):
    #         while len(reconstructions) >= 2:
    #             pdb_code1, pdb_code2 = get_max_pair(reconstructions)
    #             # print(pdb_code1, pdb_code2)
    #             left = reconstructions[pdb_code1][0]
    #             right = reconstructions[pdb_code2][0]
    #             intersection = left & right
    #             delta_1 = left - intersection
    #             delta_2 = right - intersection

    #             left_child = (delta_1, reconstructions[pdb_code1][1])
    #             right_child = (delta_2, reconstructions[pdb_code2][1])
    #             reconstructions[pdb_code1 + '_' + pdb_code2] = (
    #                 intersection, {pdb_code1: left_child, pdb_code2: right_child})

    #             del reconstructions[pdb_code1]
    #             del reconstructions[pdb_code2]

    #         def key_in_keys(key, keys):
    #             for k in keys:
    #                 if key in k: return k
    #             return None

    #         root_tree_node = list(reconstructions.items())[0]

    #         pdb_to_view_opt = {}

    #         for pdb_code in self.pdb_to_view:
    #             tree_node = root_tree_node
    #             bms = [tree_node[1][0]] if len(tree_node[1][0]) > 0 else []
    #             bm_height = 1
    #             key = key_in_keys(pdb_code, tree_node[1][1])
    #             while key:
    #                 bm_height += 1
    #                 tree_node = (key, tree_node[1][1][key])
    #                 if len(tree_node[1][0]) > 0:
    #                     tree_node[1][0].run_optimize()
    #                     bms.append(tree_node[1][0])
    #                 key = key_in_keys(pdb_code, tree_node[1][1])
    #             if len(bms) > 0:
    #                 pdb_to_view_opt[pdb_code] = bms

    #         return pdb_to_view_opt

    #     edge_bitmaps = get_optimized_bitmaps(
    #         {pdb_code: (view[1][0], {}) for pdb_code, view in self.pdb_to_view.items()})
    #     isolated_node_bitmaps = get_optimized_bitmaps(
    #         {pdb_code: (view[0][0], {}) for pdb_code, view in self.pdb_to_view.items()})

    #     # TODO: assertion code (remove) {
    #     for pdb_code, view in self.pdb_to_view.items():
    #         union_bm1 = BitMap64()
    #         for bm in isolated_node_bitmaps.get(pdb_code, BitMap64()):
    #             union_bm1 = union_bm1 | bm

    #         assert union_bm1.symmetric_difference_cardinality(view[0][0]) == 0

    #         union_bm1 = BitMap64()
    #         for bm in edge_bitmaps[pdb_code]:
    #             union_bm1 = union_bm1 | bm

    #         assert union_bm1.symmetric_difference_cardinality(view[1][0]) == 0
    #     # }

    #     # replace existing views with optimized ones
    #     self.pdb_to_view = {pdb_code: (isolated_node_bitmaps.get(pdb_code, None), edge_bitmaps.get(pdb_code, None)) for
    #                         pdb_code in self.pdb_to_view}
    #     self.pdb_to_view = {k: v for k, v in self.pdb_to_view.items() if (k, v) != (None, None)}


def main():
    import os
    import metadata
    from graphein.protein.config import ProteinGraphConfig
    from graphein.protein.edges.distance import add_hydrogen_bond_interactions, add_peptide_bonds
    from graphein.protein.graphs import construct_graph
    from graphein.protein.utils import download_pdb
    import networkx as nx
    import time

    print("main")

    pdb_codes = []

    params_to_change = {"granularity": "N", "edge_construction_functions": [add_atomic_edges, add_hydrogen_bond_interactions, add_peptide_bonds]}

    config = ProteinGraphConfig(**params_to_change)
    # print(config.model_dump())
    
    with open(f'{data}/{dataset}', 'r') as f:
        for line in f:
            pdb_codes.append(line.strip())

    protein_graphs_with_data = {}
    protein_graphs_without_data = {}

    time_begin = time.time()

    i = 0

    pdb_codes_copy = pdb_codes.copy()

    for pdb_code in pdb_codes_copy:
        print(i, pdb_code)
        i += 1
        
        if os.path.exists(f"{pdb_dir}/{pdb_code}.pdb"):
            try:
                pdb_file = os.path.abspath(f"{pdb_dir}/{pdb_code}.pdb")
            except Exception as e:
                print(f"Error reading {pdb_code}: {e}")
                pdb_codes.remove(pdb_code)
                with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                    f.write(f"Error reading {pdb_code}: {e}\n")
                continue
        else:
            try:
                pdb_file = download_pdb(pdb_code, f"{pdb_dir}/")
                if pdb_file is None:
                    print(f"Failed to download {pdb_code}")
                    pdb_codes.remove(pdb_code)
                    with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                        f.write(f"Failed to download {pdb_code}\n")
                    continue
            except Exception as e:
                print(f"Error downloading {pdb_code}: {e}")
                pdb_codes.remove(pdb_code)

                with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                    f.write(f"Error downloading {pdb_code}: {e}\n")
                continue

        graph = construct_graph(config=config, path=pdb_file)


        graph.graph.clear()
        print(graph)

        protein_graphs_with_data[pdb_code] = graph.copy()  # Store graph

        for node in graph.nodes():
            graph.nodes[node].clear()
        for u, v in graph.edges():
            graph.edges[u, v].clear()

        protein_graphs_without_data[pdb_code] = graph.copy()  # Store graph

        del graph

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
        v.graph.clear()
        v_size += asizeof.asizeof(v._node) / 1024 / 1024
        e_size += asizeof.asizeof(v._adj) / 1024 / 1024

        v_serialized += asizeof.asizeof(pickle.dumps(v._node)) / 1024 / 1024
        e_serialized += asizeof.asizeof(pickle.dumps(v._adj)) / 1024 / 1024

    time_begin = time.time()

    global_graph_obj = compress_with_composition(protein_graphs_with_data)

    time_end = time.time()
    compress_time = time_end - time_begin
    print("Time to compress:", compress_time)

    random.shuffle(pdb_codes)

    times_to_extract = []

    for pdb_code in pdb_codes:
        time_begin = time.time()
        g = global_graph_obj.extract_pdb_graph(pdb_code)
        time_end = time.time()
        times_to_extract.append(time_end - time_begin)

        try:
            assert nx.utils.nodes_equal(g._node, protein_graphs_with_data[pdb_code]._node) 
            assert nx.utils.edges_equal(g._adj, protein_graphs_with_data[pdb_code]._adj)
        except AssertionError as e:
            print(f"Error in graph extraction for {pdb_code}: {e}")
            print("Extracted graph nodes:", g.nodes(data=True))
            print("Original graph nodes:", protein_graphs_with_data[pdb_code].nodes(data=True))
            print("Extracted graph edges:", g.edges(data=True))
            print("Original graph edges:", protein_graphs_with_data[pdb_code].edges(data=True))
            with open(f"{errors_path}/{dataset}_errors.log", "a") as f:
                f.write(f"Error in graph extraction for {pdb_code}: {e}\n")

            continue

    extract_time = sum(times_to_extract)/len(times_to_extract)
    print("Time to extract:", extract_time)
    

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

    with open(f"{results_path}/{dataset}_results.txt", "w") as f:
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


#passo 1: discriminar o view     -- v
#passo 2: toy exemple   //fazer no draw io tendo uma visao em formato de desenho do objeto final
#passo 3: rodar experimentos grandes
