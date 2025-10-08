import Builder
from PDBGraphStore import PDBGraphStore
import os
import metadata
from bidict import bidict
import random
import time
import numpy as np
from pympler import asizeof
import pickle
import Builder
import edge_functions_Model as edgeModel
import operations
import networkx as nx
import traceback

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.graphs import compute_edges
from graphein.protein.utils import download_pdb

def initialize_errors_directory(current_file_path):
    error_path = os.path.abspath(f"{current_file_path}/../../errors/")

    if not os.path.exists(error_path):
        os.makedirs(error_path)

    return error_path

def create_dataset_error_file(error_path, dataset_name):
    print(f'error path: {error_path}')
    with open(f"{error_path}/{dataset_name}_errors.log", "w") as file:
        file.write(f"Errors log for {dataset_name}")

def initialize_results_directory(current_file_path):
    result_path = os.path.abspath(f"{current_file_path}/../../results/")

    if not os.path.exists(result_path):
        os.makedirs(result_path)

    return result_path

def create_dataset_result_file(result_path, dataset_name):
    with open(f"{result_path}/{dataset_name}_results.log", "w") as file:
        file.write(f"Result log for {dataset_name}\n")

def initialize_pdb_data_path(general_data_path):
    pdb_data_path = os.path.abspath(f"{general_data_path}/pdb_files/")

    if not os.path.exists(pdb_data_path):
        os.makedirs(pdb_data_path)

    return pdb_data_path

def write_result(dataset, msg, result_path, file_mode='a'):
    with open(f"{result_path}/{dataset}_results.log", file_mode) as file:
        file.write(msg)

def write_error(dataset, msg, error_path, file_mode='a'):
    with open(f"{error_path}/{dataset}_errors.log", file_mode) as file:
        file.write(msg)

def read_dataset(general_data_path, dataset_txt_name, file_mode='r'):
    pdb_codes = list()
    with open(f"{general_data_path}/{dataset_txt_name}", file_mode) as file:
        for line in file:
            if line[0] != '#':
                pdb_codes.append(line.strip().upper())

    return pdb_codes

def define_graphein_edge_funcs(func_idx=1):
    # usado no experimento 1
    # return random.sample([v for _, v in edgeModel.edge_functions_dict.items()], 3)
    return [edgeModel.edge_functions_dict[f] for f in ["aromatic", "bb_carbonyl_carbonyl", "delaunay"]] 
    
    # usado no experimento 2
    # sorted_func_list = sorted(edgeModel.edge_functions_dict.keys())
    # print(f'returning {sorted_func_list[func_idx]}')
    # print(f'returning {[edgeModel.edge_functions_dict[sorted_func_list[func_idx]]]}')
    # return [edgeModel.edge_functions_dict[sorted_func_list[func_idx]]]
    
    #usado no experimento 3
    # print(list(edgeModel.edge_functions_dict.values()))
    # return list(edgeModel.edge_functions_dict.values()) 


def define_configuration(edge_construction_funcs):
    return {
        "granularity": "CA",
        "edge_construction_functions": edge_construction_funcs
    }

def time_count(time_start):
    return time.time() - time_start

def get_pdb_file(pdb_data_path, pdb_code):
    pdb_code = pdb_code.lower()
    if os.path.exists(f"{pdb_data_path}/{pdb_code}.pdb"):
        print(f"Reading {pdb_code} from local directory")
        pdb_file = os.path.abspath(f"{pdb_data_path}/{pdb_code}.pdb")
    else:
        print(f"Downloading {pdb_code} from PDB")
        try:
            pdb_file = download_pdb(pdb_code, f"{pdb_data_path}/")
        except Exception as e:
            raise e
        
    if pdb_file == None:
        raise Exception("Error reading the pdb file")

    return pdb_file

def prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path, func_idx):
    protein_graph_with_metadata_dict = dict()
    protein_graph_without_metadata_dict = dict()

    print(func_idx)
    edge_funcs = define_graphein_edge_funcs(func_idx)
    config = ProteinGraphConfig(**define_configuration(edge_funcs))

    try:
        pdb_file = get_pdb_file(pdb_data_path=pdb_data_path, pdb_code=pdb_code)
    except Exception as e:     
        raise e    

    graph = construct_graph(config=config, path=pdb_file)
    graph.graph["pdb_code"] = pdb_code
    print(graph, "\n")

    graph_config = graph.graph["config"]
    graph_pdb_code = graph.graph["pdb_code"]
    graph.graph.clear()
    graph.graph["config"] = graph_config
    graph.graph["pdb_code"] = graph_pdb_code

    protein_graph_with_metadata_dict.setdefault(pdb_code, []).append(graph.copy())

    for node in graph.nodes():
        graph.nodes[node].clear()
        
    for u, v in graph.edges():
        graph.edges[u, v].clear()


    protein_graph_without_metadata_dict.setdefault(pdb_code, []).append(graph.copy())

    del graph

    return protein_graph_with_metadata_dict, protein_graph_without_metadata_dict

def test():
    # graphs_extracted = operations.extract_pdb_graphs_multiprocessing(pdb_store, ["1BXL", "1G5J"], ["bb_carbonyl_carbonyl"], 2)


    # for graphs in graphs_extracted:
    #     for g in graphs:
    #         print(g.graph["pdb_code"])
    #         print(g)
    #         for u, v in g.edges():
    #             print(g.edges[u, v])
    # #insert pdb into the pdb_store
    # graph_to_insert, _ = prepare_graph(pdb_data_path, "2NL9", dataset_name, error_path)

    # pdb_store.insert_pdbs(graph_to_insert)
    
    # #remove pdb from the pdb_store
    # # print(pdb_store.remove_multiple_pdbs(["2NL9"]))    

    # #merge 2 pdb_stores
    # pdb_store2 = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, node_attr_keys, edge_attr_keys)

    # #operation merge graph stores
    # gs_merged = operations.merge_graph_stores([pdb_store, pdb_store2])
    # print(pdb_store)
    # #operation split graph store into 2
    # gs1, gs2 = operations.split_graph_store(pdb_store, ["2NL9"])
    
    # print(gs1, gs2)
    pass

def experimento_1():

    #config usada:
    #granularity: CA
    #edge_construction_funcs: ["aromatic", "bb_carbonyl_carbonyl", "delaunay"]

    #datasets usados: todos (porem, alguns nao foram possiveis terminar a execucao por causa de estouro da memoria)
    
    current_file_path = os.path.dirname(os.path.realpath(metadata.__file__))
    general_data_path = os.environ.get("DATA_DIR") if os.environ.get("DATA_DIR") is not None else os.path.abspath(f"{current_file_path}/../../data/")
    dataset_txt_file_name = os.environ.get("DATASET")
    dataset_name = dataset_txt_file_name.split(".")[0]
    
    error_path = initialize_errors_directory(current_file_path=current_file_path)
    result_path = initialize_results_directory(current_file_path=current_file_path)
    pdb_data_path = initialize_pdb_data_path(general_data_path=general_data_path)

    create_dataset_error_file(error_path, dataset_name)
    create_dataset_result_file(result_path, dataset_name)

    print(f"\
          current_file_path={current_file_path}, \
          general_data_path={general_data_path}, \
          dataset_txt_file_name={dataset_txt_file_name}, \
          dataset_name={dataset_name}, \
          error_path={error_path}, \
          result_path={result_path}, \
          pdb_data_path={pdb_data_path} \
          ")

    pdb_codes = read_dataset(general_data_path=general_data_path, dataset_txt_name=dataset_txt_file_name)

    protein_graph_with_metadata_dict = {}
    protein_graph_without_metadata_dict = {}
    number_of_nodes_in_which_graph = list()
    number_of_edges_in_which_graph = list()

    time_start = time.time()

    for i, pdb_code in enumerate(pdb_codes.copy()):
        try:
            graph_with_data, graph_without_data = prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path,1)
        except Exception as e:
            msg = traceback.format_exc()
            print(msg)
            write_error(dataset_name, msg, error_path)
            continue
        protein_graph_with_metadata_dict[pdb_code] = graph_with_data[pdb_code]
        protein_graph_without_metadata_dict[pdb_code] = graph_without_data[pdb_code]

        number_of_edges_in_which_graph.append(len(protein_graph_without_metadata_dict[pdb_code][-1].edges()))
        number_of_nodes_in_which_graph.append(len(protein_graph_without_metadata_dict[pdb_code][-1].nodes()))

    msg = f'Average number of nodes: {np.mean(number_of_nodes_in_which_graph)} \
        \nAverage number of edges: {np.mean(number_of_edges_in_which_graph)} \
    '
    print(msg)

    write_result(dataset=dataset_name, msg=msg, result_path=result_path)        
        
    del number_of_edges_in_which_graph
    del number_of_nodes_in_which_graph
    
    time_to_construct = time_count(time_start=time_start)            

    msg = f'\nTime to construct graphs: {time_to_construct}\nNumber of graphs: {len(protein_graph_with_metadata_dict)}'

    write_result(dataset=dataset_name, msg=msg, result_path=result_path)

    v_size = 0
    e_size = 0
    v_serialized = 0
    e_serialized = 0

    for _, v in protein_graph_with_metadata_dict.items():
        v_size = np.sum([asizeof.asizeof(g._node)/1024/1024 for g in v])
        e_size = np.sum([asizeof.asizeof(g._adj)/1024/1024 for g in v])
        v_serialized = np.sum([asizeof.asizeof(pickle.dumps(g._node))/1024/1024 for g in v])
        e_serialized = np.sum([asizeof.asizeof(pickle.dumps(g._adj))/1024/1024 for g in v])


    time_start = time.time()
    print(v_size, v_serialized, e_size, e_serialized)

    node_to_id,\
    edge_to_id,\
    pdb_to_nodes,\
    pdb_to_edges,\
    node_attrs,\
    edge_kinds,\
    node_attr_keys,\
    edge_attr_keys, \
    edge_distances, \
    all_pdb_codes = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)

    time_to_compress = time_count(time_start=time_start)

    msg = f'\nTime to compress: {time_to_compress}'
    write_result(dataset=dataset_name, msg=msg, result_path=result_path)

    pdb_store = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_kinds, node_attr_keys, edge_attr_keys, edge_distances, all_pdb_codes)

    extract_times = []

    for pdb_code in set(protein_graph_with_metadata_dict.keys()):
        for g in protein_graph_with_metadata_dict[pdb_code]:
            time_start = time.time()
            extracted_graph = pdb_store.extract_pdb_graphs([pdb_code], [edgeModel.edge_functions_dict.inverse[func] for func in g.graph["config"].edge_construction_functions])
            time_to_extract = time_count(time_start=time_start)
            
            extract_times.append(time_to_extract)
            
            msg = f' \
                \nNumber of edges in original graph: {len(g.edges())} \
                \nNumber of nodes in original graph: {len(g.nodes())} \
                \nNumber of edges in extracted graph: {len(extracted_graph[0].edges())} \
                \nNumber of nodes in extracted graph: {len(extracted_graph[0].nodes())} \
                \n\n \
            '
            print(msg)

            try:
                assert nx.utils.nodes_equal(g._node, protein_graph_with_metadata_dict[pdb_code][0]._node)
                assert nx.utils.edges_equal(g._adj, protein_graph_with_metadata_dict[pdb_code][0]._adj)
            
            except AssertionError as e:
                msg = f'\n\
                Error in graph_extraction for {pdb_code}: {e}\
                '
                write_error(dataset=dataset_name, msg=msg, error_path=error_path)
                continue
    
    extract_time_mean = np.mean(extract_times)

    msg = f'\n\
        \nMean time to extract: {extract_time_mean} \
        \nUncompressed complete graph size: {asizeof.asizeof(protein_graph_with_metadata_dict) /1024 / 1024}\
        \nUncompressed structure graph size: {asizeof.asizeof(protein_graph_without_metadata_dict) /1024 / 1024}\
        \nUncompressed edge size: {e_size}\
        \nUncompressed node size: {v_size}\
        \nUncompressed complete graph serialized: {asizeof.asizeof(pickle.dumps(protein_graph_with_metadata_dict)) / 1024 / 1024}\
        \nUncompressed edge serialized: {e_serialized}\
        \nUncompressed node serialized: {v_serialized}\
        \n\
        \nCompressed graph complete size: {pdb_store.calculate_graph_complete_space_size()}\
        \nCompressed complete node size: {pdb_store.calculate_total_nodes_size()} \
        \nCompressed complete edge size: {pdb_store.calculate_total_edges_size()} \
        \nCompressed node attributes size: {pdb_store.node_attrs_size()}\
        \nCompressed edge attributes size: {pdb_store.edge_attrs_size()}\
        \nCompressed node attributes keys size: {pdb_store.node_attr_keys_size()}\
        \nCompressed edge attributes keys size: {pdb_store.edge_attr_keys_size()}\
        \nCompressed pdb to nodes size: {pdb_store.pdb_to_nodes_size()}\
        \nCompressed pdb to edges size: {pdb_store.pdb_to_edges_size()}\
        \nCompressed node to id size: {pdb_store.node_to_id_size()}\
        \nCompressed edge to id size: {pdb_store.edge_to_id_size()}\
        \n\
        \nCompressed graph object size: {asizeof.asizeof(pdb_store)/1024/1024}\
        \nCompressed graph complete size serialized: {asizeof.asizeof(pickle.dumps(pdb_store))/1024/1024}\
    '

    write_result(dataset=dataset_name, msg=msg, result_path=result_path)

def experimento_2():
    '''
    Avaliar a sobreposição de arestas entre as funções de arestas. 
    Por exemplo, será que existem funções de aresta que geram muitas 
    arestas iguais no grafo de uma determinada granularidade? Para isso, 
    podemos construir o PDB store incrementalmente, adicionando funções 
    de aresta uma a uma, e determinando se o tamanho do PDBStore está 
    aumentando significativamente quando incluímos uma nova função de aresta.

    obs: executado somente com o bcl_ppigremlin dataset

    get_len_edges
    get_len_nodes

    '''

    #config usada: 
    #granularity: atom
    #edge_construction_funcs: todas, menos fully_connected, pois ela adiciona todas as arestas possiveis entre os nodes
    
    #dataset usado: bcl (a principio apenas ele)

    current_file_path = os.path.dirname(os.path.realpath(metadata.__file__))
    general_data_path = os.environ.get("DATA_DIR") if os.environ.get("DATA_DIR") is not None else os.path.abspath(f"{current_file_path}/../../data/")
    dataset_txt_file_name = os.environ.get("DATASET")
    dataset_name = dataset_txt_file_name.split(".")[0]
    
    error_path = initialize_errors_directory(current_file_path=current_file_path)
    result_path = initialize_results_directory(current_file_path=current_file_path)
    pdb_data_path = initialize_pdb_data_path(general_data_path=general_data_path)

    create_dataset_error_file(error_path, dataset_name)
    create_dataset_result_file(result_path, dataset_name)

    print(f"\
          current_file_path={current_file_path}, \
          general_data_path={general_data_path}, \
          dataset_txt_file_name={dataset_txt_file_name}, \
          dataset_name={dataset_name}, \
          error_path={error_path}, \
          result_path={result_path}, \
          pdb_data_path={pdb_data_path} \
          ")

    pdb_codes = read_dataset(general_data_path=general_data_path, dataset_txt_name=dataset_txt_file_name)

    protein_graph_with_metadata_dict = {}
    number_of_nodes_in_which_graph = list()
    number_of_edges_in_which_graph = list()

    for _, pdb_code in enumerate(pdb_codes.copy()):
        try:
            #construct graph with 1 edge_func
            graph_with_data, _ = prepare_graph(pdb_data_path=pdb_data_path, pdb_code=pdb_code, dataset_name=dataset_name, error_path=error_path, func_idx=0)
        except:
            msg = traceback.format_exc()
            print(msg)
            write_error(dataset_name, msg, error_path)
            continue

        protein_graph_with_metadata_dict[pdb_code] = graph_with_data[pdb_code]

    node_to_id,\
    edge_to_id,\
    pdb_to_nodes,\
    pdb_to_edges,\
    node_attrs,\
    edge_attrs,\
    node_attr_keys,\
    edge_attr_keys = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)

    pdb_store = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, node_attr_keys, edge_attr_keys)

    print(f'graphStore with 1 edgefunc amount of edge: {pdb_store.get_len_edges()}')
    print(f'graphStore with 1 edgefunc amount of node: {pdb_store.get_len_nodes()}')

    msg = f'number of nodes: {pdb_store.get_len_nodes()}\nnumber of edges initially (with one func): {pdb_store.get_len_edges()}\
        \nmemory size initially: {asizeof.asizeof(pdb_store)/1024/1024:.2f} MB'
    write_result(dataset=dataset_name, msg=msg, result_path=result_path)

    for i in range(1, len(edgeModel.edge_functions_dict.keys())):
        graphs_to_insert = {}
        for pdb_code in pdb_codes:
            try:
                print(f'pdb_code: {pdb_code}')
                graph_to_insert, _ = prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path, i)
                
            except:
                msg = traceback.format_exc()
                print(msg)
                write_error(dataset_name, msg, error_path)
                continue
            graphs_to_insert[pdb_code] = graph_to_insert[pdb_code]

        sorted_func_list = sorted(edgeModel.edge_functions_dict.keys())
        pdb_store.insert_pdbs(graphs_to_insert)
        msg = f'\nNumber of edges after {i} insertion ({sorted_func_list[i]}): {pdb_store.get_len_edges()}\
            \nMemory size: {asizeof.asizeof(pdb_store)/1024/1024:.2f} MB'

        write_result(dataset=dataset_name, msg=msg, result_path=result_path)

def experimento_3():

    #config usada: 
    #granularity: atom
    #edge_construction_funcs: todas, menos fully_connected
    #dataset usado: bcl (a principio apenas ele)

    current_file_path = os.path.dirname(os.path.realpath(metadata.__file__))
    general_data_path = os.environ.get("DATA_DIR") if os.environ.get("DATA_DIR") is not None else os.path.abspath(f"{current_file_path}/../../data/")
    dataset_txt_file_name = os.environ.get("DATASET")
    dataset_name = dataset_txt_file_name.split(".")[0]
    
    error_path = initialize_errors_directory(current_file_path=current_file_path)
    result_path = initialize_results_directory(current_file_path=current_file_path)
    pdb_data_path = initialize_pdb_data_path(general_data_path=general_data_path)

    create_dataset_error_file(error_path, dataset_name)
    create_dataset_result_file(result_path, dataset_name)

    print(f"\
          current_file_path={current_file_path}, \
          general_data_path={general_data_path}, \
          dataset_txt_file_name={dataset_txt_file_name}, \
          dataset_name={dataset_name}, \
          error_path={error_path}, \
          result_path={result_path}, \
          pdb_data_path={pdb_data_path} \
          ")

    pdb_codes = read_dataset(general_data_path=general_data_path, dataset_txt_name=dataset_txt_file_name)

    protein_graph_with_metadata_dict = {}
    number_of_nodes_in_which_graph = list()
    number_of_edges_in_which_graph = list()

    for _, pdb_code in enumerate(pdb_codes.copy()):
        try:
            #construct graph with 1 edge_func
            graph_with_data, _ = prepare_graph(pdb_data_path=pdb_data_path, pdb_code=pdb_code, dataset_name=dataset_name, error_path=error_path, func_idx=0)
        except:
            msg = traceback.format_exc()
            print(msg)
            write_error(dataset_name, msg, error_path)
            continue

        protein_graph_with_metadata_dict[pdb_code] = graph_with_data[pdb_code]

    node_to_id,\
    edge_to_id,\
    pdb_to_nodes,\
    pdb_to_edges,\
    node_attrs,\
    edge_attrs,\
    node_attr_keys,\
    edge_attr_keys = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)

    pdb_store = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, node_attr_keys, edge_attr_keys)

    edge_funcs_to_extract = list(random.sample(edgeModel.edge_functions_dict.values(), 3))

    for _, pdb_code in enumerate(pdb_codes.copy()):
        time_start = time.time()
        extracted_graph = pdb_store.extract_pdb_graphs(pdb_codes=[pdb_code], edge_construction_functions=edge_funcs_to_extract)
        time_to_extract = time_count(time_start=time_start)
        print(extracted_graph[0].graph)

        time_start = time.time()
        config = ProteinGraphConfig(**{"granularity": "atom", "edge_construction_functions": edge_funcs_to_extract})
        pdb_file = get_pdb_file(pdb_data_path=pdb_data_path, pdb_code=pdb_code)
        graph = construct_graph(config=config, path=pdb_file)
        time_to_construct = time_count(time_start=time_start)

        msg = f'Time to construct pdb graph {pdb_code}: {time_to_construct}\
            \nTime to extract the same graph: {time_to_extract}\n\n'

        write_result(dataset=dataset_name, msg=msg, result_path=result_path)

if __name__=="__main__":
    experimento_1()