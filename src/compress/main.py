import os, metadata, time, pickle, traceback

import numpy as np
import pandas as pd
from pympler import asizeof

import Builder
from MemoryMeasuring import MemoryMeasuring
from PDBGraphStore import PDBGraphStore
import edge_functions_Model as edgeModel

from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
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

def create_dataset_result_file(result_path, experiment_fields):
    with open(f"{result_path}/results.csv", "w") as file:
        file.write(experiment_fields)

def initialize_pdb_data_path(general_data_path):
    pdb_data_path = os.path.abspath(f"{general_data_path}/pdb_files/")

    if not os.path.exists(pdb_data_path):
        os.makedirs(pdb_data_path)

    return pdb_data_path

def write_result(msg, result_path, file_mode='a'):
    with open(f"{result_path}/results.csv", file_mode) as file:
        file.write("\n")
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
    return [edgeModel.edge_functions_dict[f] for f in ["delaunay", "aromatic", "aromatic_sulphur"]]
    # , "aromatic_sulphur", "delaunay", "aromatic"

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

def measure_node_attributes_memory(graphs: dict):
    node_attrs = []
    node_attributes_memory = 0

    for _, graph_list in graphs.items():
        g = graph_list[0]
        for n in g.nodes:
            attrs = g.nodes[n]
            node_attr = {}
            for k, v in attrs.items():
                if isinstance(v, pd.Series):
                    node_attr[k] = tuple([tuple(v.tolist()), tuple(v.name)])
                elif isinstance(v, np.ndarray):
                    node_attr[k] = tuple(v)
                else:
                    node_attr[k] = v
                
                node_attrs.append(node_attr)

    node_attributes_memory = asizeof.asizeof(node_attrs)/1024/1024

    return node_attributes_memory

def measure_edge_attributes_memory(graphs: dict):
    edge_attributes_memory = 0
    edge_attrs = []

    for _, graph_list in graphs.items():
        g = graph_list[0]

        for e in g.edges:
            edge_attr = g.edges[e]

            edge_attrs.append(edge_attr)

    edge_attributes_memory += asizeof.asizeof(edge_attrs)/1024/1024

    return edge_attributes_memory

def experiment_1():

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
    number_of_nodes = 0
    number_of_edges = 0

    time_start = time.time()

    for _, pdb_code in enumerate(pdb_codes.copy()):
        try:
            graph_with_data, graph_without_data = prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path,1)
        except Exception as e:
            msg = traceback.format_exc()
            print(msg)
            write_error(dataset_name, msg, error_path)
            continue
        protein_graph_with_metadata_dict[pdb_code] = graph_with_data[pdb_code]
        protein_graph_without_metadata_dict[pdb_code] = graph_without_data[pdb_code]

        number_of_edges+=len(protein_graph_without_metadata_dict[pdb_code][-1].edges())
        number_of_nodes+=len(protein_graph_without_metadata_dict[pdb_code][-1].nodes())

    time_to_construct = time_count(time_start=time_start)
    
    result_columns = [
        "dataset"
        "Time to construct graphs",
        "number of nodes in dataset",
        "number of edges in dataset",
        "Number of graphs",
        "Time to compress",
        "Uncompressed complete graph size",
        "Uncompressed complete graph size serialized",
        "Uncompressed graph structure",
        "Compressed graph",
        "Compressed graph structure",
        "pdb_code_to_id",
        "pdb_code_to_id_serialized",
        "pdb_id_to_nodes",
        "pdb_id_to_nodes_serialized",
        "pdb_id_to_edges",
        "pdb_id_to_edges_serialized",
        "node_label_to_node_id",
        "node_label_to_node_id_serialized",
        "edge_label_to_edge_id",
        "edge_label_to_edge_id_serialized",
        "Uncompressed node attr values",
        "Uncompressed edge attr values",
        "Compressed dict attributes",
        "Total compressed node attr values",
        "Total compressed edge attr values",
        "attr_keys",
        "attr_keys_serialized",
        "edge_attr_values",
        "edge_attr_values_serialized",
        "node_attr_values",
        "node_attr_values_serialized",
        "Compressed node attributes",
        "node_global_attr_keyvalue_mapping",
        "node_global_attr_keyvalue_mapping_serialized",
        "node_local_attr_keyvalue_mapping",
        "node_local_attr_keyvalue_mapping_serialized",
        "Compressed edge attributes",
        "edge_local_attr_keyvalue_mapping",
        "edge_local_attr_keyvalue_mapping_serialized",
        "Compressed graph object size",
        "Compressed graph complete size serialized"
    ]

    result_line = []

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path, msg)
    elif os.getenv("EXCOUNT") == '-1':
        msg = ",".join(result_columns)
        write_result(msg=msg, result_path=result_path)

    result_line.append(dataset_name)
    result_line.append(time_to_construct)
    result_line.append(number_of_nodes)
    result_line.append(number_of_edges)
    result_line.append(len(protein_graph_with_metadata_dict))

    body_parts, time_to_compress = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)
    result_line.append(time_to_compress)

    pdb_store = PDBGraphStore(body_parts)
    memory = MemoryMeasuring(pdb_store)

    result_line.append(asizeof.asizeof(protein_graph_with_metadata_dict) /1024 / 1024)
    result_line.append(len(pickle.dumps(protein_graph_with_metadata_dict)) /1024 / 1024)
    result_line.append(asizeof.asizeof(protein_graph_without_metadata_dict)/1024/1024)
    result_line.append(memory.total_memory())
    result_line.append(memory.graph_structure_memory())
    result_line.append(memory.pdb_code_to_id_memory())
    result_line.append(memory.pdb_code_to_id_serialized_memory())
    result_line.append(memory.pdb_id_to_nodes_memory())
    result_line.append(memory.pdb_id_to_nodes_serialized_memory())
    result_line.append(memory.pdb_id_to_edges_memory())
    result_line.append(memory.pdb_id_to_edges_serialized_memory())
    result_line.append(memory.node_label_to_node_id_memory())
    result_line.append(memory.node_label_to_node_id_serialized_memory())
    result_line.append(memory.edge_label_to_edge_id_memory())
    result_line.append(memory.edge_label_to_edge_id_serialized_memory())
    result_line.append(measure_node_attributes_memory(protein_graph_with_metadata_dict))
    result_line.append(measure_edge_attributes_memory(protein_graph_with_metadata_dict))
    result_line.append(memory.dict_attributes_memory())
    result_line.append(memory.node_attr_values_memory() + memory.node_global_attr_keyvalue_mapping_memory() + memory.node_local_attr_keyvalue_mapping_memory())
    result_line.append(memory.edge_attr_values_memory() + memory.edge_local_attr_keyvalue_mapping_memory())
    result_line.append(memory.attr_keys_memory())
    result_line.append(memory.attr_keys_serialized_memory())
    result_line.append(memory.edge_attr_values_memory())
    result_line.append(memory.edge_attr_values_serialized_memory())
    result_line.append(memory.node_attr_values_memory())
    result_line.append(memory.node_attr_values_serialized_memory())
    result_line.append(memory.node_attributes_memory())
    result_line.append(memory.node_global_attr_keyvalue_mapping_memory())
    result_line.append(memory.node_global_attr_keyvalue_mapping_serialized_memory())
    result_line.append(memory.node_local_attr_keyvalue_mapping_memory())
    result_line.append(memory.node_local_attr_keyvalue_mapping_serialized_memory())
    result_line.append(memory.edge_attributes_memory())
    result_line.append(memory.edge_local_attr_keyvalue_mapping_memory())
    result_line.append(memory.edge_local_attr_keyvalue_mapping_serialized_memory())
    result_line.append(asizeof.asizeof(pdb_store)/1024/1024)
    result_line.append(len(pickle.dumps(pdb_store))/1024/1024)

    result_line = [f'{item:.2f}' if type(item) == float else f'{item}' for item in result_line]

    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path)
    # print(msg)

def experiment_2():
    pass

def experiment_3():
    pass

def experiment_4():
    pass

def experiment_5():
    pass

def experimen_6():
    pass

if __name__=="__main__":
    experiment_1()
