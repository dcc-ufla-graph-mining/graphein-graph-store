import Builder
from PDBGraphStore import PDBGraphStore
import os
import metadata
import random
import time
import numpy as np
from pympler import asizeof
import pickle
import Builder
import edge_functions_Model as edgeModel

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
    with open(f"{error_path}/{dataset_name}_error.log", "w") as file:
        file.write(f"Errors log for {dataset_name}")

def initialize_results_directory(current_file_path):
    result_path = os.path.abspath(f"{current_file_path}/../../results/")

    if not os.path.exists(result_path):
        os.makedirs(result_path)

    return result_path

def create_dataset_result_file(result_path, dataset_name):
    with open(f"{result_path}/{dataset_name}_error.log", "w") as file:
        file.write(f"Result log for {dataset_name}")

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

def define_graphein_edge_funcs():
    # return random.sample([v for _, v in edgeModel.edge_functions_dict.items()], 3)
    return [edgeModel.edge_functions_dict[f] for f in ["aromatic", "bb_carbonyl_carbonyl"]]

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
            return e
        
    if pdb_file == None:
        raise Exception("Error reading the pdb file")

    return pdb_file

def prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path):
    protein_graph_with_metadata_dict = dict()
    protein_graph_without_metadata_dict = dict()

    edge_funcs = define_graphein_edge_funcs()
    config = ProteinGraphConfig(**define_configuration(edge_funcs))

    try:
        pdb_file = get_pdb_file(pdb_data_path=pdb_data_path, pdb_code=pdb_code)
    except Exception as e:
        write_error(dataset_name, pdb_file, error_path)            

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

def main():
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
        graph_with_data,\
        graph_without_data = prepare_graph(pdb_data_path, pdb_code, dataset_name, error_path)
        protein_graph_with_metadata_dict[pdb_code] = graph_with_data[pdb_code]
        protein_graph_without_metadata_dict[pdb_code] = graph_with_data[pdb_code]

        number_of_edges_in_which_graph.append(len(protein_graph_without_metadata_dict[pdb_code][-1].edges()))
        number_of_nodes_in_which_graph.append(len(protein_graph_without_metadata_dict[pdb_code][-1].nodes()))

    msg = f"Average number of nodes: {np.mean(number_of_nodes_in_which_graph)} \
        \nAverage number of edges: {np.mean(number_of_edges_in_which_graph)} \
    "
    print(msg)

    write_result(dataset=dataset_name, msg=msg, result_path=result_path)        
        
    del number_of_edges_in_which_graph
    del number_of_nodes_in_which_graph
    
    time_to_construct = time_count(time_start=time_start)

    for k, v in protein_graph_with_metadata_dict.items():
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
    edge_attrs,\
    node_attr_keys,\
    edge_attr_keys = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)

    time_to_compress = time_count(time_start=time_start)

    pdb_store = PDBGraphStore(node_to_id, edge_to_id, pdb_to_nodes, pdb_to_edges, node_attrs, edge_attrs, node_attr_keys, edge_attr_keys)

    graphs_extracted = pdb_store.extract_pdb_graphs(["1BXL", "1G5J"], ["aromatic", "bb_carbonyl_carbonyl"])

    for g in graphs_extracted:
        print(g.graph["pdb_code"])
        print(g)
        for u, v in g.edges():
            print(g.edges[u, v])

    #TODO: insert pdb into the pdb_store
    # graph_to_insert, _ = prepare_graph(pdb_data_path, "2NL9", dataset_name, error_path)
    # pdb_store.insert_pdbs(graph_to_insert)
    #TODO: remove pdb from the pdb_store
    
    #TODO build the output string and then print it and write it to results file

if __name__=="__main__":
    main()