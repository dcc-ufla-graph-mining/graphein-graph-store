import os, metadata, time, traceback, random

import numpy as np
import pandas as pd
from pympler import asizeof
import networkx as nx
import pickle as pk

import Builder
from MemoryMeasuring import MemoryMeasuring
from PDBGraphStore import PDBGraphStore
import edge_functions_Model as edgeModel
from operations import *

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

def create_dataset_result_file(result_path, experiment_fields, func):
    with open(f"{result_path}/results_{func}.csv", "w") as file:
        file.write(experiment_fields)

def initialize_pdb_data_path(general_data_path):
    pdb_data_path = os.path.abspath(f"{general_data_path}/pdb_files/")

    if not os.path.exists(pdb_data_path):
        os.makedirs(pdb_data_path)

    return pdb_data_path

def write_result(msg="", result_path="./", file_mode='a', func=""):
    with open(f"{result_path}/results_{func}.csv", file_mode) as file:
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

def build_graph(return_list_of_graphs=False):
    #config usada:
    #granularity: CA
    #edge_construction_funcs: ["delaunay", "aromatic", "aromatic_sulphur"]
    
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
    number_of_edges = 0
    number_of_nodes = 0

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

    if return_list_of_graphs:
        return protein_graph_with_metadata_dict

    body_parts, time_to_compress = Builder.compress_pdb_graphs(protein_graph_with_metadata_dict)
    pdb_store = PDBGraphStore(body_parts)

    exp_1_misc = {
        "time_to_construct": time_to_construct,
        "time_to_compress": time_to_compress,
        "number_of_nodes": number_of_nodes,
        "number_of_edges": number_of_edges,
        "protein_graph_with_data": protein_graph_with_metadata_dict,
        "protein_graph_without_data": protein_graph_without_metadata_dict
    }

    return exp_1_misc, pdb_store

current_file_path = os.path.dirname(os.path.realpath(metadata.__file__))
general_data_path = os.environ.get("DATA_DIR") if os.environ.get("DATA_DIR") is not None else os.path.abspath(f"{current_file_path}/../../data/")
dataset_txt_file_name = os.environ.get("DATASET")
dataset_name = dataset_txt_file_name.split(".")[0]

error_path = initialize_errors_directory(current_file_path=current_file_path)
result_path = initialize_results_directory(current_file_path=current_file_path)
pdb_data_path = initialize_pdb_data_path(general_data_path=general_data_path)

create_dataset_error_file(error_path, dataset_name)

def experiment_1(misc, pdb_store):
    result_columns = [
        "dataset",
        "Time to construct graphs",
        "number of nodes in dataset",
        "number of edges in dataset",
        "Number of graphs",
        "Time to compress",
        "Uncompressed complete graph size", #g
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
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_1.__name__)

    result_line.append(dataset_name)
    result_line.append(misc["time_to_construct"])
    result_line.append(misc["number_of_nodes"])
    result_line.append(misc["number_of_edges"])
    result_line.append(len(misc["protein_graph_with_data"]))

    result_line.append(misc["time_to_compress"])
    
    memory = MemoryMeasuring(pdb_store)

    result_line.append(asizeof.asizeof(misc["protein_graph_with_data"]) /1024 / 1024)
    result_line.append(len(pickle.dumps(misc["protein_graph_with_data"])) /1024 / 1024)
    result_line.append(asizeof.asizeof(misc["protein_graph_without_data"])/1024/1024)
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
    result_line.append(measure_node_attributes_memory(misc["protein_graph_with_data"]))
    result_line.append(measure_edge_attributes_memory(misc["protein_graph_with_data"]))
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

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_1.__name__)

def experiment_2(pdb_store):
    result_columns = [
        "dataset",
        "avg_time_to_extract"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_2.__name__)

    pdb_codes = pdb_store.get_this_pdb_list()

    times_to_extract = []
    for code in pdb_codes:
        time_start = time.time()
        _ = pdb_store.extract_pdb(code)
        times_to_extract.append(time_count(time_start=time_start))
    
    result_line.append(f'{np.mean(times_to_extract):.2f}')

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_2.__name__)

def experiment_3(pdb_store):
    result_columns = [
        "dataset",
        "avg_time_to_extract_1",
        "avg_time_to_extract_2",
        "avg_time_to_extract_4",
        "avg_time_to_extract_8",
        "avg_time_to_extract_16",
        "avg_time_to_extract_32"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_3.__name__)
    
    pdb_codes = pdb_store.get_this_pdb_list()
    
    i = 1

    for _ in range(6):
        time_start = time.time()
        pdb_extracted = extract_pdb_graphs_multiprocessing(pdb_store=pdb_store, pdb_codes=pdb_codes,num_cpus=i)

        result_line.append(f'{time_count(time_start=time_start)/len(pdb_codes):.2f}')
        i *= 2

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_3.__name__)

def experiment_4(pdb_store):
    result_columns = [
        "dataset",
        "avg_time_to_remove"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_4.__name__)

    pdb_codes = list(pdb_store.get_this_pdb_list())
    print(f'{type(pdb_codes)}, pdb_codes: {pdb_codes}')
    pdb_to_remove = random.choice(pdb_codes)


    time_start = time.time()
    pdb_store = remove_graph_from_store([pdb_to_remove], pdb_store)
    time_to_remove = time_count(time_start=time_start)
    
    result_line.append(f'{time_to_remove:.2f}')

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_4.__name__)

def experiment_5(pdb_store):
    result_columns = [
        "dataset",
        "time_to_split"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_5.__name__)

    time_start = time.time()
    pdb_store = split_graph_store(pdb_store=pdb_store, pdb_code_list=None)
    time_to_split = time_count(time_start=time_start)
    
    result_line.append(f'{time_to_split:.2f}')

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_5.__name__)

def experiment_6():
    result_columns = [
        "dataset",
        "time_to_merge"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_6.__name__)

    graphs = build_graph(True)

    mid = len(graphs) // 2

    keys = list(graphs.keys())

    body_parts, _ = Builder.compress_pdb_graphs({k: graphs[k] for k in keys[:mid]})
    store_1 =  PDBGraphStore(body_parts=body_parts)

    body_parts, _ = Builder.compress_pdb_graphs({k: graphs[k] for k in keys[mid:]})
    store_2 = PDBGraphStore(body_parts=body_parts)


    print(store_1)
    print(store_2)

    time_start = time.time()
    super_store = merge_graph_stores([store_1, store_2])
    time_to_merge = time_count(time_start=time_start)
    
    result_line.append(f'{time_to_merge:.2f}')

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_6.__name__)

    print(super_store)

def experiment_7(protein_graphs, pdb_store):
    with open("v1.pkl", "wb") as f:
        pk.dump(protein_graphs, f)

    with open("v2.pkl", "wb") as f:
        pk.dump(pdb_store, f)

def experiment_8():
    pass

def print_(pdb_graphs, pdb_store=None):
    for pdb_code, graph_list in pdb_graphs.items():
        for n in graph_list[0].nodes():
            for k, v in graph_list[0].nodes[n].items():
                print(k, v)
            print("\n\n")

def extract_max_min_graphs(protein_graphs):
    qtd_edges = {}
    qtd_nodes = {}

    for pdb_code, graphs in protein_graphs.items():
        graph = graphs[0]

        qtd_edges[pdb_code] = len(graph.edges)
        qtd_nodes[pdb_code] = len(graph.nodes)

    max_edge_key = max(qtd_edges, key=qtd_edges.get)
    max_node_key = max(qtd_nodes, key=qtd_nodes.get)

    min_edge_key = min(qtd_edges, key=qtd_edges.get)
    min_node_key = min(qtd_nodes, key=qtd_nodes.get)

    with open(f"min_max_results/{dataset_name}.csv", "a") as f:
        f.write(f"\n{qtd_edges[min_edge_key]},{qtd_edges[max_edge_key]},{qtd_nodes[min_node_key]},{qtd_nodes[max_node_key]}")

def toy_exemple():
    def make_meiler(residue_name: str, values: []) -> pd.Series:
        return pd.Series(
            values,
            index=[f"dim_{i}" for i in range(1, 8)],
            name=residue_name,
            dtype="float64",
        )


    ala_meiler = make_meiler("ALA", [1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23])
    gln_meiler = make_meiler("GLN", [1.56, 0.18, 3.95, -0.22, 5.65, 0.35, 0.25])
    glu_meiler = make_meiler("GLU", [1.56, 0.15, 3.78, -0.64, 3.09, 0.42, 0.21])


    G_1bxl = nx.Graph()

    G_1bxl.add_node(
        "A:ALA:89:N",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="N",
        element_symbol="N",
        coords=np.array([8.240, 5.876, -2.816], dtype=np.float32),
        b_factor=0.3700000047683716,
        meiler=ala_meiler,
    )

    G_1bxl.add_node(
        "A:ALA:89:CA",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="CA",
        element_symbol="C",
        coords=np.array([8.590, 4.452, -2.554], dtype=np.float32),
        b_factor=0.3700000047683716,
        meiler=ala_meiler,
    )

    G_1bxl.add_node(
        "A:GLN:88:C",
        chain_id="A",
        residue_name="GLN",
        residue_number=88,
        atom_type="C",
        element_symbol="C",
        coords=np.array([7.761, 6.236, -3.975], dtype=np.float32),
        b_factor=0.3799999952316284,
        meiler=gln_meiler,
    )

    G_1bxl.add_node(
        "A:ALA:89:C",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="C",
        element_symbol="C",
        coords=np.array([7.328, 3.582, -2.581], dtype=np.float32),
        b_factor=0.3499999940395355,
        meiler=ala_meiler,
    )

    G_1bxl.add_node(
        "A:ALA:89:CB",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="CB",
        element_symbol="C",
        coords=np.array([9.253, 4.341, -1.180], dtype=np.float32),
        b_factor=0.3700000047683716,
        meiler=ala_meiler,
    )

    G_1bxl.add_edge(
        "A:ALA:89:N",
        "A:ALA:89:CA",
        kind={"covalent"},
        bond_length=1.4896038743929203,
        distance=1.4896038743929203,
    )

    G_1bxl.add_edge(
        "A:ALA:89:N",
        "A:GLN:88:C",
        kind={"covalent"},
        bond_length=1.3047304477692825,
        distance=1.3047304477692825,
    )

    G_1bxl.add_edge(
        "A:ALA:89:CA",
        "A:ALA:89:C",
        kind={"covalent"},
        bond_length=1.533060154638553,
        distance=1.533060154638553,
    )

    G_1bxl.add_edge(
        "A:ALA:89:CA",
        "A:ALA:89:CB",
        kind={"covalent"},
        bond_length=1.5296293756226533,
        distance=1.5296293756226533,
    )


    G_1g5j = nx.Graph()

    G_1g5j.add_node(
        "A:ALA:89:N",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="N",
        element_symbol="N",
        coords=np.array([-16.303, -7.686, -10.552], dtype=np.float32),
        b_factor=0.4399999976158142,
        meiler=ala_meiler,
    )

    G_1g5j.add_node(
        "A:ALA:89:CA",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="CA",
        element_symbol="C",
        coords=np.array([-16.895, -6.506, -9.935], dtype=np.float32),
        b_factor=0.4000000059604645,
        meiler=ala_meiler,
    )

    G_1g5j.add_node(
        "A:GLU:48:C",
        chain_id="A",
        residue_name="GLU",
        residue_number=48,
        atom_type="C",
        element_symbol="C",
        coords=np.array([-15.052, -7.681, -11.000], dtype=np.float32),
        b_factor=0.46000000834465027,
        meiler=glu_meiler,
    )

    G_1g5j.add_node(
        "A:ALA:89:C",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="C",
        element_symbol="C",
        coords=np.array([-16.100, -6.084, -8.705], dtype=np.float32),
        b_factor=0.36000001430511475,
        meiler=ala_meiler,
    )

    G_1g5j.add_node(
        "A:ALA:89:CB",
        chain_id="A",
        residue_name="ALA",
        residue_number=89,
        atom_type="CB",
        element_symbol="C",
        coords=np.array([-18.345, -6.776, -9.563], dtype=np.float32),
        b_factor=0.44999998807907104,
        meiler=ala_meiler,
    )

    G_1g5j.add_edge(
        "A:ALA:89:N",
        "A:ALA:89:CA",
        kind={"covalent"},
        bond_length=1.457241665525958,
        distance=1.457241665525958,
    )

    G_1g5j.add_edge(
        "A:ALA:89:N",
        "A:GLU:48:C",
        kind={"covalent"},
        bond_length=1.3288072023326782,
        distance=1.3288072023326782,
    )

    G_1g5j.add_edge(
        "A:ALA:89:CA",
        "A:ALA:89:C",
        kind={"covalent"},
        bond_length=1.5241424747361823,
        distance=1.5241424747361823,
    )

    G_1g5j.add_edge(
        "A:ALA:89:CA",
        "A:ALA:89:CB",
        kind={"covalent"},
        bond_length=1.5211118260851224,
        distance=1.5211118260851224,
    )


    pdb_graphs = {
        "1bxl": [G_1bxl],
        "1g5j": [G_1g5j],
    }

    # print(pdb_graphs)

    # for pdb_code, graphs in pdb_graphs.items():
    #     graph = graphs[0]
    #     print(pdb_code)
    #     for n in graph.nodes():
    #         print(graph.nodes[n]['chain_id'])
    #         # 
    #     for e in graph.edges():
    #         print(graph.edges[e])

    body_parts, _ = Builder.compress_pdb_graphs(pdb_graphs)
    pdb_store = PDBGraphStore(body_parts)

    pdb_store.print_attr()

def experiment_o1(pdb_graphs):
    result_columns = [
        "dataset",
        "avg_time_to_extract"
    ]

    result_line = []

    result_line.append(dataset_name)

    if os.getenv("EXCOUNT") == '0':
        msg = ",".join(result_columns)
        create_dataset_result_file(result_path=result_path,experiment_fields=msg, func=experiment_o1.__name__)
    
    pdb_codes = pdb_graphs.keys()

    times_to_extract = []
    for code in pdb_codes:
        time_start = time.time()
        _ = pdb_graphs[code]
        times_to_extract.append(time_count(time_start=time_start))
    
    result_line.append(f'{np.mean(times_to_extract):.10f}')

    print(result_line)
    msg = ",".join(result_line)

    write_result(msg=msg, result_path=result_path, file_mode='a', func=experiment_o1.__name__)

if __name__=="__main__":
    exp_1_misc, pdb_store = build_graph()
    
    # pdb_graphs = build_graph(True)
    # experiment_o1(pdb_graphs)

    # experiment_1(exp_1_misc, pdb_store)
    # del exp_1_misc
    # experiment_2(pdb_store)
    # experiment_3(pdb_store)
    # experiment_4(pdb_store)
    # experiment_5(pdb_store)
    # experiment_6()
    experiment_7(exp_1_misc["protein_graph_with_data"], pdb_store)
    # experiment_8()
    # print_(pdb_graphs)
    # with open(f"min_max_results/{dataset_name}.csv", "w") as f:
    #     f.write("min_edges,max_edges,min_nodes,max_nodes")
    # extract_max_min_graphs(exp_1_misc["protein_graph_with_data"])

    # toy_exemple()