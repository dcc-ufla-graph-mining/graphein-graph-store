from typing import Dict

import networkx as nx
import numpy as np

from subdue import Subdue, Graph

import bisect

from pyroaring import BitMap


def compress_with_subdue(protein_graphs, **subdue_params):
    union_graph, pdb_code_to_id, graph_offsets = prepare_graphs(protein_graphs)

    #id_to_pdb_code = dict()
    #for pdb_code in pdb_code_to_id:
    #    gid = pdb_code_to_id[pdb_code]
    #    id_to_pdb_code[gid] = pdb_code

    #for u in union_graph.nodes:
    #    print(u, union_graph.nodes[u])

    #print(list(enumerate(graph_offsets)))
    reduced_graph, subdue_patterns = run_subdue(union_graph, **subdue_params)
    patterns = []
    graph_ids_set = []
    for iteration in range(len(subdue_patterns)):
        for pattern_id in range(len(subdue_patterns[iteration])):
            print(f"Iteration={iteration}, PatternId={pattern_id}")
            subgraphs = subdue_patterns[iteration][pattern_id]

            # get pattern format
            pattern = nx.Graph()
            for e in subgraphs[0]['edges']:
                u = union_graph.nodes[e[0]]['label']
                v = union_graph.nodes[e[1]]['label']
                pattern.add_edge(u, v)
            patterns.append(pattern)

            print(pattern.edges)

            #subgraph_patterns = [get_pattern_from_subgraph(s, union_graph) for s in subgraphs]
            #for p1 in subgraph_patterns:
            #    for p2 in subgraph_patterns:
            #        assert nx.utils.graphs_equal(p1, p2), f'{p1} {p2} {pattern}'

            # get a list of graph IDs the pattern occurs
            graph_ids = BitMap()
            for subgraph in subgraphs:
                sample_vertex = int(subgraph['nodes'][0])
                graph_id = bisect.bisect_right(graph_offsets, sample_vertex)
                #print(graph_id, id_to_pdb_code[graph_id], subgraph['nodes'])
                assert graph_id not in graph_ids, f"gid={graph_id} gids={graph_ids} p={pattern.edges}"
                graph_ids.add(graph_id)
            graph_ids_set.append(graph_ids)

    reduced_graph = subdue_to_nxgraph(reduced_graph)

    # some pattern bitmaps may be exactly the same, in this case, merge such pattern groups
    graph_bitmap_to_pattern_ids = {}
    for pattern_id in range(len(patterns)):
        k = frozenset(graph_ids_set[pattern_id])
        pattern_ids = graph_bitmap_to_pattern_ids.get(k, [])
        pattern_ids.append(pattern_id)
        graph_bitmap_to_pattern_ids[k] = pattern_ids

    pattern_groups = [v for v in graph_bitmap_to_pattern_ids.values() if len(v) > 1]
    del graph_bitmap_to_pattern_ids

    for pattern_group in pattern_groups:
        print(f"Merging group {pattern_group}")
        target_id = pattern_group[0]
        target_pattern = patterns[target_id]

        for pattern_id in pattern_group[1:]:
            other_pattern = patterns[pattern_id]
            target_pattern.update(nodes=other_pattern.nodes, edges=other_pattern.edges)
            patterns[pattern_id] = None
            graph_ids_set[pattern_id] = None

    print(f"NumPatternsBeforeFilter {len(patterns)} {len(graph_ids_set)}")
    patterns = [p for p in patterns if p is not None]
    graph_ids_set = [bm for bm in graph_ids_set if bm is not None]
    print(f"NumPatternsAfterFilter {len(patterns)} {len(graph_ids_set)}")

    return pdb_code_to_id, reduced_graph, patterns, graph_ids_set, graph_offsets


def subdue_to_nxgraph(subdue_graph: Graph):
    reduced_graph = nx.Graph()

    for vertex in subdue_graph.vertices.values():
        vertex_id = int(vertex.id)
        label = vertex.attributes['label']
        reduced_graph.add_node(vertex_id, label=label)

    for edge in subdue_graph.edges.values():
        u = int(edge.source.id)
        v = int(edge.target.id)
        reduced_graph.add_edge(u, v)

    return reduced_graph


def get_pattern_from_subgraph(subgraph, union_graph):
    # get pattern format
    pattern = nx.Graph()
    for e in subgraph['edges']:
        u = union_graph.nodes[e[0]]['label']
        v = union_graph.nodes[e[1]]['label']
        pattern.add_edge(u, v)
    return pattern


def prepare_graphs(protein_graphs: Dict[str, nx.graph]):
    union_graph = nx.Graph()
    pdb_code_to_id = dict()
    graph_offsets = np.empty(len(protein_graphs) + 1, dtype=int)
    graph_offsets[0] = 0
    for pdb_code in protein_graphs:
        g = protein_graphs[pdb_code]
        pdb_code_id = len(pdb_code_to_id) + 1
        pdb_code_to_id[pdb_code] = pdb_code_id
        graph_offsets[pdb_code_id] = graph_offsets[pdb_code_id - 1] + g.number_of_nodes()
        node_ids = {}
        i = 0
        for u in g.nodes:
            uid = str(node_ids.get(u, union_graph.number_of_nodes()))
            union_graph.add_node(uid, label=u)
            node_ids[u] = uid

        for u, v in g.edges:
            union_graph.add_edge(node_ids[u], node_ids[v])

    return union_graph, pdb_code_to_id, graph_offsets


def run_subdue(union_graph, **params):
    out = Subdue.nx_subdue(union_graph, **params)
    return out


class PDBGraphStore:
    def __init__(self, pdb_code_to_id, reduced_graph, patterns, graph_ids_set, graph_offsets):
        self.pdb_code_to_id = pdb_code_to_id
        self.reduced_graph = reduced_graph
        self.patterns = patterns
        self.graph_ids_set = graph_ids_set
        self.graph_offsets = graph_offsets

    def extract_pdb_graph(self, pdb_code):
        graph_id = self.pdb_code_to_id[pdb_code]
        from_vertex = self.graph_offsets[graph_id - 1]
        to_vertex = self.graph_offsets[graph_id]
        sg = nx.subgraph(self.reduced_graph, [u for u in range(from_vertex, to_vertex)])
        nodes = [sg.nodes[u]['label'] for u in sg.nodes]
        edges = [(sg.nodes[u]['label'], sg.nodes[v]['label']) for u, v in sg.edges]
        pdb_subgraph = nx.Graph()
        pdb_subgraph.add_nodes_from(nodes)
        pdb_subgraph.add_edges_from(edges)

        # get patterns
        graphs_to_compose = []
        for i in range(len(self.patterns)):
            if graph_id in self.graph_ids_set[i]:
                graphs_to_compose.append(self.patterns[i])

        graphs_to_compose.append(pdb_subgraph)

        return nx.compose_all(graphs_to_compose)


if __name__ == "__main__":
    import os
    import compress
    from graphein.protein.config import ProteinGraphConfig
    from graphein.protein.edges.atomic import add_atomic_edges
    from graphein.protein.edges.distance import add_hydrogen_bond_interactions, add_peptide_bonds
    from graphein.protein.graphs import construct_graph
    from graphein.protein.utils import download_pdb
    import networkx as nx
    import pickle

    # from subdue.Subdue import nx_subdue

    file_path = os.path.dirname(os.path.realpath(compress.__file__))
    print(file_path)

    params_to_change = {"granularity": "atom", "edge_construction_functions": [add_atomic_edges]}
    #params_to_change = {"granularity": "CA",
    #                    "edge_construction_functions": [add_peptide_bonds, add_hydrogen_bond_interactions]}

    config = ProteinGraphConfig(**params_to_change)
    print(config.model_dump())

    # List of PDB codes
    # pdb_codes = ["1CRN", "4HHB", "2MNR"]  # Example PDB codes
    pdb_codes = []
    with open(f'{file_path}/../../data/soybean_ppigremlin.txt', 'r') as f:
        for line in f:
            pdb_codes.append(line.strip())

    # Function to construct graphs from PDB codes
    protein_graphs = {}
    protein_graphs_with_data = {}
    i = 0
    for pdb_code in pdb_codes:
        print(i, pdb_code)
        i += 1
        try:
            pdb_file = download_pdb(pdb_code, out_dir=f"{file_path}/../../data/pdb_files")  # Download PDB file
            graph = construct_graph(config=config, path=pdb_file)
            protein_graphs_with_data[pdb_code] = graph  # Store graph
            graph_without_data = nx.Graph()
            graph_without_data.add_nodes_from(graph.nodes)
            graph_without_data.add_edges_from(graph.edges)
            protein_graphs[pdb_code] = graph_without_data
        except:
            continue

    pdb_code_to_id, reduced_graph, patterns, graph_ids_set, graph_offsets = compress_with_subdue(protein_graphs,
                                                                                                    node_attributes=[
                                                                                                        'label'], \
                                                                                                    numBest=1,
                                                                                                    iterations=100,
                                                                                                    minSize=0,
                                                                                                    maxSize=5, limit=0, \
                                                                                                    prune=False,
                                                                                                    verbose=False,
                                                                                                    overlap="vertex")

    pdb_graph_store = PDBGraphStore(pdb_code_to_id, reduced_graph, patterns, graph_ids_set, graph_offsets)
    import time

    for pdb_code in protein_graphs:
        start_time = time.time()
        g = protein_graphs[pdb_code]
        elapsed_time = time.time() - start_time
        print(f"g {elapsed_time:f} seconds")

        start_time = time.time()
        h = pdb_graph_store.extract_pdb_graph(pdb_code)
        elapsed_time = time.time() - start_time
        print(f"h {elapsed_time:f} seconds")

        assert id(g) != id(h) and nx.utils.graphs_equal(g, h)

    print(reduced_graph)
    print(len(patterns), [str(p) for p in patterns])
    print(len(graph_ids_set), graph_ids_set)

    import sys


    def networkx_space_cost(g):
        edge_mem = sum([sys.getsizeof(e) for e in g.edges])
        node_mem = sum([sys.getsizeof(n) for n in g.nodes])
        return (edge_mem + node_mem) / 1024 / 1024


    def dict_bitmap_space_cost(d):
        return sys.getsizeof(d) / 1024 / 1024


    patterns_cost = sum([networkx_space_cost(p) for p in patterns])
    reduced_graph_cost = networkx_space_cost(reduced_graph)
    graph_ids_set_cost = dict_bitmap_space_cost(graph_ids_set)
    graph_offsets_cost = sys.getsizeof(graph_offsets) / 1024 / 1024
    compressed_storage = patterns_cost + reduced_graph_cost + graph_ids_set_cost + graph_offsets_cost
    print("CompressedStorage", compressed_storage)

    uncompressed_storage = 0
    for pdb_code in protein_graphs:
        uncompressed_storage += (sys.getsizeof(pdb_code) / 1024 / 1024) + networkx_space_cost(protein_graphs[pdb_code])
    print("UncompressedStorage", uncompressed_storage)

    print("ReductionPercentage", (uncompressed_storage - compressed_storage) / uncompressed_storage)

    print("SizeSerializedCompressed", len(pickle.dumps(pdb_graph_store)))
    print("SizeSerializedUncompressed", len(pickle.dumps(protein_graphs)))
