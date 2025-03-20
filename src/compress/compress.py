import bisect
import pickle
import random
import time
from typing import Dict

import networkx as nx
import numpy as np
from bidict import bidict
from graphein.protein import add_atomic_edges
from pympler import asizeof
from pyroaring import BitMap, BitMap64

from subdue import Subdue, Graph


def compress_with_composition(protein_graphs):
    edge_to_pdbs = {}
    isolated_node_to_pdbs = {}
    pdb_to_view = {}
    for pdb_code, g in protein_graphs.items():
        pdb_to_view[pdb_code] = ([BitMap64()], [BitMap64()])
        for u in nx.isolates(g):
            pdbs = isolated_node_to_pdbs.get(u, [])
            pdbs.append(pdb_code)
            isolated_node_to_pdbs[u] = pdbs
        for e in g.edges:
            pdbs = edge_to_pdbs.get(e, [])
            pdbs.append(pdb_code)
            edge_to_pdbs[e] = pdbs

    edge_id = 0
    edge_to_id = {}
    for e in edge_to_pdbs:
        edge_to_id[e] = edge_id
        pdbs = edge_to_pdbs[e]
        for pdb_code in pdbs:
            pdb_to_view[pdb_code][1][0].add(edge_id)
        edge_id += 1

    del edge_to_pdbs

    node_id = 0
    isolated_node_to_id = {}
    for u in isolated_node_to_pdbs:
        isolated_node_to_id[u] = node_id
        pdbs = isolated_node_to_pdbs[u]
        for pdb_code in pdbs:
            pdb_to_view[pdb_code][0][0].add(node_id)
        node_id += 1

    del isolated_node_to_pdbs

    isolated_node_to_id = bidict(isolated_node_to_id)
    edge_to_id = bidict(edge_to_id)

    # TODO: assertion code (remove) {
    for pdb_code, view in pdb_to_view.items():
        original_graph = protein_graphs[pdb_code]

        isolated_nodes = [isolated_node_to_id.inverse[node_id] for node_id in view[0][0]]
        edges = [edge_to_id.inverse[edge_id] for edge_id in view[1][0]]
        extracted_graph = nx.Graph()
        extracted_graph.update(edges=edges, nodes=isolated_nodes)

        assert nx.utils.graphs_equal(extracted_graph, original_graph)
    # }

    return PDBGraphStoreBitmap(isolated_node_to_id, edge_to_id, pdb_to_view)


def compress_with_subdue(protein_graphs, **subdue_params):
    union_graph, pdb_code_to_id, graph_offsets = prepare_graphs(protein_graphs)

    # id_to_pdb_code = dict()
    # for pdb_code in pdb_code_to_id:
    #    gid = pdb_code_to_id[pdb_code]
    #    id_to_pdb_code[gid] = pdb_code

    # for u in union_graph.nodes:
    #    print(u, union_graph.nodes[u])

    # print(list(enumerate(graph_offsets)))
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

            # subgraph_patterns = [get_pattern_from_subgraph(s, union_graph) for s in subgraphs]
            # for p1 in subgraph_patterns:
            #    for p2 in subgraph_patterns:
            #        assert nx.utils.graphs_equal(p1, p2), f'{p1} {p2} {pattern}'

            # get a list of graph IDs the pattern occurs
            graph_ids = BitMap()
            for subgraph in subgraphs:
                sample_vertex = int(subgraph['nodes'][0])
                graph_id = bisect.bisect_right(graph_offsets, sample_vertex)
                # print(graph_id, id_to_pdb_code[graph_id], subgraph['nodes'])
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


class PDBGraphStoreSubdue:
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


class PDBGraphStoreBitmap:

    def __init__(self, isolated_node_to_id, edge_to_id, pdb_to_view):
        self.isolated_node_to_id = isolated_node_to_id
        self.edge_to_id = edge_to_id
        self.pdb_to_view = pdb_to_view

    def extract_pdb_graph(self, pdb_code):
        view = self.pdb_to_view.get(pdb_code, None)

        def union_bitmaps(bms):
            union_bm = BitMap64()
            if bms is None: return union_bm

            for bm in bms: union_bm = union_bm | bm
            return union_bm

        isolated_nodes = [self.isolated_node_to_id.inverse[node_id] for node_id in union_bitmaps(view[0])]
        edges = [self.edge_to_id.inverse[edge_id] for edge_id in union_bitmaps(view[1])]
        extracted_graph = nx.Graph()
        extracted_graph.update(edges=edges, nodes=isolated_nodes)
        return extracted_graph

    def run_tree_compression(self):
        # build pairs jaccard
        def get_max_pair(pdb_to_bitmaps):
            best_pair = None
            best_jacc = None
            for pdb_code1, bms1 in pdb_to_bitmaps.items():
                bm_1 = bms1[0]
                for pdb_code2, bms2 in pdb_to_bitmaps.items():
                    if pdb_code1 >= pdb_code2: continue
                    bm_2 = bms2[0]
                    jacc = bm_1.intersection_cardinality(bm_2)
                    if best_jacc is None or jacc > best_jacc:
                        best_jacc = jacc
                        best_pair = (pdb_code1, pdb_code2)
            return best_pair

        def get_optimized_bitmaps(reconstructions):
            while len(reconstructions) >= 2:
                pdb_code1, pdb_code2 = get_max_pair(reconstructions)
                # print(pdb_code1, pdb_code2)
                left = reconstructions[pdb_code1][0]
                right = reconstructions[pdb_code2][0]
                intersection = left & right
                delta_1 = left - intersection
                delta_2 = right - intersection

                left_child = (delta_1, reconstructions[pdb_code1][1])
                right_child = (delta_2, reconstructions[pdb_code2][1])
                reconstructions[pdb_code1 + '_' + pdb_code2] = (
                    intersection, {pdb_code1: left_child, pdb_code2: right_child})

                del reconstructions[pdb_code1]
                del reconstructions[pdb_code2]

            def key_in_keys(key, keys):
                for k in keys:
                    if key in k: return k
                return None

            root_tree_node = list(reconstructions.items())[0]

            pdb_to_view_opt = {}

            for pdb_code in self.pdb_to_view:
                tree_node = root_tree_node
                bms = [tree_node[1][0]] if len(tree_node[1][0]) > 0 else []
                bm_height = 1
                key = key_in_keys(pdb_code, tree_node[1][1])
                while key:
                    bm_height += 1
                    tree_node = (key, tree_node[1][1][key])
                    if len(tree_node[1][0]) > 0:
                        tree_node[1][0].run_optimize()
                        bms.append(tree_node[1][0])
                    key = key_in_keys(pdb_code, tree_node[1][1])
                if len(bms) > 0:
                    pdb_to_view_opt[pdb_code] = bms

            return pdb_to_view_opt

        edge_bitmaps = get_optimized_bitmaps(
            {pdb_code: (view[1][0], {}) for pdb_code, view in self.pdb_to_view.items()})
        isolated_node_bitmaps = get_optimized_bitmaps(
            {pdb_code: (view[0][0], {}) for pdb_code, view in self.pdb_to_view.items()})

        # TODO: assertion code (remove) {
        for pdb_code, view in self.pdb_to_view.items():
            union_bm1 = BitMap64()
            for bm in isolated_node_bitmaps.get(pdb_code, BitMap64()):
                union_bm1 = union_bm1 | bm

            assert union_bm1.symmetric_difference_cardinality(view[0][0]) == 0

            union_bm1 = BitMap64()
            for bm in edge_bitmaps[pdb_code]:
                union_bm1 = union_bm1 | bm

            assert union_bm1.symmetric_difference_cardinality(view[1][0]) == 0
        # }

        # replace existing views with optimized ones
        self.pdb_to_view = {pdb_code: (isolated_node_bitmaps.get(pdb_code, None), edge_bitmaps.get(pdb_code, None)) for
                            pdb_code in self.pdb_to_view}
        self.pdb_to_view = {k: v for k, v in self.pdb_to_view.items() if (k, v) != (None, None)}


if __name__ == "__main__":
    import os
    import compress
    from graphein.protein.config import ProteinGraphConfig
    from graphein.protein.edges.distance import add_hydrogen_bond_interactions, add_peptide_bonds
    from graphein.protein.graphs import construct_graph
    from graphein.protein.utils import download_pdb
    import networkx as nx

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
    # protein_graphs_with_data = {}
    i = 0
    for pdb_code in pdb_codes:
        print(i, pdb_code)
        i += 1
        try:
            pdb_file = download_pdb(pdb_code, out_dir=f"{file_path}/../../data/pdb_files")  # Download PDB file
            graph = construct_graph(config=config, path=pdb_file)
            # protein_graphs_with_data[pdb_code] = graph  # Store graph
            graph_without_data = nx.Graph()
            graph_without_data.add_nodes_from(graph.nodes)
            graph_without_data.add_edges_from(graph.edges)
            protein_graphs[pdb_code] = graph_without_data
        except:
            continue

    print("UncompressedBytes", asizeof.asizeof(protein_graphs) / 1024 / 1024)
    print("UncompressedBytesSer", len(pickle.dumps(protein_graphs)) / 1024 / 1024)
    pdb_graph_store_bitmap = compress_with_composition(protein_graphs)
    print("CompressedBytes", asizeof.asizeof(pdb_graph_store_bitmap) / 1024 / 1024)
    print("CompressedBytesSer", len(pickle.dumps(pdb_graph_store_bitmap)) / 1024 / 1024)
    pdb_graph_store_bitmap.run_tree_compression()
    print("CompressedAndOptimizedBytes", asizeof.asizeof(pdb_graph_store_bitmap) / 1024 / 1024)
    print("CompressedAndOptimizedBytesSer", len(pickle.dumps(pdb_graph_store_bitmap)) / 1024 / 1024)

    # TODO: assertion code (remove) {

    start_time = time.time()
    for pdb_code, protein_graph in protein_graphs.items():
        continue
    elapsed_time = time.time() - start_time
    print("RuntimeIterationUncompressed", elapsed_time)

    start_time = time.time()
    for pdb_code, protein_graph in protein_graphs.items():
        computed_graph = pdb_graph_store_bitmap.extract_pdb_graph(pdb_code)
        assert nx.utils.graphs_equal(computed_graph, protein_graph)
    elapsed_time = time.time() - start_time
    print("RuntimeIterationCompressed", elapsed_time)

    # }
