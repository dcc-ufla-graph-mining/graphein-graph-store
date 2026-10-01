from __future__ import annotations
from typing import TYPE_CHECKING
import networkx as nx
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from pkg.PDBGraphStore import PDBGraphStore


class Extract:
    """
    Auxiliary class responsible for extracting PDB graphs from a PDBGraphStore.
    Works in coordination with the controlling (mother) PDBGraphStore instance.
    """

    MEILER = {
        "ALA": (1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23),
        "ARG": (2.34, 0.29, 6.13, -1.01, 10.74, 0.36, 0.25),
        "ASN": (1.60, 0.13, 2.95, -0.60, 6.52, 0.21, 0.22),
        "ASP": (1.60, 0.11, 2.78, -0.77, 2.95, 0.25, 0.20),
        "CYS": (1.77, 0.13, 2.43, 1.54, 6.35, 0.17, 0.41),
        "GLN": (1.56, 0.18, 3.95, -0.22, 5.65, 0.35, 0.25),
        "GLU": (1.56, 0.15, 3.78, -0.64, 3.09, 0.42, 0.21),
        "GLY": (0.00, 0.00, 0.00, 0.00, 6.07, 0.13, 0.15),
        "HIS": (2.99, 0.23, 4.66, 0.13, 7.69, 0.27, 0.30),
        "ILE": (4.19, 0.19, 4.00, 1.80, 6.04, 0.30, 0.45),
        "LEU": (2.59, 0.19, 4.00, 1.70, 6.04, 0.39, 0.31),
        "LYS": (1.89, 0.22, 4.77, -0.99, 9.99, 0.32, 0.27),
        "MET": (2.35, 0.22, 4.43, 1.23, 5.71, 0.38, 0.32),
        "PHE": (2.94, 0.29, 5.89, 1.79, 5.67, 0.30, 0.38),
        "PRO": (2.67, 0.00, 2.72, 0.72, 6.80, 0.13, 0.34),
        "SER": (1.31, 0.06, 1.60, -0.04, 5.70, 0.20, 0.28),
        "THR": (3.03, 0.11, 2.60, 0.26, 5.60, 0.21, 0.36),
        "TRP": (3.21, 0.41, 8.08, 2.25, 5.94, 0.32, 0.42),
        "TYR": (2.94, 0.30, 6.47, 0.96, 5.66, 0.25, 0.41),
        "VAL": (3.67, 0.14, 3.00, 1.22, 6.02, 0.27, 0.49),
        "UNK": (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    }

    ONE_HOT = [
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
        'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'
    ]

    def __init__(self, store: PDBGraphStore | None = None):
        self.store = store

    @property
    def body_parts(self) -> dict | None:
        return self.store.get_body_parts() if self.store else None

    @staticmethod
    def edge_label_undirected(edge_label: tuple) -> tuple:
        return tuple(sorted(edge_label))

    def _get_target_store_and_body_parts(self, store: PDBGraphStore | None) -> tuple[PDBGraphStore, dict]:
        target_store = store or self.store
        if target_store is None:
            raise ValueError("A PDBGraphStore instance must be provided.")
        body_parts = target_store.get_body_parts()
        if body_parts is None:
            raise ValueError("PDBGraphStore body_parts dictionary is not initialized.")
        return target_store, body_parts

    # --- Residue Embeddings Submethods ---

    def get_meiler_by_residue(self, residue_name: str) -> pd.Series:
        """
        Retorna o Meiler de um resíduo específico.
        """
        m = pd.Series(
            self.MEILER[residue_name],
            name=residue_name,
            index=[f'dim_{x}' for x in [1, 2, 3, 4, 5, 6, 7]],
            dtype='float64'
        )
        return m

    def get_one_hot_by_residue(self, residue_name: str) -> np.ndarray:
        """
        Retorna o one_hot array de um resíduo específico.
        """
        unk = np.zeros(20, dtype=np.uint8)
        if residue_name.upper() != 'UNK':
            unk[self.ONE_HOT.index(residue_name.upper())] = 1
        return unk

    # --- Node Reconstruction Submethods ---

    def reconstruct_node_global_attrs(self, node_label: str, extracted_graph: nx.Graph, store: PDBGraphStore | None = None):
        target_store, _ = self._get_target_store_and_body_parts(store)
        attr_aux = node_label.split(":")
        chain_id = attr_aux[0]
        residue_name = attr_aux[1]
        residue_number = attr_aux[2]

        granularity = target_store.get_granularity()

        if granularity == 'atom':
            atom_type = attr_aux[3]
        else:
            atom_type = granularity

        element_symbol = atom_type[0]

        extracted_graph.nodes[node_label]['chain_id'] = chain_id
        extracted_graph.nodes[node_label]['residue_name'] = residue_name
        extracted_graph.nodes[node_label]['residue_number'] = residue_number
        extracted_graph.nodes[node_label]['atom_type'] = atom_type
        extracted_graph.nodes[node_label]['element_symbol'] = element_symbol

        config_str = str(target_store.get_config())
        if 'meiler' in config_str:
            meiler = self.get_meiler_by_residue(residue_name)
            extracted_graph.nodes[node_label]['meiler'] = meiler

        if 'amino_acid_one_hot' in config_str:
            amino_acid_one_hot = self.get_one_hot_by_residue(residue_name)
            extracted_graph.nodes[node_label]['amino_acid_one_hot'] = amino_acid_one_hot

    def reconstruct_node_attrs(self, node_id: int, pdb_id: int, extracted_graph: nx.Graph, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        node_attributes = body_parts["node_attr_keyvalue_mapping"][(pdb_id, node_id)]
        node_label = body_parts["node_label_to_node_id"].inverse[node_id]

        coords_id = node_attributes[:3]
        coords_values = [body_parts["node_coords_values"].inverse[x] for x in coords_id]

        print(f'coords_values: {coords_values}')

        coords = np.array(coords_values, dtype=np.float32)
        b_factor = np.float64(body_parts["node_b_factor_values"].inverse[node_attributes[3]])

        print(f'b_factor: {b_factor}')

        extracted_graph.nodes[node_label]['coords'] = coords
        extracted_graph.nodes[node_label]['b_factor'] = b_factor

    def reconstruct_nodes(self, extracted_graph: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        target_store, body_parts = self._get_target_store_and_body_parts(store)
        for node_label in extracted_graph.nodes:
            node_id = body_parts["node_label_to_node_id"][node_label]
            self.reconstruct_node_global_attrs(node_label, extracted_graph, target_store)
            self.reconstruct_node_attrs(node_id, pdb_id, extracted_graph, target_store)

    # --- Edge Reconstruction Submethods ---

    def reconstruct_edge_kinds(self, attributes: list, g: nx.Graph, edge_label: tuple, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        attr_key = "kind"
        kind_list = []

        for i in range(len(attributes)):
            attr_value_id = attributes[i]
            kind_value = body_parts["edge_attr_values"].inverse[attr_value_id]
            kind_list.append(kind_value)

        kinds = set(kind_list)
        g.edges[edge_label][attr_key] = kinds

    def reconstruct_edge_distance(self, attributes: list, g: nx.Graph, edge_label: tuple, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        attr_key = "distance"
        attr_value_id = attributes[0]

        distance_value = body_parts["edge_attr_values"].inverse[attr_value_id]
        g.edges[edge_label][attr_key] = distance_value

    def reconstruct_edges(self, extracted_graph: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        target_store, body_parts = self._get_target_store_and_body_parts(store)
        for edge_label in extracted_graph.edges:
            e = self.edge_label_undirected(edge_label)

            n1_idx = body_parts["node_label_to_node_id"][e[0]]
            n2_idx = body_parts["node_label_to_node_id"][e[1]]

            edge_id = body_parts["edge_label_to_edge_id"][(n1_idx, n2_idx)]
            attributes = body_parts["edge_attr_keyvalue_mapping"][(pdb_id, edge_id)]

            self.reconstruct_edge_kinds(attributes[1:], extracted_graph, edge_label, target_store)
            self.reconstruct_edge_distance(attributes[:1], extracted_graph, edge_label, target_store)

    # --- Main Extract Orchestration ---

    def extract(self, pdb_to_extract: str, store: PDBGraphStore | None = None) -> nx.Graph:
        target_store, body_parts = self._get_target_store_and_body_parts(store)
        pdb = pdb_to_extract.lower()
        extracted_graph = nx.Graph()
        pdb_id = body_parts["pdb_code_to_id"][pdb]

        nodes = [body_parts["node_label_to_node_id"].inverse[node_id] for node_id in body_parts["pdb_id_to_nodes"][pdb_id]]
        edges = []

        for edge_id in body_parts["pdb_id_to_edges"][pdb_id]:
            (n1_idx, n2_idx) = body_parts["edge_label_to_edge_id"].inverse[edge_id]
            n1 = body_parts["node_label_to_node_id"].inverse[n1_idx]
            n2 = body_parts["node_label_to_node_id"].inverse[n2_idx]
            edges.append((n1, n2))

        extracted_graph.add_nodes_from(nodes)
        extracted_graph.add_edges_from(edges)

        self.reconstruct_nodes(extracted_graph, pdb_id, target_store)
        self.reconstruct_edges(extracted_graph, pdb_id, target_store)

        extracted_graph.graph['config'] = target_store.get_config()
        extracted_graph.graph['pdb_code'] = pdb

        return extracted_graph

    def __call__(self, pdb_to_extract: str, store: PDBGraphStore | None = None) -> nx.Graph:
        return self.extract(pdb_to_extract, store)


# Alias
Extractor = Extract
