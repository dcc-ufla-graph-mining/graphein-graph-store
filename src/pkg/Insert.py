from __future__ import annotations
from typing import TYPE_CHECKING
import networkx as nx
from pyroaring import BitMap64

if TYPE_CHECKING:
    from pkg.PDBGraphStore import PDBGraphStore


class Insert:
    """
    Auxiliary class responsible for inserting PDB graphs into a PDBGraphStore.
    Works in coordination with the controlling (mother) PDBGraphStore instance.
    """

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

    # --- Structure Construction Submethods ---

    def construct_node_structure(self, g: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        for node_label in g.nodes:
            if node_label not in body_parts["node_label_to_node_id"]:
                body_parts["node_label_to_node_id"][node_label] = len(body_parts["node_label_to_node_id"])

            node_id = body_parts["node_label_to_node_id"][node_label]
            body_parts["pdb_id_to_nodes"][pdb_id].add(node_id)

    def construct_edge_structure(self, g: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        for e in g.edges:
            edge_label = self.edge_label_undirected(e)

            n1_idx = body_parts["node_label_to_node_id"][edge_label[0]]
            n2_idx = body_parts["node_label_to_node_id"][edge_label[1]]

            e_ = (n1_idx, n2_idx)

            if e_ not in body_parts["edge_label_to_edge_id"]:
                body_parts["edge_label_to_edge_id"][e_] = len(body_parts["edge_label_to_edge_id"])

            edge_id = body_parts["edge_label_to_edge_id"][e_]
            body_parts["pdb_id_to_edges"][pdb_id].add(edge_id)

    def construct_structure_attributes(self, g: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        self.construct_node_structure(g, pdb_id, store)
        self.construct_edge_structure(g, pdb_id, store)

    # --- Node Attributes Processing Submethods ---

    def process_node_attr_values(self, node: dict, store: PDBGraphStore | None = None) -> list:
        _, body_parts = self._get_target_store_and_body_parts(store)
        node_attr_list = []

        # coords is a np.ndarray of float32
        coords = node["coords"]

        for coord in coords:
            if coord not in body_parts["node_coords_values"]:
                body_parts["node_coords_values"][coord] = len(body_parts["node_coords_values"])

            coord_id = body_parts["node_coords_values"][coord]
            node_attr_list.append(coord_id)

        # b_factor is a float64
        b_factor = node['b_factor']

        if b_factor not in body_parts["node_b_factor_values"]:
            body_parts["node_b_factor_values"][b_factor] = len(body_parts["node_b_factor_values"])

        b_factor_id = body_parts["node_b_factor_values"][b_factor]
        node_attr_list.append(b_factor_id)

        if len(node_attr_list) != 4:
            print("ERROR processing node_attr_values in insert")
            raise ValueError("ERROR")

        return node_attr_list

    def process_node_attrs(self, pdb_id: int, node_id: int, node: dict, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        node_attr_list = self.process_node_attr_values(node, store)
        body_parts["node_attr_keyvalue_mapping"][(pdb_id, node_id)] = node_attr_list

    def process_nodes(self, g: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        for n in g.nodes:
            node_id = body_parts["node_label_to_node_id"][n]
            self.process_node_attrs(pdb_id, node_id, g.nodes[n], store)

    # --- Edge Attributes Processing Submethods ---

    def process_edge_distances(self, distance: float, store: PDBGraphStore | None = None) -> list:
        _, body_parts = self._get_target_store_and_body_parts(store)
        distance_keyvalue_mapping_list = []

        attr_key = "distance"
        attr_value = distance

        if attr_value not in body_parts["edge_attr_values"]:
            body_parts["edge_attr_values"][attr_value] = len(body_parts["edge_attr_values"])

        attr_value_id = body_parts["edge_attr_values"][attr_value]
        distance_keyvalue_mapping_list.append(attr_value_id)

        return distance_keyvalue_mapping_list

    def process_edge_kinds(self, kinds: set, store: PDBGraphStore | None = None) -> list:
        _, body_parts = self._get_target_store_and_body_parts(store)
        kind_keyvalue_mapping_list = []

        attr_key = "kind"
        attr_key_id = body_parts["edge_attr_keys"].index(attr_key)

        for kind in kinds:
            attr_value = kind

            if attr_value not in body_parts["edge_attr_values"]:
                body_parts["edge_attr_values"][attr_value] = len(body_parts["edge_attr_values"])

            attr_value_id = body_parts["edge_attr_values"][attr_value]
            kind_keyvalue_mapping_list.append(attr_value_id)

        return kind_keyvalue_mapping_list

    def process_edge_attrs(self, pdb_id: int, edge_id: int, edge: dict, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        edge_attr_keyvalue_mapping = []

        edge_attr_keyvalue_mapping.extend(self.process_edge_distances(edge["distance"], store))
        edge_attr_keyvalue_mapping.extend(self.process_edge_kinds(edge["kind"], store))

        body_parts["edge_attr_keyvalue_mapping"][(pdb_id, edge_id)] = edge_attr_keyvalue_mapping

    def process_edges(self, g: nx.Graph, pdb_id: int, store: PDBGraphStore | None = None):
        _, body_parts = self._get_target_store_and_body_parts(store)
        for e in g.edges:
            edge = self.edge_label_undirected(e)

            n1_idx = body_parts["node_label_to_node_id"][edge[0]]
            n2_idx = body_parts["node_label_to_node_id"][edge[1]]

            edge_id = body_parts["edge_label_to_edge_id"][(n1_idx, n2_idx)]
            self.process_edge_attrs(pdb_id, edge_id, g.edges[e], store)

    # --- Main Insert Orchestration ---

    def insert(self, pdb_to_insert: dict, store: PDBGraphStore | None = None):
        """
        input: dict[str: nx.Graph]
        """
        target_store, body_parts = self._get_target_store_and_body_parts(store)

        if not target_store.get_config():
            k = list(pdb_to_insert.keys())[0]
            config = pdb_to_insert[k].graph['config']
            target_store.set_config(config)

        for pdb_code, pdb_graph in pdb_to_insert.items():
            pdb_code = pdb_code.lower()

            if pdb_code in body_parts["pdb_code_to_id"]:
                print(f'pdb {pdb_code} is already stored')
                continue
            body_parts["pdb_code_to_id"][pdb_code] = len(body_parts["pdb_code_to_id"])

            pdb_id = body_parts["pdb_code_to_id"][pdb_code]
            if pdb_id not in body_parts["pdb_id_to_edges"]:
                body_parts["pdb_id_to_edges"][pdb_id] = BitMap64()
            if pdb_id not in body_parts["pdb_id_to_nodes"]:
                body_parts["pdb_id_to_nodes"][pdb_id] = BitMap64()

            self.construct_structure_attributes(pdb_graph, pdb_id, target_store)
            self.process_nodes(pdb_graph, pdb_id, target_store)
            self.process_edges(pdb_graph, pdb_id, target_store)

    def __call__(self, pdb_to_insert: dict, store: PDBGraphStore | None = None):
        return self.insert(pdb_to_insert, store)


# Alias
Inserter = Insert
