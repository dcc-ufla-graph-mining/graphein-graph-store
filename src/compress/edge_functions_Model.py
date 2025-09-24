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

from bidict import bidict

edge_functions_dict = bidict({
    "aromatic": add_aromatic_interactions,
    "aromatic_sulphur": add_aromatic_sulphur_interactions,
    "bb_carbonyl_carbonyl": add_backbone_carbonyl_carbonyl_interactions,
    "cation_pi": add_cation_pi_interactions,
    "delaunay": add_delaunay_triangulation,
    "disulfide": add_disulfide_interactions,
    "fully_connected": add_fully_connected_edges,
    "hbond": add_hydrogen_bond_interactions,
    "hydrophobic": add_hydrophobic_interactions,
    "ionic": add_ionic_interactions,
    "knn": add_k_nn_edges,
    "peptide_bond": add_peptide_bonds,
    # "distance_window": add_distance_window,
    # "pi_stacking": add_pi_stacking_interactions,
    # "t_stacking": add_t_stacking,
    # "salt_bridge": add_salt_bridges, 
})