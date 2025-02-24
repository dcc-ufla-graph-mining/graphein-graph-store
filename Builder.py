from ACTORS_SETUP import *

ray.cluster_resources()

@ray.remote
class Builder:
    '''
    This is the class that will use the graphein api to 
    fetch the raw pdb file and transform it into a nx.Graph
    object
    '''
    def __init__(self, id):
        '''
        graph_type: the type of graph to be created (ppi or pdb)
        pdb_list: the list of pdb or another notation of protein codes to 
                create the graph from it.
        '''
        self.id = id

    def build(self, params):
        print("building nx graph...")

        self.params = params
        self.pdb_list = self.params["pdb_codes"]
        self.type = self.params["graph_type"]

        if self.type == "pdb":
            if len(self.pdb_list) > 1:
                self.pdb_is_list = True
            else:
                self.pdb_is_list = False
            self.graph = self._build_pdb()
        elif self.type == "ppi":
            self.graph = self._build_ppi()
        else:
            raise ValueError("O valor do parametro deve ser pdb ou ppi")
        
        print(self.graph)
        return self.graph

    def _build_ppi(self):
        func_list = list()
        edge_construction_funcs=self.params["params_to_change"]["edge_construction_functions"]

        for func in edge_construction_funcs:
            func_list.append(func.__name__)

        config = PPIGraphConfig()
        graph = compute_ppi_graph(config=config, protein_list=self.pdb_list, edge_construction_funcs=edge_construction_funcs)
        self.mounting_ref_str(config=config.dict(), edge_construction_functions=func_list)
        
        return ray.put(graph)

    def _build_pdb(self):
        #params_to_change pode ser um objeto serializavel
        params_to_change = self.params["params_to_change"]
        params_to_change["edge_construction_functions"] = self.unpack_function_params(params_to_change["edge_construction_functions"])
        config = ProteinGraphConfig(**params_to_change)

        is_list = False

        if self.pdb_is_list:
            graph = construct_graphs_mp(pdb_code_it=self.pdb_list, config=config, num_cores=6)
        else:
            graph = construct_graph(config=config, pdb_code=self.pdb_list[0])
            graph = [graph]

        self.mounting_ref_str(config=config.dict(), edge_construction_functions=None)

        return [ray.put(g) for g in graph]
    
    def unpack_function_params(self, edge_functions):
        function_definition = {
            "add_atomic_edges": add_atomic_edges,
            "add_bond_order": add_bond_order,
            "add_ring_status": add_ring_status,
            "assign_bond_states_to_dataframe": assign_bond_states_to_dataframe,
            "assign_covalent_radii_to_dataframe": assign_covalent_radii_to_dataframe,
            "identify_bond_type_from_mapping": identify_bond_type_from_mapping,
        }

        functions = list()

        for func in edge_functions:
            functions.append(function_definition[func])

        return functions

    def mounting_ref_str(self, config, edge_construction_functions):

        #modificar para exportar a config pra um dict
        ref_str = list()

        for i in range(len(self.pdb_list)):
            self.pdb_list[i] = self.pdb_list[i].upper()
            ref_str.append(f"{config}#{self.pdb_list[i]}#{self.type}#{edge_construction_functions}")
