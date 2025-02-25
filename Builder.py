from ACTORS_SETUP import *

ray.cluster_resources()

@ray.remote
class Builder:
    def __init__(self, id):
        self.id = id

    def build(self, config, kind, pdb_list, storage)-> list:
        print("building nx graph...")
        self.config = config
        self.pdb_list = pdb_list
        self.type = kind
        self.storage = storage

        print(324)

        if self.type == "pdb":
            self._build_pdb()
        # elif self.type == "ppi":
        #     self.graph = self._build_ppi()
        else:
            raise ValueError("O valor do parametro deve ser pdb ou ppi")
    

    def _build_pdb(self):
        if len(self.pdb_list) > 1:
            graph = construct_graphs_mp(pdb_code_it=self.pdb_list, config=self.config, num_cores=6)
        else:
            graph = construct_graph(config=self.config, pdb_code=self.pdb_list[0])
            graph = [graph]

        for i in range(len(graph)):
            ray.get(self.storage.put.remote(self.mounting_ref_pdb_str(self.pdb_list[i]), graph[i]))

    # def _build_ppi(self):
    #     func_list = list()
    #     edge_construction_funcs=self.params["params_to_change"]["edge_construction_functions"]

    #     for func in edge_construction_funcs:
    #         func_list.append(func.__name__)

    #     config = PPIGraphConfig()
    #     graph = compute_ppi_graph(config=config, protein_list=self.pdb_list, edge_construction_funcs=edge_construction_funcs)
    #     self.mounting_ref_str(config=config.dict(), edge_construction_functions=func_list)
        
    #     return ray.put(graph)

    
    
    

    def mounting_ref_pdb_str(self, pdb=None):
        config = self.config.dict()

        config["edge_construction_functions"] = [func.__name__ for func in config["edge_construction_functions"]] if config["edge_construction_functions"] else None
        config["node_metadata_functions"] = [func.__name__ for func in config["node_metadata_functions"]] if config["node_metadata_functions"] else None
        config["edge_metadata_functions"] = [func.__name__ for func in config["edge_metadata_functions"]] if config["edge_metadata_functions"] else None
        config["graph_metadata_functions"] = [func.__name__ for func in config["graph_metadata_functions"]] if config["graph_metadata_functions"] else None
        
        ref_str = f"{config}#{pdb}"
        return ref_str
