from ACTORS_SETUP import *

@ray.remote
class Parser:
    '''
    This is the class that will transform the nx.Graph object
    into fractal input files
    '''
    def __init__(self, id):
        self.id = id

    def parse(self, graph, dir_name, GLOBAL_STORAGE):
        self.dir_name = dir_name
        self._makedir()
        graphdir = self.dir_name
        print("generating fractal input files...")

        print(self.dir_name)        

        print(graph.edges(data=True))

        if graph == None :
            raise ValueError("graph does not exist")
        
        nvertices = graph.number_of_nodes()
        nedges = graph.number_of_edges()

        vmap = {}

        for u in sorted(list(graph.nodes())):   
            vmap[u] = len(vmap)

        remapped_graph = nx.relabel_nodes(graph, vmap)
        
        vlabels = nx.get_node_attributes(remapped_graph, 'protein_id')

        for i in range(nvertices):
            try: vlabels[i]
            except KeyError: vlabels[i] = ""

        
        elabels = [
            data["kind"] if "kind" in data else ""
            for _, _, data in remapped_graph.edges(data=True)
        ]

        elabels = [list(data)[0] for data in elabels]

        '''
        for i in range(nedges): 
            print(elabels[i])
        '''

        for i in range(nedges):
            try:elabels[i]
            except KeyError: elabels[i] = ""


        with open(f"{graphdir}/metadata", "w") as f:
            f.write(f"{nvertices} {nedges}\n")

        with open(f"{graphdir}/adjlists", "w") as f:
            edge_id_map = {}
            for u in range(nvertices):
                adjlist = sorted(list(remapped_graph.neighbors(u)))
                for i in range(len(adjlist)):
                    v = adjlist[i]
                    e = edge_id_map.get((min(u, v), max(u, v)), len(edge_id_map))
                    edge_id_map[(min(u, v), max(u, v))] = e
                    if i > 0: f.write(" ")
                    f.write(f"{v},{e}")
                f.write("\n")

        with open(f"{graphdir}/vlabels", "w") as f:
            for u in range(nvertices):
                vlabel = vlabels[u]
                
                f.write(f"{vlabel}\n")


        with open(f"{graphdir}/elabels", "w") as f:
            for u in range(nedges):
                elabel = elabels[u]

                if elabel == "SINGLE":
                    elabel = 1
                elif elabel == "DOUBLE":
                    elabel = 2
                elif elabel == "covalent":
                    elabel = 3
                else:
                    elabel = 0

                f.write(f"{elabel}\n")

        return 0

    def _makedir(self):
        #name = "/app/" + self.dir_name
        name = "/home/heliohsilva/projects/ray-fractal/graphein-on-ray/" + self.dir_name
        directory = os.path.expanduser(name)  # Expands '~' to the full home path
        os.makedirs(directory, exist_ok=True)

        self.dir_name = directory