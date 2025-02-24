from ACTORS_SETUP import *
\
@ray.remote
class Executor:
    '''
    This is the class that will execute some fractal algorithm
    '''
    def __init__(self, id):
        self.id = id

        
    def setup(self):
        print("inicializing fractal")
        builder = pyfractal.DefaultSparkBuilder()
        builder = builder.master("local[8]")
        builder = builder.config("spark.driver.memory","2g")
        builder = builder.appName("FractalQuickstartApp")
        spark = builder.getOrCreate()

        self.fc = FractalContext(spark)
        return 0

    def run(self, file, is_labeled=True, algorithm="clique"):
        self.file = os.path.abspath("/home/heliohsilva/projects/ray-fractal/graphein-on-ray/" + file)
        #file = os.path.abspath("/app/" + file)

        if not os.path.exists(self.file):  raise FileNotFoundError(self.file)

        if is_labeled:
            fg = self.fc.vertex_labeled_graph(self.file)
        else:
            fg = self.fc.unlabeled_graph(self.file)

        match algorithm:
            case "clique":
                print(self._clique(fg))
            case "enum":
                self._enumeration_tree_pattern_induced(fg)

            case _: raise ValueError("The algorithm values available for now are: clique, enum")

        return 0

    def _clique(self, fg):
        cliques = fg.cliques(2).collect()

        if cliques:
            subgraphs_titles = {s:f"" for s in cliques}
            fig = vizutil.draw_graphs_in_grid(subgraphs_titles, ncols=10, figsize=(30,7), with_labels=True, labeled=False, node_size=500)

            fig.tight_layout()
            fig.savefig(f"{self.file}/graphein-cliques-unlabeled.pdf", bbox_inches="tight")
        else:
            print("there are no cliques")

        return "success"


    def _enumeration_tree_pattern_induced(self, fg):
        tailed_triangle = nx.from_edgelist([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)])
        subgraphs = fg.pfractoid(tailed_triangle).extend(4).subgraphs_networkx().collect()
        vizutil.draw_vertex_enumeration_tree_from_subgraphs(subgraphs, figsize=(20,5), font_size=5, node_size=50)
        plt.savefig(f"{self.file}/graphein-enumeration-tree-pattern-induced.pdf", bbox_inches="tight")

