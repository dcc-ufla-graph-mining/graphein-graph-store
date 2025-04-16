from pyspark.sql import SparkSession
import pyfractal
from pyfractal import vizutil
from pyfractal.model import FractalContext

print("inicializing fractal")
builder = pyfractal.DefaultSparkBuilder()
builder = builder.master("local[8]")
builder = builder.config("spark.driver.memory","2g")
builder = builder.appName("FractalQuickstartApp")
spark = builder.getOrCreate()

fc = FractalContext(spark)

fg = fc.unlabeled_graph("test1/")
motif_count = fg.motif_counting(5)
subgraphs_titles = {m:f"count = {c}" for m,c in motif_count}
fig = vizutil.draw_graphs_in_grid(subgraphs_titles, ncols=10, figsize=(10,2), with_labels=False, node_size=700, labeled=False)
fig.tight_layout()
fig.savefig("toygraph-motifs-unlabeled.pdf", bbox_inches="tight")