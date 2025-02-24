import networkx as nx
import matplotlib.pyplot as plt
import os

#protein tool
from graphein.protein.edges.atomic import *
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graphs_mp
from graphein.protein.graphs import construct_graph


#ppi tool
from graphein.ppi.graphs import compute_ppi_graph
from graphein.ppi.edges import add_string_edges, add_biogrid_edges
from graphein.ppi.config import PPIGraphConfig


#fractal
from pyspark.sql import SparkSession
import pyfractal
from pyfractal import vizutil
from pyfractal.model import FractalContext

import ray

from ACTORS_SETUP import *

ray.init(num_cpus=4, object_store_memory=4*1024*1024*1024)


"""
(
       'AP3B1', 
       'SLC44A2', 
       {
              'kind': 
              {
                     'string'
              }
       }
)

uma aresta ppi tem o formato acima

"""

#print(g.nodes(data=True))

"""
(
       'TMED5', 
       {
              'protein_id': 'TMED5'
       }
)

um node de um ppi graph tem o formato acima

"""