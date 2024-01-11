# %%
from subprocess import run

run("lamin init --storage pbmcs --schema bionty", shell=True)

import anndata as ad
import lamindb as ln
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import scanpy as sc

ln.track(ln.Transform(name="My notebook"))


# %%
gene_bt = lb.Gene.bionty(organism="human")

# %%
