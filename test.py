# %%
from subprocess import run

run("lamin init --storage pbmcs --schema bionty", shell=True)

import anndata as ad
import lamindb as ln
import lnschema_bionty as lb
import numpy as np
import pandas as pd
import scanpy as sc
from lamin_utils import logger

ln.track(transform := ln.Transform(name="My notebook"))


# %%
adata = sc.read_10x_mtx(
    "filtered_gene_bc_matrices/hg19/", var_names="gene_ids", cache=True
)
print("Before save:", adata.var["gene_symbols"].dtype)
ln.Artifact.from_anndata(
    adata, description="raw pbmcs", field=lb.Gene.ensembl_gene_id, organism="human"
).save()
adata.var["gene_symbols"] = adata.var["gene_symbols"].astype(str)
print("After save:", adata.var["gene_symbols"].dtype)

# %%
# Giving up on those unvalidated ensembl ids
genes = lb.Gene.inspect(adata.var.index, lb.Gene.ensembl_gene_id, organism="human")
adata = adata[:, genes.validated]
# %%

# %%


# gene_mapping = {gene.ensembl_gene_id: gene.symbol for gene in genes}
# ln.save(genes)
# adata.var["ensembl_symbol"] = adata.var_names.map(gene_mapping.get)

# %%
sc.pl.highest_expr_genes(
    adata,
    n_top=20,
    gene_symbols="gene_symbols",
)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# %%
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
# %%
artifact = ln.Artifact.from_anndata(
    adata,
    description="filtered pbmcs",
    field=lb.Gene.ensembl_gene_id,
    organism="human",
)


# %%
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# %%
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
artifact = ln.Artifact.from_anndata(
    adata,
    description="normalized_and_logged",
    field=lb.Gene.ensembl_gene_id,
    organism="human",
)
# %%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.leiden(adata)
sc.tl.paga(adata)
sc.tl.umap(adata)
# %%

sc.tl.rank_genes_groups(adata, "leiden", method="t-test")
# %%
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, gene_symbols="gene_symbols")
# %%
marker_genes = [
    "IL7R",
    "CD79A",
    "MS4A1",
    "CD8A",
    "CD8B",
    "LYZ",
    "CD14",
    "LGALS3",
    "S100A8",
    "GNLY",
    "NKG7",
    "KLRB1",
    "FCGR3A",
    "MS4A7",
    "FCER1A",
    "CST3",
    "PPBP",
]

# %%
pd.DataFrame(adata.uns["rank_genes_groups"]["names"]).map(
    adata.var["gene_symbols"].get
).head(10)
# %%
new_cluster_names = [
    "CD4 T",
    "CD14 Monocytes",
    "B",
    "CD8 T",
    "NK",
    "FCGR3A Monocytes",
    "Dendritic",
    "Megakaryocytes",
]

lb.CellType.public().search("helper T cell").head(100)

# %%
# # %%

# # %%
#
# # %%


effector_t_cell = lb.CellType.filter(uid="3nfZTVV4").one()
effector_t_cell

cell_types = lb.CellType.from_values(adata.obs["cell_type"])
ln.save(cell_types)
tissues = lb.Tissue.from_values(adata.obs["tissue"])
ln.save(tissues)
diseases = lb.Disease.from_values(adata.obs["disease"])
ln.save(diseases)

# register artifact and annotate with features & labels
artifact = ln.Artifact.from_anndata(
    adata,
    description="anndata with obs",
    field=lb.Gene.ensembl_gene_id,
    organism=("human"),
)
features = ln.Feature.lookup()
artifact.labels.add(cell_types, features.cell_type)
artifact.labels.add(tissues, features.tissue)
artifact.labels.add(diseases, features.disease)
# %%
from typing import Any, Iterable, NamedTuple, Protocol, TypeVar

FT = TypeVar("FT", bound=tuple[str, ...], covariant=True)


class NamedTupleProtocol(Protocol[FT]):
    _fields: FT


from collections import namedtuple


def foo(name: str, fields: Iterable[str]) -> NamedTupleProtocol:
    return namedtuple(name, fields)(1, 2)


type(foo("a", ("a", "b")))
# %%
sc.pl.highest_expr_genes(
    adata,
)
# %%
# ln.save(genes)
# obs_features = ln.Feature.from_df(adata.obs)
# ln.save(obs_features)

# # validate and register labels
# cell_types = lb.CellType.from_values(adata.obs["cell_type"])
# ln.save(cell_types)
# tissues = lb.Tissue.from_values(adata.obs["tissue"])
# ln.save(tissues)
# diseases = lb.Disease.from_values(adata.obs["disease"])
# ln.save(diseases)

# # register artifact and annotate with features & labels
# artifact = ln.Artifact.from_anndata(
#     adata,
#     description="anndata with obs",
#     field=lb.Gene.ensembl_gene_id,
#     organism=(  # optionally, globally set organism via lb.settings.organism = "human"
#         "human"
#     ),
# )
# artifact.save()
# features = ln.Feature.lookup()
# artifact.labels.add(cell_types, features.cell_type)
# artifact.labels.add(tissues, features.tissue)
# artifact.labels.add(diseases, features.disease)


# %%
