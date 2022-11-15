
import math
from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData


def filter_min_genes_cells(adata:AnnData, min_genes:int=200, min_cells=3)->AnnData:
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    return adata

def calculate_QC(adata: AnnData)->AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata

def draw_QC_plot(adata:AnnData, geo_number:str):
    Path('figures/qc').mkdir(parents=True, exist_ok=True)
    with plt.rc_context():
        fig = sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show = False
        )
        fig.savefig(f"figures/qc/{geo_number}.png")
    

def filter_by_QC(
    adata:AnnData,
    n_genes_by_counts = 6000,
    total_counts = 1000,
    pct_count_mt = 20,
):
    adata.obs['outlier_mt'] = adata.obs.pct_counts_mt > pct_count_mt
    adata.obs['outlier_total'] = adata.obs.total_counts > total_counts
    adata.obs['outlier_ngenes'] = adata.obs.n_genes_by_counts > n_genes_by_counts

    print('%u cells with high %% of mitochondrial genes' % (sum(adata.obs['outlier_mt'])))
    print('%u cells with large total counts' % (sum(adata.obs['outlier_total'])))
    print('%u cells with large number of genes' % (sum(adata.obs['outlier_ngenes'])))

    adata = adata[~adata.obs['outlier_mt'], :]
    adata = adata[~adata.obs['outlier_total'], :]
    adata = adata[~adata.obs['outlier_ngenes'], :]
    sc.pp.filter_genes(adata, min_cells=1)

    return adata


def normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata


def logarithmize(adata):
    sc.pp.log1p(adata)
    return adata


def hvg(adata):
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5
    )
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    return adata


def scale(adata):
    sc.pp.scale(adata, max_value=10)
    return adata


def pca(adata):
    sc.tl.pca(adata, svd_solver="arpack")
    return adata


def add_name(adata, geo_number):
    anndata.obs["batch"] = geo_number
    return anndata


def integrate_harmony(anndata:AnnData) -> AnnData:
    sc.external.pp.harmony_integrate(
        anndata,
        "batch",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        max_iter_harmony=100,
    )
    anndata.obsm["X_pca"] = anndata.obsm["X_pca_harmony"]
    return anndata

def compute_neighbors(anndata:AnnData, n_neighbors = 10, n_pcs=50)->AnnData:
    sc.pp.neighbors(anndata, n_neighbors=n_neighbors, n_pcs=n_pcs)    
    return anndata

def compute_cluster(anndata:AnnData, resolution:float = 0.5)->AnnData:
        sc.tl.leiden(anndata, resolution=resolution)
        return anndata

def compute_umap(anndata: AnnData)->AnnData:
    sc.tl.umap(anndata)
    return anndata

def save_umap(anndata: AnnData, n_neighbors = 10, n_pcs=50, resolution:float = 0.5):
    Path("figures/umap").mkdir(parents=True, exist_ok=True)
    leiden_umap = sc.pl.umap(anndata, color='leiden', show=False, title = f"n_neighbors = {n_neighbors}, n_pcs = {n_pcs}, resolution = {resolution}")
    fig = leiden_umap.get_figure()
    fig.savefig(f"figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png")

def save_anndata(anndata: AnnData, save_path: Path):
    save_path.parent.mkdir(parents=True, exist_ok=True)
    anndata.write(save_path)

