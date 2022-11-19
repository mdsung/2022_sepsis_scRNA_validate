from pathlib import Path

import scanpy as sc
from anndata import AnnData


def normalize(anndata):
    sc.pp.normalize_total(anndata, target_sum=1e4)
    return anndata


def logarithmize(anndata):
    sc.pp.log1p(anndata)
    return anndata


def hvg(anndata):
    sc.pp.highly_variable_genes(
        anndata, min_mean=0.0125, max_mean=3, min_disp=0.5
    )
    anndata.raw = anndata
    anndata = anndata[:, anndata.var.highly_variable]
    return anndata


def scale(anndata):
    sc.pp.scale(anndata, max_value=10)
    return anndata


def pca(anndata):
    sc.tl.pca(anndata, svd_solver="arpack")
    return anndata


def add_name(anndata, geo_number):
    anndata.obs["batch"] = geo_number
    return anndata


def integrate_harmony(anndata: AnnData) -> AnnData:
    sc.external.pp.harmony_integrate(
        anndata,
        "batch",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        max_iter_harmony=100,
    )
    anndata.obsm["X_pca"] = anndata.obsm["X_pca_harmony"]
    return anndata


def compute_neighbors(anndata: AnnData, n_neighbors=10, n_pcs=50) -> AnnData:
    sc.pp.neighbors(anndata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return anndata


def compute_cluster(anndata: AnnData, resolution: float = 0.5) -> AnnData:
    sc.tl.leiden(anndata, resolution=resolution)
    return anndata


def compute_umap(anndata: AnnData) -> AnnData:
    sc.tl.umap(anndata)
    return anndata


def save_umap(
    anndata: AnnData, n_neighbors=10, n_pcs=50, resolution: float = 0.5
):
    Path("figures/umap").mkdir(parents=True, exist_ok=True)
    leiden_umap = sc.pl.umap(
        anndata,
        color="leiden",
        show=False,
        title=f"n_neighbors = {n_neighbors}, n_pcs = {n_pcs}, resolution = {resolution}",
    )
    fig = leiden_umap.get_figure()
    fig.savefig(f"figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png")


def save_anndata(anndata: AnnData, save_path: Path):
    save_path.parent.mkdir(parents=True, exist_ok=True)
    anndata.write(save_path)
