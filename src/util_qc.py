import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr
from anndata import AnnData


def filter_min_genes_cells(
    adata: AnnData, min_genes: int = 200, min_cells=3
) -> AnnData:
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    return adata


def calculate_QC(adata: AnnData) -> AnnData:
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs["percent_mt2"] = (
        np.sum(adata[:, adata.var["mt"]].X, axis=1).A1
        / np.sum(adata.X, axis=1).A1
    )
    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1

    return adata


def draw_QC_plot(adata: AnnData, flag: str):
    Path("figures/qc").mkdir(parents=True, exist_ok=True)
    with plt.rc_context():
        fig = sc.pl.violin(
            adata,
            [
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt",
                "pct_counts_ribo",
                "pct_counts_hb",
            ],
            jitter=0.4,
            multi_panel=True,
            groupby="batch",
            rotation=45,
            show=False,
        )
        plt.savefig(f"figures/qc/QC-feature-{flag}.png")


def filter_cell_by_QC(
    adata: AnnData,
    n_genes_by_counts=6000,
    total_counts=1000,
    pct_count_mt=20,
    pct_count_ribo=5,
) -> AnnData:

    # due to we remove doublet cells with other library
    # adata = _filter_by_total_counts(adata, total_counts)
    # adata = _filter_by_ngenes(adata, n_genes_by_counts)
    adata = _filter_by_mito_genes(adata, pct_count_mt)
    adata = _filter_by_ribo_genes(adata, pct_count_ribo)
    print(f"Remaining cells {adata.n_obs}")
    # sc.pp.filter_genes(adata, min_cells=1)
    return adata


def _filter_by_total_counts(adata: AnnData, total_counts) -> AnnData:
    adata.obs["outlier_total"] = adata.obs.total_counts > total_counts
    print(f"{sum(adata.obs['outlier_total'])} cells with large total counts")
    return adata[~adata.obs["outlier_total"], :]


def _filter_by_ngenes(adata: AnnData, n_genes_by_counts) -> AnnData:
    adata.obs["outlier_ngenes"] = (
        adata.obs.n_genes_by_counts > n_genes_by_counts
    )
    print(
        f"{sum(adata.obs['outlier_ngenes'])} cells with large number of genes"
    )
    return adata[~adata.obs["outlier_ngenes"], :]


def _filter_by_mito_genes(adata: AnnData, pct_count_mt: int) -> AnnData:
    adata.obs["outlier_mt"] = adata.obs.pct_counts_mt > pct_count_mt
    print(
        f"{sum(adata.obs['outlier_mt'])} cells with high % of mitochondrial genes"
    )
    return adata[~adata.obs["outlier_mt"], :]


def _filter_by_ribo_genes(adata: AnnData, pct_count_ribo: int) -> AnnData:
    adata.obs["outlier_ribo"] = adata.obs.pct_counts_ribo < pct_count_ribo
    print(
        f"{sum(adata.obs['outlier_ribo'])} cells with high % of ribosomal genes"
    )
    return adata[~adata.obs["outlier_ribo"], :]


def filter_gene_by_QC(adata):
    malat1 = adata.var_names.str.startswith("MALAT1")

    mito_genes = adata.var_names.str.startswith("MT-")
    hb_genes = adata.var_names.str.contains("^HB[^(P)]")

    remove = np.add(mito_genes, malat1)
    remove = np.add(remove, hb_genes)
    keep = np.invert(remove)

    adata = adata[:, keep]

    print(adata.n_obs, adata.n_vars)
    return adata


def evaluate_sex_related_genes(adata: AnnData):
    annot = sc.queries.biomart_annotations(
        "hsapiens",
        [
            "ensembl_gene_id",
            "external_gene_name",
            "start_position",
            "end_position",
            "chromosome_name",
        ],
    ).set_index("external_gene_name")

    chrY_genes = adata.var_names.intersection(
        annot.index[annot.chromosome_name == "Y"]
    )

    adata.obs["percent_chrY"] = (
        np.sum(adata[:, chrY_genes].X, axis=1).A1
        / np.sum(adata.X, axis=1).A1
        * 100
    )

    adata.obs["XIST-counts"] = adata.X[
        :, adata.var_names.str.match("XIST")
    ].toarray()

    # sc.pl.scatter(adata, x="XIST-counts", y="percent_chrY", color="batch")

    with plt.rc_context():
        sc.pl.violin(
            adata,
            ["XIST-counts", "percent_chrY"],
            jitter=0.4,
            groupby="batch",
            rotation=45,
        )
        plt.savefig(f"figures/qc/sex_related_gene.png")


def evaluate_cell_cycle_genes(adata: AnnData):
    target_url = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
    cell_cycle_genes = [
        line.decode("utf-8").strip()
        for line in urllib.request.urlopen(target_url)
    ]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]

    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    copied_adata = adata.copy()
    # normalize to depth 1000
    sc.pp.normalize_per_cell(copied_adata, counts_per_cell_after=1e4)
    # logaritmize
    sc.pp.log1p(copied_adata)
    # scale
    sc.pp.scale(copied_adata)
    sc.tl.score_genes_cell_cycle(copied_adata, s_genes=s_genes, g2m_genes=g2m_genes)
    with plt.rc_context():
        sc.pl.violin(
            copied_adata,
            ["S_score", "G2M_score"],
            jitter=0.4,
            groupby="batch",
            rotation=45,
            show=False,
        )
        plt.savefig(f"figures/qc/cell_cycle.png")



def calculate_doublet(adata: AnnData):
    scrub = scr.Scrublet(adata.X)
    (
        adata.obs["doublet_scores"],
        adata.obs["predicted_doublets"],
    ) = scrub.scrub_doublets()
    # scrub.plot_histogram()
    sum(adata.obs["predicted_doublets"])
    adata.obs["doublet_info"] = adata.obs["predicted_doublets"].astype(str)
    
    with plt.rc_context():
        sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, groupby = 'doublet_info', rotation=45, show=False)
        plt.savefig(f"figures/qc/doublet.png")
    
    
    # sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # adata = adata[:, adata.var.highly_variable]
    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # sc.pp.scale(adata, max_value=10)
    # sc.tl.pca(adata, svd_solver='arpack')
    # sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    # sc.tl.umap(adata)
    
    # umap_fig = sc.pl.umap(adata, color=['doublet_scores','doublet_info','batch'])
    # umap_fig.savefig(f"figures/qc/doublet.png")
    
    return adata


def filter_doublet(adata: AnnData):
    adata = adata[adata.obs["doublet_info"] == "False", :]
    print(f"Remaining cells {adata.n_obs}")
    return adata
