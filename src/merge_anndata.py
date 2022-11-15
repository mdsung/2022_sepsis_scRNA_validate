import sys
from pathlib import Path

import anndata
import pandas as pd
import scanpy as sc
from anndata import AnnData

from src.util_anndata import (
    compute_cluster,
    compute_neighbors,
    compute_umap,
    filter_by_QC,
    hvg,
    integrate_harmony,
    logarithmize,
    normalize,
    pca,
    save_anndata,
    save_umap,
    scale,
)


def load_qc_data(path:Path):
    return pd.read_csv(path)

def load_anndata_from_h5ad(geo_number: str) -> AnnData:
    anndata = sc.read_h5ad(f"data/processed/anndata/{geo_number}.h5ad")
    anndata.obs_names_make_unique()
    return anndata

def main():
    qc_table_path = Path(sys.argv[1]) if sys.argv[1] is not None else Path("figures/qc/qc.csv")
    output_path = Path(sys.argv[2]) if sys.argv[2] is not None else Path("data/processed/anndata/sepsis.h5ad")
    
    qc = load_qc_data(qc_table_path)
    geo_numbers = qc.geo.to_list()

    anndata_list = []
    for geo_number in geo_numbers:
        adata = load_anndata_from_h5ad(geo_number)
        
        n_genes_by_counts = qc.n_genes_by_counts[qc.geo==geo_number].item()
        total_counts = qc.total_counts[qc.geo==geo_number].item()
        pct_count_mt = qc.pct_count_mt[qc.geo==geo_number].item()
        
        adata = filter_by_QC(adata, n_genes_by_counts = n_genes_by_counts, total_counts = total_counts, pct_count_mt = pct_count_mt)

        anndata_list.append(adata)
    
    total_anndata = anndata.concat(anndata_list, index_unique=None)
    total_anndata.obs_names_make_unique()
    total_anndata = normalize(total_anndata)
    total_anndata = logarithmize(total_anndata)
    total_anndata = hvg(total_anndata)
    total_anndata = scale(total_anndata)
    total_anndata = pca(total_anndata)
    total_anndata = integrate_harmony(total_anndata)
    save_anndata(total_anndata, output_path)
    
    # total_anndata = compute_neighbors(total_anndata, n_neighbors = 10, n_pcs = 50)
    # total_anndata = compute_cluster(total_anndata, resolution = 0.5)
    # total_anndata = compute_umap(total_anndata)
    # save_umap(total_anndata)


if __name__ == "__main__":
    main()
