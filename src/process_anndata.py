import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from src.util_anndata import (
    hvg,
    integrate_harmony,
    logarithmize,
    normalize,
    pca,
    save_anndata,
    scale,
)
from src.util_qc import (
    calculate_doublet,
    calculate_QC,
    draw_QC_plot,
    evaluate_cell_cycle_genes,
    evaluate_sex_related_genes,
    filter_cell_by_QC,
    filter_doublet,
    filter_gene_by_QC,
    filter_min_genes_cells,
)


def load_qc_data(path:Path):
    return pd.read_csv(path)

def load_anndata_from_h5ad(filepath:Path) -> AnnData:
    anndata = sc.read_h5ad(filepath)
    return anndata

def main():
    qc_table_path = Path("figures/qc/qc.csv")
    input_path = Path("data/processed/anndata/merged.h5ad")
    qc_table_path = Path(sys.argv[1]) if sys.argv[1] is not None else Path("figures/qc/qc.csv")
    input_path = Path(sys.argv[2]) if sys.argv[2] is not None else Path("data/processed/anndata/merged.h5ad")
    output_path = Path(sys.argv[3]) if sys.argv[3] is not None else Path("data/processed/anndata/merged_processed.h5ad")
    
    qc = load_qc_data(qc_table_path)
    
    total_anndata = load_anndata_from_h5ad(input_path)
    total_anndata.obs_names_make_unique()
    total_anndata = calculate_QC(total_anndata)
    draw_QC_plot(total_anndata, "before")
    
    total_anndata = filter_min_genes_cells(total_anndata)

    adata_list = []
    for batch in total_anndata.obs.batch.unique():
        print('Processing batch: ', batch)
        adata = total_anndata[total_anndata.obs.batch == batch]
        
        n_genes_by_counts = qc.n_genes_by_counts[qc.geo==batch].item()
        total_counts = qc.total_counts[qc.geo==batch].item()
        pct_count_mt = qc.pct_count_mt[qc.geo==batch].item()
        pct_count_ribo = 5
                
        adata = filter_cell_by_QC(adata, n_genes_by_counts = n_genes_by_counts, total_counts = total_counts, pct_count_mt = pct_count_mt, pct_count_ribo = pct_count_ribo)
        adata_list.append(adata)
        
    total_anndata = anndata.concat(adata_list, index_unique=None)
    total_anndata = filter_gene_by_QC(total_anndata)
    print('Total number of cells: ', total_anndata.n_obs)
    print('Total number of genes: ', total_anndata.n_vars)
    
    draw_QC_plot(total_anndata, "after")
    
    evaluate_sex_related_genes(total_anndata)
    evaluate_cell_cycle_genes(total_anndata)
    total_anndata = calculate_doublet(total_anndata)
    # total_anndata = filter_doublet(total_anndata)
    
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
