import sys

import scanpy as sc

from src.util_anndata import compute_cluster, compute_neighbors, compute_umap, save_umap
from src.util_qc import filter_doublet


def main():
    n_neighbors = int(sys.argv[1])
    n_pcs = int(sys.argv[2])
    resolution = float(sys.argv[3])
    
    total_anndata = sc.read_h5ad("data/processed/anndata/merged_processed.h5ad")

    total_anndata = filter_doublet(total_anndata)

    total_anndata = compute_neighbors(total_anndata, n_neighbors = n_neighbors, n_pcs = n_pcs)
    total_anndata = compute_cluster(total_anndata, resolution = resolution)
    total_anndata = compute_umap(total_anndata)
    save_umap(total_anndata, n_neighbors = n_neighbors, n_pcs = n_pcs, resolution = resolution)
    
if __name__ == "__main__":
    main()