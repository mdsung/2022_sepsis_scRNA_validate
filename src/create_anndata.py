import sys
from pathlib import Path

import scanpy as sc

from src.util_anndata import (
    add_name,
    calculate_QC,
    draw_QC_plot,
    filter_min_genes_cells,
    save_anndata,
)


def load_anndata_from_mtx(geo_number: str):
    adata = sc.read_10x_mtx(
        f"data/raw/{geo_number}/",  # the directory with the `.mtx` file
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=True,
    )
    adata.var_names_make_unique()
    return adata

def process_anndata(adata, geo_number):
    adata = add_name(adata, geo_number)
    return adata

def main():
    geo_number = sys.argv[1]
    output = sys.argv[2]
    if output is None:
        output = f"data/processed/anndata/{geo_number}.h5ad"
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    adata = load_anndata_from_mtx(geo_number)
    adata = process_anndata(adata, geo_number)
    adata = filter_min_genes_cells(adata)
    adata = calculate_QC(adata)
    draw_QC_plot(adata, geo_number)
    save_anndata(adata, Path(f"data/processed/anndata/{geo_number}.h5ad"))

if __name__ == "__main__":
    main()
