import sys
from pathlib import Path

import anndata as ad
import scanpy as sc
from anndata import AnnData

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

def main():
    geo_numbers = sys.argv[1:]
    
    adata_list = []
    for geo_number in geo_numbers:
        adata = load_anndata_from_mtx(geo_number)
        adata = add_name(adata, geo_number)
        adata_list.append(adata)

    adata = ad.concat(adata_list)
    output_path = Path(f"data/processed/anndata/merged.h5ad")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    save_anndata(adata, output_path)

if __name__ == "__main__":
    main()
