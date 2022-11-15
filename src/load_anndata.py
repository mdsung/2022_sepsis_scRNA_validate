import math

import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt


def load_metadata():
    return pd.read_csv('data/processed/metadata.csv')


def load_anndata(geo_number:str):
    adata = sc.read_10x_mtx(
        f"data/raw/{geo_number}/",  # the directory with the `.mtx` file
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=True,
    )
    adata.var_names_make_unique() 
    return adata

def process_anndata(adata, patient_id):
    adata = _filter_min_genes_cells(adata)
    adata = _filter_by_QC(adata)
    adata = _normalize(adata)
    adata = _logarithmize(adata)
    adata = _hvg(adata)
    adata = _scale(adata)
    adata = _pca(adata)
    

    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show = False
        )
        plt.savefig(f'data/processed/qc_violin_{patient_id}.png')
    
    return adata

def _filter_min_genes_cells(adata, min_genes = 200, min_cells = 3 ):
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    return adata
# sc.pl.highest_expr_genes(adata, n_top=20)

def _filter_by_QC(adata, n_genes_by_counts = (200, 6000), total_counts = (1000, math.inf), pct_count_mt = (0, 20)):
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"   
    )
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs.n_genes_by_counts < n_genes_by_counts[1]) & (adata.obs.n_genes_by_counts >n_genes_by_counts[0]), :]
    adata = adata[adata.obs.total_counts > total_counts[0], :]
    adata = adata[adata.obs.pct_counts_mt < pct_count_mt[1], :]
    return adata

def _normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata

def _logarithmize(adata):
    sc.pp.log1p(adata)
    return adata

def _hvg(adata):
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    return adata
    
def _scale(adata):
    sc.pp.scale(adata, max_value=10)
    return adata

def _pca(adata):
    sc.tl.pca(adata, svd_solver='arpack')
    return adata

def main():
    metadata = load_metadata()
    metadata = metadata[metadata.group.isin(['Non-survivor', 'Survivor'])]

    for _, (geo_number, Group, timepoint, group) in metadata.iterrows():
        adata = load_anndata(geo_number)
        adata = process_anndata(adata, Group)
        adata.write(f'data/processed/{geo_number}-{timepoint}-{group}.h5ad')
        
if __name__ == "__main__":
    main()