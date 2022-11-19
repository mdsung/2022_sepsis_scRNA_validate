import pandas as pd

rule download_raw:
    input:
        script="download_raw_data.sh"
    output: 
        "data/raw/metadata.txt"
    shell:
        "bash {input.script}"

rule preprocess_metadata:
    input:
        metadata="data/raw/metadata.txt",
        script="src/preprocess_metadata.py"
    output:
        "data/processed/metadata.csv"
    shell:
        "python {input.script} {input.metadata}"

def get_geo_numbers()->list:
    metadata = pd.read_csv('data/processed/metadata.csv')
    metadata = metadata[metadata.group.isin(['Non-survivor', 'Survivor'])]
    return metadata.GEO.tolist()

geo_numbers = get_geo_numbers()

rule create_anndata:
    input: 
        expand("data/raw/{geo_number}/barcodes.tsv.gz", geo_number = geo_numbers),
        expand("data/raw/{geo_number}/features.tsv.gz", geo_number = geo_numbers),
        expand("data/raw/{geo_number}/matrix.mtx.gz", geo_number = geo_numbers),
        script = "src/create_anndata.py",
    output:
        "data/processed/anndata/merged.h5ad",
    shell:
        "python {input.script} {geo_numbers}"

rule merge_anndata:
    input:
        "data/processed/anndata/merged.h5ad",
        qc="figures/qc/qc.csv",
        script="src/process_anndata.py"
    output:
        "data/processed/anndata/merged_processed.h5ad"
    shell:
        "python {input.script} {input.qc} {input[0]} {output}"

rule draw_umap:
    input:
        "data/processed/anndata/merged_processed.h5ad",
        script="src/draw_umaps.py"
    output:
        "figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png"
    shell:
        "python {input.script} {wildcards.n_neighbors} {wildcards.n_pcs} {wildcards.resolution}"

rule draw_umaps:
    input: 
        expand("figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png", 
            n_neighbors=[20, 25, 30, 35, 40],
            n_pcs=[30, 40, 50],
            resolution=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5])