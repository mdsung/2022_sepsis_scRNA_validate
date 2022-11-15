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
        "data/raw/{geo_number}/barcodes.tsv.gz",
        "data/raw/{geo_number}/features.tsv.gz",
        "data/raw/{geo_number}/matrix.mtx.gz",
        script = "src/create_anndata.py",
    output:
        "data/processed/anndata/{geo_number}.h5ad",
        "figures/qc/{geo_number}.png"
    shell:
        "python {input.script} {wildcards.geo_number} {output[0]}"

rule merge_anndata:
    input:
        expand("data/processed/anndata/{geo_number}.h5ad", geo_number=geo_numbers),
        qc="figures/qc/qc.csv",
        script="src/merge_anndata.py"
    output:
        "data/processed/anndata/sepsis.h5ad"
    shell:
        "python {input.script} {input.qc} {output}"

rule draw_umap:
    input:
        "data/processed/anndata/sepsis.h5ad",
        script="src/draw_umaps.py"
    output:
        "figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png"
    shell:
        "python {input.script} {wildcards.n_neighbors} {wildcards.n_pcs} {wildcards.resolution}"

rule draw_umaps:
    input: 
        expand("figures/umap/umap_{n_neighbors}_{n_pcs}_{resolution}.png", 
            n_neighbors=[5, 10, 15, 20],
            n_pcs=[10, 20, 30, 40, 50],
            resolution=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5])