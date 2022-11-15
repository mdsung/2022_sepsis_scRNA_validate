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
