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
        "data/raw/metadata.txt"
    output:
        "data/processed/metadata.csv"
    shell:
        "python src/preprocess_metadata.py"

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
    output:
        "data/processed/{geo_number}.h5ad"
    shell:
        "python src/create_anndata.py {wildcards.geo_number}"


rule merge_anndata:
    input:
        expand("data/processed/{geo_number}.h5ad", geo_number=geo_numbers)