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


