#!/bin/bash

# Download raw data from the web
## Study: Dynamic changes in human single cell transcriptional signatures during fatal sepsis
## GSE167363

strings=("GSM5102900" 
        "GSM5102901" 
        "GSM5102902"
        "GSM5102903" 
        "GSM5102904" 
        "GSM5102905" 
        "GSM5511351" 
        "GSM5511352" 
        "GSM5511353" 
        "GSM5511354" 
        "GSM5511355"
        "GSM5511356")

wget -e robots=off -U mozilla 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167363&format=file' -O data/raw/GSE167363_RAW.tar
tar xvf data/raw/GSE167363_RAW.tar -C data/raw/
rm data/raw/GSE167363_RAW.tar

for string in "${strings[@]}"; do
    mkdir -p data/raw/${string};
    mv data/raw/${string}*.gz data/raw/${string}/;
    mv data/raw/${string}/*_matrix.mtx.gz data/raw/${string}/matrix.mtx.gz;
    mv data/raw/${string}/*_features.tsv.gz data/raw/${string}/features.tsv.gz;
    mv data/raw/${string}/*_barcodes.tsv.gz data/raw/${string}/barcodes.tsv.gz;
done

# Create Patient metadata file
echo -e "GSM5102900	HC1" >> data/raw/metadata.txt
echo -e "GSM5102901	HC2" >> data/raw/metadata.txt
echo -e "GSM5102902	NS LS T0" >> data/raw/metadata.txt
echo -e "GSM5102903	NS LS T6" >> data/raw/metadata.txt
echo -e "GSM5102904	S1 T0" >> data/raw/metadata.txt
echo -e "GSM5102905	S1 T6" >> data/raw/metadata.txt
echo -e "GSM5511351	NS ES_T0" >> data/raw/metadata.txt
echo -e "GSM5511352	NS ES_T6" >> data/raw/metadata.txt
echo -e "GSM5511353	S2_T0" >> data/raw/metadata.txt
echo -e "GSM5511354	S2_T6" >> data/raw/metadata.txt
echo -e "GSM5511355	S3_T0" >> data/raw/metadata.txt
echo -e "GSM5511356	S3_T6" >> data/raw/metadata.txt

