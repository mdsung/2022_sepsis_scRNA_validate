#!/bin/bash

# Download raw data from the web
## Study: Dynamic changes in human single cell transcriptional signatures during fatal sepsis
## GSE167363

wget -e robots=off -U mozilla 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167363&format=file' -O data/raw/GSE167363_RAW.tar
tar xvf data/raw/GSE167363_RAW.tar -C data/raw/
rm data/raw/GSE167363_RAW.tar

## Create Patient metadata file
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

