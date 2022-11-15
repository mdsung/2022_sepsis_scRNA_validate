#!/bin/bash

# Download raw data from the web
## Study: Dynamic changes in human single cell transcriptional signatures during fatal sepsis
## GSE167363

wget -e robots=off -U mozilla 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167363&format=file' -O data/raw/GSE167363_RAW.tar
tar xvf data/raw/GSE167363_RAW.tar -C data/raw/
rm data/raw/GSE167363_RAW.tar

