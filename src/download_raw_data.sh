#!/bin/bash

# Download raw data from the web
## Study: An immune-cell signature of bacterial sepsis (Patient PBMCs)
## https://singlecell.broadinstitute.org/single_cell/study/SCP548/an-immune-cell-signature-of-bacterial-sepsis-patient-pbmcs#/
## using bulk download - need the authentication code to download. 
## the authentication code is valid for 30 minutes.
## When click the buld download, you can check download curl commands in the browser console.

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP548&auth_code=eW7JLBAD&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt

mv SCP548/ data/raw/