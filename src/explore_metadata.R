library(tidyverse)
library(here)

data <- read_tsv("/home/mindong/2022_sepsis_scRNA_validate/data/raw/SCP548/metadata/scp_meta_updated.txt") %>% filter(NAME !='TYPE')
data 
