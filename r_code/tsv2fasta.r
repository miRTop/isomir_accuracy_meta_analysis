library(seqinr)
library(tidyverse)
setwd(here::here()) # set working directory in project folder

mirx = rio::import("input/nbt.4183-S3.xlsx", skip=2) %>% 
    janitor::clean_names()
s = lapply(mirx$sequence,function(s) {
    s
    })
seqinr::write.fasta(sequences = s,
                    names = mirx$equimolar_sequence_id,
                    file.out = "data/tewari.fasta")
