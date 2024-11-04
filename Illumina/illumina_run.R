rm(list = ls())

library(parallel)
library(plyr)
library(dplyr)
library(stringr)


PATH = "~/GFP_Illumina/rawdata/"

files <- list.files(pattern = 'fq$', recursive = T, path = PATH)

files <- unique(gsub('_[12].fq','',files))

mclapply(files, function(x){
  
  fq1 = paste0(PATH, x, '_1.fq')
  fq2 = paste0(PATH, x, '_2.fq')
  
  output_dir = '~/GFP_Illumina/GFP.barcode/'
  
  output_base = str_match(x,"/([sS]*\\d+[-2]*)[-_]FKDL")
  output_base = output_base[2]
  print(output_base)
  system(command = paste0("~/GFP_Illumina/GFP_index_illumina_barcode_calling.py ",
                          fq1, ' ', fq2, ' ', output_dir, ' ', output_base))
    
  }, mc.cores = 100)
