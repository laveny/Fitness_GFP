#!/bin/bash



## 1. find barcode from illumina rawdata

Rscript GFP_illumina_run.R

## 2. calculate fitness

Rscript GFP_fitness.R 
