#!/bin/bash



## 1. find barcode from illumina rawdata

Rscript illumina_run.R

## 2. calculate fitness

Rscript cal_fitness.R 
