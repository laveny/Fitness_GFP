library(plyr)
library(dplyr)
library(parallel)
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)
library(ggpointdensity)
library(viridis)
library(ggpubr)
library(gridExtra)
library(rstatix)
library(purrr)
library(cowplot)
library(readxl)

##

rm(list = ls())

##

SaveBase = '~/GFP_Illumina/Analysis/'

if(! dir.exists(SaveBase)){
  
  dir.create(SaveBase,recursive = T)
}


### Files
## PacBio Analysis result for AG and TG mutation library

File.list = list()
File.list$AG = c("AG.havesim.txt", "AG.PASS_pro.txt")  
File.list$TG = c("TG.havesim.txt", "TG.PASS_pro.txt")



## 1. Load in Raw data

setwd("~/GFP_Illumina/GFP.barcode/")

files = list.files(pattern = "\\d.txt")
files = grep(pattern = '-',invert = T, value = T, files)
length(files)
files

results <- mclapply(files, function(x){
  
  Name = sub('.txt','',x)
  fn <- data.table::fread(x, header = F)
  colnames(fn) <- c('barcode',Name)
  return(fn)
}, mc.cores = 56)

fn.raw <- Reduce(function(x,y) merge(x, y, by = 'barcode', all = TRUE), results)

#### SaveOut ####
save(fn.raw, file = paste0(SaveBase,"fn.raw.Rdata"))


## 2. PrePare OD data.frame

## it is a data.frame with each OD represents how many cell densities

od.f <- read_excel('OD_cell_count.xlsx')  

colnames(od.f) <- c('od','cells')

head(od.f)

od.f <- od.f[od.f$od <= 1.58,]

# Using Poly = 10
fit.model <- lm(cells ~ poly(od, 10),data = od.f)

## it is the table of OD during competitive assay

OD.exp.f <- read_excel("/culture_OD.xlsx", sheet = 2)

## generation for repeat1 
OD.exp.f <- OD.exp.f %>% mutate(od = REP1)

OD.exp.f$cells_cons <- predict(fit.model, OD.exp.f)

OD.exp.f$cells = OD.exp.f$cells_cons *  OD.exp.f$VOL

OD.exp.f$REP1_generation = log2(OD.exp.f$cells / OD.exp.f$TRANS)

OD.exp.f <- OD.exp.f %>% group_by(Sample) %>%
    arrange(Stimes) %>%
    mutate(SUM_REP1_generation = cumsum(REP1_generation)) %>%
    ungroup()


## generation for repeat2
OD.exp.f <- OD.exp.f %>% mutate(od = REP2)

OD.exp.f$cells_cons <- predict(fit.model, OD.exp.f)

OD.exp.f$cells = OD.exp.f$cells_cons *  OD.exp.f$VOL

OD.exp.f$REP2_generation = log2(OD.exp.f$cells / OD.exp.f$TRANS)

OD.exp.f <- OD.exp.f %>% group_by(Sample) %>%
    arrange(Stimes) %>%
    mutate(SUM_REP2_generation = cumsum(REP2_generation)) %>%
    ungroup()


## generation for repeat3
OD.exp.f <- OD.exp.f %>% mutate(od = REP3)

OD.exp.f$cells_cons <- predict(fit.model, OD.exp.f)

OD.exp.f$cells = OD.exp.f$cells_cons *  OD.exp.f$VOL

OD.exp.f$REP3_generation = log2(OD.exp.f$cells / OD.exp.f$TRANS)

OD.exp.f <- OD.exp.f %>% group_by(Sample) %>%
    arrange(Stimes) %>%
    mutate(SUM_REP3_generation = cumsum(REP3_generation)) %>%
    ungroup()

## generation for mean of three repeats
OD.exp.f <- OD.exp.f %>% mutate(od = (REP1 + REP2 + REP2) /3)

OD.exp.f$cells_cons <- predict(fit.model, OD.exp.f)

OD.exp.f$cells = OD.exp.f$cells_cons *  OD.exp.f$VOL

OD.exp.f$generation = log2(OD.exp.f$cells / OD.exp.f$TRANS)

OD.exp.f <- OD.exp.f %>% group_by(Sample) %>%
    arrange(Stimes) %>%
    mutate(SUM_generation = cumsum(generation)) %>%
    ungroup()
                 
#### SaveOut ####
save(OD.exp.f, file = paste0(SaveBase,'OD.Rdata'))


## 3. Functions
SUBSET_SAMPLE = function(Sample.fn, TODO_SAM, fn.raw){
        to_contains = unlist(Sample.fn %>% filter(grepl(TODO_SAM, new_Sample)) %>% select(ID))
       pattern <- paste0("^([Ss](", paste(to_contains, collapse = "|"), ")|", paste(to_contains, collapse = "|"), ")$")


        fn <- fn.raw %>% select('barcode',matches(pattern))

        if(ncol(fn) != 14 + 1){

            print('error for selection')
        }

        return(fn)
}

READ_Cluster_file = function(Cluster.file){
  
  cluster.fn <- data.table::fread(Cluster.file, header = F)
  colnames(cluster.fn) <- c('barcode','major_barcode')
  head(cluster.fn)
  return(cluster.fn)
  
}

READ_geno_file = function(geno.file){
  
  geno.fn <- data.table::fread(geno.file, header=F)
  colnames(geno.fn) <- c('bar','geno','Ratio','Pro','stop_codon')
  geno.fn <- geno.fn %>% select(-Ratio)
  head(geno.fn)
  return(geno.fn)
  
}


FUN_Cluster = function(fn, geno.file, Cluster.file, sim = TRUE){
  
  geno.fn = READ_geno_file(geno.file)
  
  cluster.fn = READ_Cluster_file(Cluster.file)
  
  geno.fn <- merge(geno.fn, cluster.fn, by.x = 'bar', by.y = 'major_barcode', all.x = T)
  ## in this file, bar is major_bar, barcode is barcode, it contains bar, geno, Pro, stop_codon and barcode
  
  
  fn[is.na(fn)] = 0
  print('### -------- 1. Selecting Barcodes Tested in PacBio Results -------- ###')
  
  if (sim) {
    fn <- fn %>% filter(barcode %in% geno.fn$barcode)
    print(paste0('\tIncluding similar barcodes : n_barcode is : ',nrow(fn)))
    fn = merge(fn, geno.fn, all.x = T, by.x = 'barcode', by.y = 'barcode')
    fn = fn %>% select(-barcode)%>% rename(c('barcode' = 'bar')) 
    ## sum by major_barcode
    dt <- as.data.table(fn)
    result <- dt[, lapply(.SD, function(x) if (is.numeric(x)) sum(x) else unique(x)[1]), by = barcode]
    fn <- as.data.frame(result)
    print(paste0('\tAfter merging similar bar : n_barcode is : ',nrow(fn)))
    
  } else {
    
    fn <- fn %>% filter(barcode %in% cluster.fn$major_barcode)
    fn = merge(fn, geno.fn, all.x = T, by.x = 'barcode', by.y = 'barcode')
    fn = fn %>% select(-bar)
    print(paste0('\tNOT Including similar barcodes : n_barcode is : ',nrow(fn)))
    
  }
  
  return(fn)
  
}

FUN_READSSUM = function(fn, keywords){
  
  to_sum = grep(keywords, colnames(fn), value = T)
  
  fn <- fn %>% mutate(!! paste0(keywords,'_readsSUM') := rowSums(select(., all_of(to_sum))))
  
  return(fn)
}

FUN_cal_SUM = function(fn, to_sum){
  
  
  for(i in to_sum){
    
    fn <- FUN_READSSUM(fn, i)
  }
  
  return(fn)
  
}

FUN_Drop0Bar = function(fn){
  
  
  print('### -------- 2. Dropping Barcodes Not Appear in T0 Repeats [ T0 readsSUM = 0 ] -------- ###')
  print(paste0("\tBefore Drop Barcode D0 SUM = 0, Number of barcode is : ",nrow(fn)))
  fn = fn %>% filter(D0_readsSUM > 0)
  print(paste0("\tAfter Drop Barcode D0 SUM = 0 is, Number of barcode : ",nrow(fn)))
  
  return(fn)
}


FUN_CUTOFF_cal <- function(list_cutoff = seq(0,2000,20), fn, var){
  
  result.f <- data.frame(CUTOFF = list_cutoff)
  
  result = lapply(list_cutoff, function(x){
    
    return(nrow(fn[fn[[var]] > x,]))
    
  })
  
  result.f$count = unlist(result)
  
  return(result.f)
  
}

FUN_CUTOFF_Draw = function(fn, var, list_cutoff = seq(0,2000,20), list_lines = c(10,50,100,200,500,800,1000), CUT = 100){
  
  MEAN = mean(fn[[var]])
  MEDIAN = median(fn[[var]])
  MAX = max(fn[[var]])
  MIN = min(fn[[var]])
  
  if('barcode' %in% colnames(fn)){ 
    
    Ylabs = 'Number of Barcodes'
    Title = 'WT Minimal Reads for each Barcode'
    
  } else { 
    
    Ylabs = 'Number of Genos'
    Title = 'Variants Minimal Reads for Barcodes SUM'
    
  }
  
  p <-FUN_CUTOFF_cal(list_cutoff, fn, var) %>% 
    ggplot(aes(x = CUTOFF, y = count))+
    geom_point(size = .3)+
    theme_classic()+
    geom_vline(xintercept = list_lines, color = '#339DD4', linetype = 'dotted')+
    geom_vline(xintercept = CUT , color = '#EC7F6B', linetype = 'dotted')+
    labs(x = 'CUTOFF(>x)', y = Ylabs, title = Title,
         subtitle = paste(var, paste(list_lines, collapse = ' '), sep = ' : '), 
         caption = paste0('Max = ', MAX,' ; Median = ', MEDIAN,' ; Mean = ', MEAN,' ; Min = ',MIN))+
    theme(plot.title = element_text(size = 7, hjust = .5),plot.subtitle = element_text(size = 6, hjust = .5),
          axis.title = element_text(size = 5), axis.text = element_text(size = 4),
          plot.caption = element_text(size = 4))
  
  return(p)
  
}

FUN_Plot_Save = function(p, Path, Height, Width, pdf = TRUE){
  
  if (pdf) {ggsave(p, file = paste0(Path,'.pdf'),units = 'cm', height = Height, width = Width) }
  
  jpeg(paste0(Path,'.jpg'), res = 800, units = 'cm', height = Height, width = Width)
  print(p)
  dev.off()
  print(p)
  return(paste0('Saving to ', Path ," : done"))
  
}

FILTER_WT = function(fn, SaveBase, SIM = T){
  
  ## 1. filter D0 readsSUM
  
  ## 1.1 draw cutoff plot
  fn.wt <- fn %>% filter(geno == 'WT')
  p.Wt <- FUN_CUTOFF_Draw(fn.wt,'D0_readsSUM')
  if(SIM){FUN_Plot_Save(p.Wt, Path = paste0(SaveBase,'_CUTOFF_WT_D0_SIM'), 8, 12)} else {FUN_Plot_Save(p.Wt, Path = paste0(SaveBase,'_CUTOFF_WT_D0_UNSIM'), 8, 12)}
  
  ## 1.2 filter D0_readsSUM , remove those < 100
  bar.100 = fn.wt[fn.wt$D0_readsSUM >=100,]$barcode
  
  ## 2. filter D7 readsSUM
  p.Wt <- FUN_CUTOFF_Draw(fn.wt,'D7_readsSUM', list_cutoff = seq(0,100,10), list_lines = c(1,3,5,7, 10, 20, 50), CUT = 5)
  if(SIM){FUN_Plot_Save(p.Wt, Path = paste0(SaveBase,'_CUTOFF_WT_D7_SIM'), 8, 12)} else {FUN_Plot_Save(p.Wt, Path = paste0(SaveBase,'_CUTOFF_WT_D7_UNSIM'), 8, 12)}
  
  bar.5 = fn.wt[fn.wt$D7_readsSUM >=5,]$barcode
  
  
  ## 3. filter D7/D0 ratio  
  fn.wt <- fn.wt %>% mutate(Ratio = D7_readsSUM/D0_readsSUM)
  MEAN = mean(fn.wt$Ratio, na.rm = T)
  SD = sd(fn.wt$Ratio, na.rm = T)
  print(paste0("WT barcodes's MEAN : ", MEAN, " ; SD : ", SD))
  
  fn.wt <- fn.wt %>% mutate(group = ifelse(Ratio > MEAN - SD & Ratio < MEAN + SD, 'PASS', 'FAIL'))
  
  bar.sd = fn.wt[fn.wt$group == 'PASS',]$barcode
  table(fn.wt$group)
  ## 3.1 draw distribution of Ratio
  p <- fn.wt %>% filter(Ratio > 0 & Ratio < 6) %>% ggplot(aes(Ratio)) +
    geom_histogram(aes(fill = group),alpha = 0.5, position = 'identity', bins = 1000) + 
    theme_classic()+
    theme(plot.title = element_text(size = 7, hjust = .5),plot.subtitle = element_text(size = 6, hjust = .5),
          axis.title = element_text(size = 5), axis.text = element_text(size = 4),
          plot.caption = element_text(size = 4))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values = c('FAIL' = '#83868C', "PASS" = '#339DD4'))+
    scale_color_manual(values = c('FAIL' = '#83868C', "PASS" = '#339DD4'))
  
  if(SIM){FUN_Plot_Save(p, Path = paste0(SaveBase,'Ratio_WT_SIM'), 8, 12)} else {FUN_Plot_Save(p, Path = paste0(SaveBase,'Ratio_WT_UNSIM'), 8, 12)}
  
  
  ## 4. intersect
  
  bar.PASS = intersect(bar.100, bar.5)
  bar.PASS = intersect(bar.PASS, bar.sd)
  
  return(list(D0.PASS = bar.100, D7.PASS = bar.5, SD.PASS = bar.sd, FINAL = bar.PASS))
  
}

FILTER_VAR = function(fn, SaveBase, SIM = T){
  
  ## 1. filter D0 readsSUM
  
  ## 1.1 draw cutoff plot
  fn.var <- fn %>% filter(geno != 'WT')
  p.var <- FUN_CUTOFF_Draw(fn.var,'D0_readsSUM')
  if(SIM){FUN_Plot_Save(p.var, Path = paste0(SaveBase,'CUTOFF_var_D0_SIM'), 8, 12)} else {FUN_Plot_Save(p.var, Path = paste0(SaveBase,'CUTOFF_var_D0_UNSIM'), 8, 12)}
  
  ## 1.2 filter D0_readsSUM , remove those < 100
  geno.100 = fn.var[fn.var$D0_readsSUM >=100,]$geno
  
  return(geno.100)
  
  
}

fun_cal_geno = function(fn){
  fn <- fn %>% select(-barcode)
  dt <- as.data.table(fn)
  result <- dt[, lapply(.SD, function(x) if (is.numeric(x)) sum(x) else unique(x)[1]), by = geno]
  fn <- as.data.frame(result)
  return(fn)
}


fun_cal_Ratio = function(fn, TODO_SAM){
    for(x in c('D0')) {
        fn[[paste0(x, '_readsSUM_CPM')]] =  fn[[paste0(x, '_readsSUM')]] / sum(fn[[paste0(x, '_readsSUM')]]) * 10^6
        for(i in c(1,2)){
        fn[[paste0(x, '_rep',i,'_CPM')]] = fn[[paste0(TODO_SAM,'_',x, '_rep',i)]] / sum(fn[[paste0(TODO_SAM,'_',x, '_rep',i)]]) * 10^6     
    }                                        
             
  }
                                                    
   for( x in c('D1','D3','D5','D7')) {
        fn[[paste0(x, '_readsSUM_CPM')]] =  fn[[paste0(x, '_readsSUM')]] / sum(fn[[paste0(x, '_readsSUM')]]) * 10^6
        fn[[paste0(x, '_Ratio')]] =  fn[[paste0(x, '_readsSUM_CPM')]] / fn[["D0_readsSUM_CPM"]]                                      
        for(i in c(1,2,3)){
            fn[[paste0(x, '_rep',i,'_CPM')]] = fn[[paste0(TODO_SAM,'_',x, '_rep',i)]] / sum(fn[[paste0(TODO_SAM,'_',x, '_rep',i)]]) * 10^6
            fn[[paste0(x, '_rep',i,'_Ratio')]] = fn[[paste0(x, '_rep',i,'_CPM')]] / fn[["D0_readsSUM_CPM"]] 
        }                                             
  }                                                 

    fn = as.data.frame(fn)
    return(fn)        
}                                       


fun_cal_fitness = function(fn, OD.exp.f, TODO_SAM){
    
    fn <- fun_cal_Ratio(fn, TODO_SAM)
    
    OD.f <- OD.exp.f %>% filter(Sample == TODO_SAM)

    for( x in c('D1','D3','D5','D7')) {

        generation = as.numeric(OD.f %>% filter(Time == x) %>% select(SUM_generation))
                            
        REF = fn[fn$geno == "WT",][[paste0(x, '_Ratio')]]
        fn[[paste0(x, '_rg')]] = (fn[[paste0(x, '_Ratio')]]/REF) ^ (1/generation)
                                          
        for(i in c(1,2,3)){
            REF = fn[fn$geno == "WT",][[paste0(x, '_rep',i,'_Ratio')]]
            generation = as.numeric(OD.f %>% filter(Time == x) %>% select(paste0('SUM_REP',i,'_generation')))
            fn[[paste0(x, '_rep',i, '_rg')]] = (fn[[paste0(x, '_rep',i,'_Ratio')]]/REF) ^ (1/generation)
                                             
        }                                  

    }

    return(fn)        
}                                                                                                 





## 4. Main Function for fitness calculation


SAM_List = c('AGS','TGS','AGY','TGY','AUU','AUY','AUS','TUU','TUY','TUS')

mclapply(SAM_List, function(TODO_SAM){
        SaveBase = paste0(SaveBase,TODO_SAM,'/')
        if (! dir.exists(SaveBase)){
            
            dir.create(SaveBase)
        }
            ## SUBSET fn.raw
        sample.sub <- Sample.fn %>% filter(grepl(TODO_SAM, official_ID))
        sub.f = SUBSET_SAMPLE(sample.sub, TODO_SAM, fn.raw)
        sub.f = sub.f[sub.f$barcode != "",]
        colnames(sub.f) = gsub('^[Ss]','',colnames(sub.f))
        rename_list <- setNames(sample.sub$official_ID[sample.sub$ID %in% colnames(sub.f)], 
                                sample.sub$ID[sample.sub$ID %in% colnames(sub.f)])
        rename_list[grepl(TODO_SAM, rename_list)]
        sub.f <- sub.f %>% plyr::rename(rename_list)
    
            ## get geno.file and cluster.file
        Cluster.file = File.list[[sub('[SUY]$','',TODO_SAM)]][1]
        geno.file = File.list[[sub('[SUY]$','',TODO_SAM)]][2]


        ## cluster

        fn = FUN_Cluster(sub.f, geno.file, Cluster.file, sim = T)

        stats = list()
        Raw_number_of_reads <- colSums(sub.f[sapply(sub.f, is.numeric)], na.rm = TRUE)
        stats[[1]] = as.data.frame(Raw_number_of_reads)
    
        number_of_reads_PacBio_SIM <- colSums(fn[sapply(fn, is.numeric)], na.rm = TRUE)
        stats[[2]] = as.data.frame(number_of_reads_PacBio_SIM)

    
        ## Cal readsSUM and Drop barcodes with T0 readsSUM == 0
        to_sum = c('D0','D1','D3','D5','D7')
        fn <- FUN_cal_SUM(fn, to_sum)

        fn <- FUN_Drop0Bar(fn)
  
        save(fn, file = paste0(SaveBase, TODO_SAM,'_fn_BeforeFill.Rdata'))
    
        ## filter WT
        ## 1. D0 readsSUM >= 100
        ## 2. D7 readsSUM >= 5
        ## 3. D7 readsSUM / D0 readsSUM within MEAN ± SD
        fn.wt.filter <- FILTER_WT(fn = fn, SaveBase = paste0(SaveBase, TODO_SAM),SIM = T)    
        fn.WT.PASS = fn.wt.filter$FINAL

    
        ## filter WT barcode
        fn <- fn %>% filter(geno != 'WT' | barcode %in% fn.WT.PASS)

        # filter Mutant
        ## D0 readsSUM >= 100
        fn_gene = fun_cal_geno(fn)
    
        fn.VAR.PASS = FILTER_VAR(fn = fn_gene, SaveBase, SIM = T)
    
        fn = fn %>% filter(geno %in% c(fn.VAR.PASS, 'WT'))
    
        ## cal Ratio for per barcode
        fn.bar_Ratio <- fun_cal_Ratio(fn, TODO_SAM)
        save(fn.bar_Ratio, file = paste0(SaveBase, TODO_SAM, '_fn.bar.Ratio.Rdata'))
    
        ## cal fn_gene
        
        fn_gene = fun_cal_geno(fn)
        save(fn_gene,ile = paste0(SaveBase, TODO_SAM, '_fn_gene.Rdata'))
        
    
        number_of_reads_PASS <- colSums(fn[sapply(fn, is.numeric)], na.rm = TRUE)
        number_of_barcode <- sapply(fn, function(x) sum(x > 0))
        stats[[3]] = as.data.frame(number_of_reads_PASS)
        stats[[4]] = as.data.frame(number_of_barcode)
    
    
        # cal fitness
        fn_fitness = fun_cal_fitness(fn_gene, OD.exp.f,TODO_SAM)
        fn_fitness <- fn_fitness %>% mutate(n_mut = ifelse(geno == 'WT', 0, lapply(geno, function(x) length(unlist(strsplit(x, ' ')))) )) %>% mutate(n_mut = as.integer(n_mut))
       
        save(fn_fitness, file = paste0(SaveBase, TODO_SAM, '_fn_fitness.Rdata'))
        
        ## stats
        number_of_Single <- sapply(fn_fitness %>% filter(n_mut == 1), function(x) sum(x > 0))
        number_of_Double <- sapply(fn_fitness %>% filter(n_mut == 2), function(x) sum(x > 0))
        number_of_Triple<- sapply(fn_fitness %>% filter(n_mut == 3), function(x) sum(x > 0))
        stats[[5]] = as.data.frame(number_of_Single)
        stats[[6]] = as.data.frame(number_of_Double)
        stats[[7]] = as.data.frame(number_of_Triple)
       
        save(stats, file = paste0(SaveBase, TODO_SAM, '_stats.Rdata'))

}, mc.cores = length(SAM_List))
