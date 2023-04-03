## Author: Surbhi Sona ##

library(Seurat)
library(SeuratObject)
library(Matrix)

library(openxlsx)
library(readxl)
library(readr)
library(scales)
library(stringr)
library(data.table)

library(gtools)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(lattice)
library(viridis)

library(Hmisc)
library(dplyr)
library(tidyr)
library(purrr)

# The testing is done per group (total 3) based on BCG treatment status. 


##### Functions ##### 

save_norm_exp_mat1<-function(obj,assay,fname)
{
  raw<-obj@assays[[assay]]@counts %>% as.data.frame.matrix()
  expr<-apply(raw,2,function(x) {x/sum(x) *10000})
  expr<-cbind(rownames(expr),expr)
  rownames(expr)<-NULL
  colnames(expr)<-c('Gene',colnames(expr)[-1])
  expr<-as.data.frame(expr)
  write.table(x = expr,
              file =fname, sep = '\t', 
              row.names = F, col.names = T, quote = F)
}

save_norm_exp_mat<-function(obj,assay,fname)
{
  dat<-obj@assays[[assay]]@data %>% as.data.frame.matrix()
  dat<-cbind(rownames(dat),dat)
  rownames(dat)<-NULL
  colnames(dat)<-c('Gene',colnames(dat)[-1])
  dat<-as.data.frame(dat)
  write.table(x = dat,
              file =fname, sep = '\t', 
              row.names = F, col.names = T, quote = F)
}

save_metadata<-function(obj,meta,fname)
{
  metadat<-as.data.frame(obj@meta.data[,meta])
  metadat<-cbind(rownames(obj@meta.data),metadat)
  rownames(metadat)<-NULL
  colnames(metadat)<-c('Cell',meta)
  
  # Save meta data and table
  
  write.table(x = metadat,
              file =fname, sep = '\t', row.names = F, 
              col.names = T, quote = F)
}

metadata_cutoff_based_subset<-function(obj,meta='cell_label',cutoff=10)
{
  tab<-table(obj@meta.data[,meta]) %>% as.data.frame()
  k<-which(tab$Freq>cutoff)
  keep<-as.character(tab$Var1[k])
  Idents(obj)<-meta
  sub_obj<-subset(obj,idents=keep)
  return(sub_obj)
}


### 1. Read obj and subset by group ###
out<-'/home/sonas/cellphonedb/data/scRNA/'

obj<-readRDS('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/scRNA/DATA/dta_cancer_final.rds')

### 2. Subet by treatment group ###
Idents(obj)<-'Groups'
sub_nw<-subset(obj,ident='Naive_w') # 32909 cells
sub_nwo<-subset(obj,ident='Naive_wo') # 103534 cells
sub_rec<-subset(obj,ident='Recurrent') # 113786 cells

# Combine both naive groups
sub_n<-subset(obj,ident=c('Naive_w','Naive_wo')) # 136443 cells

### 3. Now subset by compartment ### - Do not run

## Naive_w ##
#Idents(sub_nw)<-'low_res_compartment'
#sub_nw_ui<-subset(sub_nw,ident=c('immune','uroepithelial')) # 3379 and 29491 cells respectively

## Naive_wo ##
#Idents(sub_nwo)<-'low_res_compartment'
#sub_nwo_ui<-subset(sub_nwo,ident=c('immune','uroepithelial')) #  20164 and  82945 cells respectively

## Recurrent ##
#Idents(sub_rec)<-'low_res_compartment'
#sub_rec_ui<-subset(sub_rec,ident=c('immune','uroepithelial')) #  31468 and  80402 cells respectively


### 4. check cell count ###

## Naive_w ##
cell_tab<-table(sub_nw@meta.data$cell.names) %>% 
  as.data.frame()

colnames(cell_tab)<-c('cell_type','count')
write.table(cell_tab,paste0(out,'uro_imm_cell_count_Naive_w.txt'), sep='\t',row.names=FALSE)


## Naive_wo ##
cell_tab<-table(sub_nwo@meta.data$cell.names) %>% 
  as.data.frame()

colnames(cell_tab)<-c('cell_type','count')
write.table(cell_tab,paste0(out,'uro_imm_cell_count_Naive_wo.txt'), sep='\t',row.names=FALSE)


## Recurrent ##
cell_tab<-table(sub_rec@meta.data$cell.names) %>% 
  as.data.frame()

colnames(cell_tab)<-c('cell_type','count')
write.table(cell_tab,paste0(out,'uro_imm_cell_count_Recurrent.txt'), sep='\t',row.names=FALSE)


## Naive combined ##
cell_tab<-table(sub_n@meta.data$cell.names) %>% 
  as.data.frame()

colnames(cell_tab)<-c('cell_type','count')
write.table(cell_tab,paste0(out,'uro_imm_cell_count_Naive.txt'), sep='\t',row.names=FALSE)

### 5. Now create cpdb inputs ###

## a) Naive_w ##
prefix<-'Naive_w'

#sub<-metadata_cutoff_based_subset(obj=sub_nw,meta='cell.names',cutoff=10)
sub<-metadata_cutoff_based_subset(obj=sub_nw,meta='cell.names',cutoff=100)

# Save metadata #
save_metadata(sub,meta='cell.names',fname=paste0(out,prefix,'_metadata.txt'))

# Save norm expression matrix #

#save_norm_exp_mat(sub,assay='SCT',fname=paste0(out,prefix,'_counts.txt'))

writeMM(sub@assays$SCT@data, file = paste0(out,prefix,'/matrix.mtx'))

write(x = colnames(sub@assays$SCT@data), file = paste0(out,prefix,"/barcodes.tsv"))

write(x = rownames(sub@assays$SCT@data), file = paste0(out,prefix,"/features.tsv"))



## b) Naive_wo ##

prefix<-'Naive_wo'
sub<-metadata_cutoff_based_subset(obj=sub_nwo,meta='cell.names',cutoff=10) # 

### Save metadata ###
save_metadata(sub,meta='cell.names',fname=paste0(out,prefix,'/',prefix,'_metadata.txt'))

### Save norm expression matrix ###
#save_norm_exp_mat(sub,assay='RNACleaned',fname=paste0(out,prefix,'_counts.txt'))

#save_norm_exp_mat(sub,assay='SCT',fname=paste0(out,prefix,'_counts.txt'))

writeMM(sub@assays$SCT@data, file = paste0(out,prefix,'/matrix.mtx'))

write(x = colnames(sub@assays$SCT@data), file = paste0(out,prefix,"/barcodes.tsv"))

write(x = rownames(sub@assays$SCT@data), file = paste0(out,prefix,"/features.tsv"))

mat<-readMM(paste0(out,prefix,'/matrix.mtx'))

## c) Recurrent ##
prefix<-'Recurrent'

sub<-metadata_cutoff_based_subset(obj=sub_rec,meta='cell.names',cutoff=10) #  cells remain

### Save metadata ###
save_metadata(sub,meta='cell.names',fname=paste0(out,prefix,'/',prefix,'_metadata.txt'))

### Save norm expression matrix ###
#save_norm_exp_mat(sub,assay='RNACleaned',fname=paste0(out,prefix,'_counts.txt'))

writeMM(sub@assays$SCT@data, file = paste0(out,prefix,'/matrix.mtx'))
write(x = colnames(sub@assays$SCT@data), file = paste0(out,prefix,"/barcodes.tsv"))

write(x = rownames(sub@assays$SCT@data), file = paste0(out,prefix,"/features.tsv"))


# d) Naive combined

prefix<-'Naive'

#sub<-metadata_cutoff_based_subset(obj=sub_n,meta='cell.names',cutoff=10) #  cells remain
sub<-metadata_cutoff_based_subset(obj=sub_n,meta='cell.names',cutoff=100)

### Save metadata ###
save_metadata(sub,meta='cell.names',fname=paste0(out,prefix,'/',prefix,'_metadata.txt'))

### Save norm expression matrix ###
#save_norm_exp_mat(sub,assay='RNACleaned',fname=paste0(out,prefix,'_counts.txt'))

writeMM(sub@assays$SCT@data, file = paste0(out,prefix,'/matrix.mtx'))

write(x = colnames(sub@assays$SCT@data), file = paste0(out,prefix,"/barcodes.tsv"))

write(x = rownames(sub@assays$SCT@data), file = paste0(out,prefix,"/features.tsv"))






