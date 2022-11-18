#!/usr/bin/env Rscript

# Author: Surbhi Sona

##### Load dependencies #####

library(dplyr)

source('1_copykat_subset_genomic_coord_per_clade_functions.R')

args<-commandArgs(trailingOnly = TRUE)

mat_path=args[1]
out=args[2]
sname=args[3]
var=args[4]
hclust_path=args[5]
k=as.numeric(args[6])
clus=as.numeric(args[7])

gi_ref_file='/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/copyKat_220kb_genomic_fragments_table.txt'

gi_ref<-read.table(gi_ref_file, header=TRUE,sep='\t')


##### User arguments #####

# Range of cutoff
#cutoff_g<-seq(from=0.015,to=0.055,by=0.005)
#cutoff_l<-cutoff_g*(-1)

#cutoff_g<-seq(from=0.025,to=0.070,by=0.005)
cutoff_g<-0.03
cutoff_l<-cutoff_g*(-1)

# tinglab paths
#out<-'/Volumes/tingalab/Surbhi/PROJECTS_tinglab_drive/scRNA_Projects/BLADDER/copyKat/genes/gain_loss_subset/'

#rec_mat2<-readRDS('/Volumes/tingalab/Surbhi/PROJECTS_tinglab_drive/scRNA_Projects/BLADDER/copyKat/matrix/discovery_subsampled_k14_recurrent_subset_matrix.rds')

# HPC paths
#out<-'/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/genes/gain_loss_subset/'

#rec_mat2<-readRDS('/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/matrix/discovery_subsampled_k14_recurrent_subset_matrix.rds')

mat<-readRDS(mat_path)

hcc<-readRDS(hclust_path)

#Replace . in barcodes to - (if present)
  dot<-grep('\\.',hcc$labels)
  if(length(dot)!=0)
{
     hcc$labels<-gsub('\\.','-',hcc$labels)
   }


hc.umap<-cutree(hcc,k=k)
#clus<-c(1,4,5,6,10) # recurrence-associated
#clus<-c(2,3,7,8,9) # recurrentce-non-associated


idx<-which( hc.umap %in% clus)
bc<-names(hc.umap)[idx]

# Subset cells based on barcode
cell_idx<-which(colnames(mat) %in% bc)
clade_mat<-mat[,cell_idx]
#clade_mat2<-cbind(mat[,1:3],clade_mat)

#nm<-paste0(sname,'_clade',clus[i])
nm<-paste0(sname,'_clade',clus)

if(var=='gain')
{#Gain
cutoff<-cutoff_g
}else if(var=='loss')
  { cutoff<-cutoff_l}else
        {cat('Specify correct variation either gain or loss \n')}



lapply(1:length(cutoff), function(x){subset_by_cutoff(rec_mat=clade_mat,usr_cutoff=cutoff[x],var=var,out=out,sname=nm) })







