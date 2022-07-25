#!/usr/bin/env Rscript

# Author: Surbhi Sona

##### Load dependencies #####

library(argparser)
library(dplyr)
library(gtools)

##### User arguments #####

parser<-arg_parser(name="copykat_subset_genomic_coord.R",description="workflow for subsetting gains/loss genomic intervals subset")

parser<-add_argument(
  parser,
  arg='--subset_matrix_path',
  short = '-m',
  type="character",
  default='./',
  help="Enter path to copykat CNA output subset matrix")

parser<-add_argument(
  parser,
  arg='--var',
  short = '-v',
  type="character",
  default='gain',
  help="Specify which coordinates to subset - gain or loss. Default is gain")

parser<-add_argument(
  parser,
  arg='--gi_ref',
  short = '-r',
  type="character",
  default='',
  help="Enter the name of copykat genomic interval reference")


parser<-add_argument(
  parser,
  arg='--chr_size_file',
  short = '-c',
  type="character",
  default='/home/sonas/copykat/misc/hg38.chrom.sizes2.txt',
  help="Enter path and name of chr_size file")

parser<-add_argument(
  parser,
  arg='--prefix',
  short = '-f',
  type="character",
  default='subset',
  help="Enter file prefix for output file. Default is prefix")

parser<-add_argument(
  parser,
  arg='--cell_cutoff',
  short = '-t',
  type="character",
  default='75,50,25',
  help="Enter cell count cutoff separated by , . Default is 75,50,25")

parser<-add_argument(
  parser,
  arg='--out',
  short = '-o',
  type="character",
  default='./',
  help="Enter the path to output dir")



args <- parse_args(parser)



##### Read Inputs #####

## Read copykat CNA subset matrix ##
rec_mat2<-readRDS(args$subset_matrix_path)

## Read gi ref ##
gi_ref<-read.table(args$gi_ref_file, header=TRUE,sep='\t')

##### Functions #####

## This function subsets the cells for each genomic interval that have a copykat CNA value > or < cutoff ##

# mat: copykat output matrix or it's subset
# cutoff: gain or loss cutoff
# var: specify whether copy number variation is 'gain' or 'loss'

cell_count_by_cutoff<-function(mat,cutoff,var)
{

# For each row of input matrix, this function calculates number of cells having values> or < cutoff
  
cell_count<-vector()

if(var=='loss')
{
for(i in 1:nrow(mat))
  {
  len<-length(which(mat[i,]<= (cutoff)))
  cell_count<-append(cell_count,len)
  }
}else if(var=='gain')
{
for(i in 1:nrow(mat))
  {
    len<-length(which(mat[i,]>= (cutoff)))
    cell_count<-append(cell_count,len)
  }
}else
{
  cat('please specify either loss or gain \n')
}

return(cell_count)
}

# Function to apply cell% cutoff

## This function filters genomic intervals based on the given cell% cutoff ##

# cell_count: per genomic interval count for cells that have CNA cutoff > or < than given threshold
# tot_cells: total cells in the given copykat matrix
# cutoff: what cell% cutoff to use to filter genomic intervals

subset_indx<-function(cell_count,tot_cells,cell_pct_cutoff=0.75)
{

idx_list<-list()

for(i in 1:length(pct_cutoff))
{
idx_list[[i]]<-which(cell_count>=(cell_pct_cutoff[i]*tot_cells))
}

names(idx_list)<-paste0(cutoff,'_pct')

return(idx_list)
}

## This function is wrapper function that calls cell_count_by_cutoff() and subset_indx() to create subset of genomic intervals by copykat CNA values cutoff and cell count % cutoff ##

# mat2: copykat output matrix or its subset
# cna_cutoff: CNA cutoff
# cell_pct_cutoff: cell count % cutoff
# var: 'gain' or 'loss'
# out: output directory
# prefix: output file prefix

subset_by_cutoff<-function(mat2,cna_cutoff,cell_pct_cutoff,var,out,prefix)
{
    
# Remove genomic interval columns
mat<-mat2[,-(1:3)]
cols<-ncol(mat)
 
 if(var=='gain')
{
cc<-cell_count_by_cutoff(mat,cutoff=cna_cutoff,var)
g_idx<-subset_indx(cell_count=cc,tot_cells=cols,cell_pct_cutoff=cell_pct_cutoff)

#gain<-lapply(1:length(g_idx),function(x){return(mat2[g_idx[[x]],1:3])})
gain<-lapply(1:length(g_idx), function(x){return(gi_ref[g_idx[[x]],])})

saveRDS(object = gain,paste0(out,prefix,'_cutoff',cna_cutoff,'_gain.rds'))
}else if(var=='loss')
{
 cc<-cell_count_by_cutoff(mat,cutoff=cna_cutoff,var)
 l_idx<-subset_indx(cell_count=cc,tot_cells=cols,cell_pct_cutoff=cell_pct_cutoff)
 
 #loss<-lapply(1:length(l_idx),function(x){return(mat2[l_idx[[x]],1:3])})
 loss<-lapply(1:length(l_idx), function(x){return(gi_ref[g_idx[[x]],])})

 saveRDS(object = loss,paste0(out,prefix,'_cutoff',usr_cutoff,'_loss.rds'))
}else
 {
  cat('you need to specify variation as either gain or loss \n')
 }
}


# Range of cutoff
cutoff_g<-seq(from=0.015,to=0.055,by=0.005)
cutoff_l<-cutoff_g*(-1)

if(arg$var=='gain')
{#Gain
lapply(1:length(cutoff_g), function(x){subset_by_cutoff(rec_mat2,cna_cutoff=cutoff_g[x], cell_pct_cutoff=args$cell_cutoff, var='gain',out=args$out,prefix=args$prefix) })
}else if(args$var=='loss')
{
#Loss
lapply(1:length(cutoff_l), function(x){subset_by_cutoff(rec_mat2,cna_cutoff=cutoff_l[x],cell_pct_cutoff=args$cell_cutoff,var='loss',out=args$out,prefix=args$prefix) })
}else
{cat('Specify correct variation either gain or loss \n')}

