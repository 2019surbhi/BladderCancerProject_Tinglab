#!/usr/bin/env Rscript

# Author: Surbhi Sona

##### Load dependencies #####

library(argparser)
library(dplyr)
library(ggplot2)
library(gtools)
library(readr)
library(readxl)
library(xlsx)


##### User arguments #####

parser<-arg_parser(name="choosing_gi_cutoffs.R",description="workflow for printing genomic interval cutoff")

parser<-add_argument(
  parser,
  arg='--out',
  short = '-o',
  type="character",
  default='',
  help="Enter the output path")

parser<-add_argument(
  parser,
  arg='--prefix',
  short = '-f',
  type="character",
  default='',
  help="Enter file name prefix for output files")

parser<-add_argument(
  parser,
  arg='--threshold_path',
  short = '-t',
  type="character",
  default='',
  help="Enter the path for threshold dir")

parser<-add_argument(
  parser,
  arg='--var',
  short = '-v',
  type="character",
  default='',
  help="Specify whether variation is 'gain' or 'loss' ")

parser<-add_argument(
  parser,
  arg='--chr_size_file',
  short = '-c',
  type="character",
  default='/home/sonas/copykat/misc/hg38.chrom.sizes2.txt',
  help="Enter path and name of chr_size file")




args <- parse_args(parser)



##  Choose final combination of cutoffs ##

#path<-'/Volumes/tingalab/Surbhi/PROJECTS_tinglab_drive/scRNA_Projects/BLADDER/copyKat/genes/gain_loss_subset/primary_recurrent/thresholds/'

args<-list()

args$threshold_path<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/thresholds/'
args$var<-'gain'
args$out<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/thresholds_assessment/'

clades<-1:10
for(i in 1:length(clades))
{
args$prefix<-paste0('primary_recurrent_clade',clades[i])

fargs$prefix<-paste0('primary_recurrent_clade',clades[i])
files<-list.files(args$threshold_path,pattern=args$var,full.names = TRUE)

pattern<-paste0(args$prefix,'_')
f<-files[grep(pattern,basename(files))]

nm<-parse_number(gsub(pattern,'',basename(f))
thresh<-lapply(f,readRDS)
names(thresh)<-nm

df_final<-data.frame()
for(i in 1:length(thresh))
{
 coord<-lapply(thresh[[i]],nrow) %>% unlist()
 df<-data.frame("genomic_segments"=coord,"percentile"=names(coord),"cutoff"=rep(names(thresh)[i], length(coord)))
 
 if(is.null(df_final)==TRUE)
 {
     df_final<-df
 }else
 {
     df_final<-rbind(df_final,df)
 }
}

#label<-seq(0.015,0.055,by=0.005)
#df_final$cutoff<-rep(label,5)

col<-c('#7CAE00','dodgerblue4','#F8766D','red','#00BFC4')


g1<-ggplot(df_final,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=percentile)) +
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',subtitle = paste0('[',args$var,']')) +
  #scale_y_continuous(breaks = round(seq(min(df_final$genomic_segments), max(df_final$genomic_segments), by = 500),1)) +
  #scale_x_continuous(breaks = round(seq(min(df_final$cutoff), max(df_final$cutoff), by = 0.005),digits = 3)) +
  scale_color_manual(values=col)+
  theme_bw()

ggsave(paste0(args$out, args$prefix, 'genomic_interval_cnt_by_cutoff_',args$var,'.png'),width=8,height = 8,units='in',g1)

write.xlsx(x = df_final,file = paste0(args$out,args$prefix, '_',var,'_cell_cnt_by_cutoff.xlsx'))


# Remove 25 p

rem<-grep('25',df_final$percentile)
df_final_sub<-df_final[-rem,]

col2<-c('dodgerblue4','#F8766D','red','#00BFC4')

g2<-ggplot(df_final_sub,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=percentile)) +
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',subtitle = paste0('[',args$var,']')) +
  #scale_y_continuous(breaks = round(seq(min(df_final$genomic_segments), max(df_final$genomic_segments), by = 500),1)) +
  #scale_x_continuous(breaks = round(seq(min(df_final$cutoff), max(df_final$cutoff), by = 0.005),digits = 3)) +
  scale_color_manual(values=col2)+
  theme_bw()

ggsave(paste0(args$out, args$prefix, 'genomic_interval_cnt_by_cutoff_',args$var,'_25p_removed.png'),width=8,height = 8,units='in',g2)
}
