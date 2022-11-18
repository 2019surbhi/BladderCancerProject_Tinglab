#!/usr/bin/env Rscript

mylib<-'/home/sonas/R/x86_64-pc-linux-gnu-library/4.1/'

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(devtools)
library(dplyr)
library(dendextend)
library(tidyverse)
library(copykat, lib.loc=mylib)

## User parameters ##
args <- commandArgs(trailing = TRUE)
arr=args[1] # array of sample index: eg. 6,7,35
sname<-args[2] #sample name to attach to heatmap output or say 'merged' or 'recurrent'
cut<-as.numeric(args[3]) # input for cutree (k or number of clusters)
mode<-args[4] # enter 'saved' to start load matrix and hclust obj, saved_hclust to load hclust obj and anything else to compute everything

# Default inputs
sfile<-'/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/sample_table.txt'
meta<-'/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/cluster_metadata.rds'

arrayidx<-as.numeric(unlist(strsplit(arr,split=',')))

outdir<-'/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/input/'
out<-'/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/'

seurat_obj<-'/home/sonas/beegfs/data/v11.pass7p.uroepithelial.clean.v2.rds'

## 1. Read Inputs ## 

if(mode=='saved')
{

cat('Reading matrix  \n')
tumor.mat<-readRDS(paste0(out,'matrix/',sname,'_matrix.rds'))

cat('Reading saved hclust obj: ', paste0(out,'hclust/',sname,'_Hclust.rds'), '\n')
  # Read hclust obj
  hcc<-readRDS(paste0(out,'hclust/',sname,'_Hclust.rds'))
  
  #Replace . in barcodes to - (if present)
  dot<-grep('\\.',hcc$labels)
  if(length(dot)!=0)
   {
     hcc$labels<-gsub('\\.','-',hcc$labels)	
   }
}else{
cat('Reading seurat obj \n')

so<-readRDS(seurat_obj)
metadata<-as.data.frame(so@meta.data, check.names=FALSE)

# make dataframe of IDs for array indexing 
nm<-data.frame(unique(so$orig.ident))

# Procure sampleID using arrayindx 
sampleID<-nm[arrayidx,]
rm(so)

# Read copyKat files
cat('Reading CopyKat output files \n')

filename<-paste0(outdir,sampleID,"_CNA.txt")
dat_lst<-lapply(filename,read.delim,sep="\t")
merged_dat<-do.call(cbind,dat_lst)
colnames(merged_dat)<-gsub("\\.", "-", colnames(merged_dat))

# Subset only tumor cells
tumor.cells <- row.names(metadata)
tumor.mat <- merged_dat[, which(colnames(merged_dat) %in% tumor.cells)]

rm(tumor.cells)

# Add chr annotations
tumor.mat<-cbind(merged_dat[,1:3],tumor.mat)

cat('Saving final merged matrix \n')
saveRDS(tumor.mat,paste0(out,'matrix/',sname,'_matrix.rds'))

## 2. Hierarchical Clustering ##
if(mode=='saved_hclust')
{
cat('Reading saved hclust obj: ', paste0(out,'hclust/',sname,'_Hclust.rds'), '\n')
  # Read hclust obj
  hcc<-readRDS(paste0(out,'hclust/',sname,'_Hclust.rds'))

}else{
   cat('Performing hierarchical clustering \n')

   # Perform clustering
   hcc<-hclust(parallelDist::parDist(t(tumor.mat),threads=20, method = "euclidean"), method = "ward.D2")

   #Save
   saveRDS(hcc, paste0(out,'hclust/',sname,'_Hclust', ".rds"))
  }
}

# Cut the hclust tree
cat('cutting tree using specified k \n')

hc.umap<- cutree(hcc,k=cut)
#names(hc.umap)<-gsub('\\.','-',names(hc.umap))

# Save hcluster assignment as table

#Remove non chr annotation entries (if present)
check<-c('chrom','chrompos','abspos')
rem<-match(check,names(hc.umap))
rem<-rem[!is.na(rem)] # Remove NA values

# skip if these strings are not found
if(length(rem>0)) #check if there is at least 1 match
 {
  
  hc.umap<-hc.umap[-rem]
 }

cat('creating df \n')
df<-data.frame(group=hc.umap, cells=names(hc.umap))
rownames(df)<-NULL

s<-strsplit(df$cells,split='_')
sID<-sapply(s,'[',2)

df$sample<-sID

tab<-table(df$group,df$sample) %>% as.data.frame.matrix()

#Read sample info table
cat('Read sample annotation file \n')
stable<-read.delim(sfile,header=TRUE,sep='\t')

# replace sample indx with sampleID
sarray<-colnames(tab)
indx<-match(sarray,stable$sindex)
sample_name<-stable$samples[indx]
colnames(tab)<-sample_name
group<-rownames(tab)
tab<-cbind(group,tab)
rownames(tab)<-NULL

cat('saving proportions table \n')
write.csv(tab,paste0(out,'prop_tables/',sname,'_cut',cut,'_sample_prop_by_hclust_group.csv'),row.names = FALSE)

rm(tab)

# dendrogram
cat('creating dendrogram \n')
dend<-as.dendrogram(hcc)
rm(hcc)

## 3. Colored annotations ##

# hclust group color functions
rbPal4 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2"))
rbPal5 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:1])
rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:5])
rbPal7 <- colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired"))

cols<-c("#ff0000","#ffa500","#008000","#0000ff","#4b0082",
           "#ffff00","#ee82ee","#6d2852","#9cb452","#9a5b7f",
           "#43345b","#244c6c","#e44c77","#5c8ca4","#13f7fb",
           "#11c963","#176046","#ce6155","#edcba7","#66241d")
my_col<-colorRampPalette(cols)

my_col2<-colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))

# Heatmap colors #
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(n = 999)

# chr annotation color #
chr <- as.numeric(tumor.mat$chrom) %% 2+1
rbPal1<-colorRampPalette(c('black','grey'))
CHR<-rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

# color breaks # 
col_breaks = c(seq(-1,-0.4,length=50),seq(-0.4,-0.2,length=150),seq(-0.2,0.2,length=600),seq(0.2,0.4,length=150),seq(0.4, 1,length=50))

# remove chr annotations from matrix
check<-c('chrom','chrompos','abspos')
rem<-match(check,colnames(tumor.mat))
rem<-rem[!is.na(rem)] # Remove NA values

# skip if these strings are not found
if(length(rem>0)) #check if there is at least 1 match
 {
  
  tumor.mat<-tumor.mat[,-(1:3)]
 }


# Now extract sampleID
bc<-colnames(tumor.mat)

# Order cell barcodes based on hclust cell order
cat('order cells by hclust cell order \n')
o<-match(names(hc.umap),bc)

cat(length(o), 'cells matched \n')
ordered_bc<-bc[o]

cat('Adding scallop gene legend \n')
gene_tab<-readRDS(meta)

select<-match(ordered_bc,rownames(gene_tab))
select<-select[complete.cases(select)]
genes_sub<-gene_tab[select,1] %>% as.data.frame()
genes_sub$cells<-rownames(gene_tab)[select]
colnames(genes_sub)<-c('clusterID','cells')

# Now match cell bc
gene_prog_label<-rep(0,ncol(tumor.mat))

indx<-match(genes_sub$cells,ordered_bc)
gene_prog_label[indx]<-genes_sub$clusterID

rm(gene_tab)

#If the ordered barcodes match the cell barcodes if hc.umap
cat('checking barcode order in hc.umap \n')
identical(names(hc.umap),ordered_bc)
df2<-data.frame('gene_program'=gene_prog_label,'cluster'=hc.umap)
tab2<-table(df2$cluster,df2$gene_program) %>% as.data.frame.matrix()

cat('saving gene prog proportions table \n')
write.csv(tab2,paste0(out,'prop_tables/',sname,'_cut',cut,'_gene_prog_cell_prop_by_hclust_group.csv'),row.names = TRUE)

rm(tab2)

cat('Splitting barcodes to get sampleID \n')
bc_lst<-strsplit(ordered_bc,split='_')
sindx<-sapply(bc_lst,FUN = '[',2)
sample_label<-stable$samples[match(sindx,stable$sindex)]
patient_label<-stable$patient[match(sindx,stable$sindex)]
primary_recurrent_label<-stable$primary_recurrent[match(sindx,stable$sindex)]

ns<-unique(sample_label) %>% length()
np<-unique(patient_label) %>% length()
ng<-unique(gene_prog_label) %>% length()

# Create color annotations
gene_prog<-my_col2(ng)[as.numeric(factor(gene_prog_label))]
sample<-my_col(ns)[as.numeric(factor(sample_label))]
patient<-rbPal7(np)[as.numeric(factor(patient_label))]
primary_recurrent<-c('green','red')[as.numeric(factor(primary_recurrent_label))]

d<-length(unique(hc.umap))
hgrp<-rbPal4(d)[as.numeric(factor(hc.umap))]

names(gene_prog)<-gene_prog_label
names(sample)<-sample_label
names(patient)<-patient_label
names(primary_recurrent)<-primary_recurrent_label
names(hgrp)<-hc.umap


## 4. Print Heatmap ##

# Create legends and legends colors
gcol<-gene_prog[!duplicated(gene_prog)]
g<-names(gcol)

scol<-sample[!duplicated(sample)]
s<-names(scol)

pcol<-patient[!duplicated(patient)]
p<-names(pcol)

rcol<-primary_recurrent[!duplicated(primary_recurrent)]
r<-names(rcol)

hcol<-hgrp[!duplicated(hgrp)]
hg<-names(hcol)

# Dendrogram colored branches
#dend<-dendsort(dend,isReverse = TRUE)

# match cutree labels across colored branches of dendrogram and hclust 
clust.cutree <- dendextend:::cutree(dend, k=cut, order_clusters_as_data = FALSE)
idx <- order(names(clust.cutree))
clust.cutree <- clust.cutree[idx]
df.merge <- merge(hc.umap,clust.cutree,by='row.names')
df.merge.sorted <- df.merge[order(df.merge$y),]
labs<-unique(df.merge.sorted$x)
dend_colrd <- color_branches(dend, k = cut, groupLabels = labs)

#dend_colrd <- color_branches(dend, k = cut)

rm(dend)
rm(df.merge.sorted)
rm(df.merge)
rm(idx)
rm(clust.cutree)
rm(hc.umap)

cat('printing heatmap \n')
h<-55
#h<-25
options(bitmapType='cairo')


# Combine annotations and print heatmap

sep<-rep('white',ncol(tumor.mat))

if(length(unique(sample_label))>length(unique(patient_label)))
{
  if(length(unique(primary_recurrent_label))>1)
  {
    annotations<-rbind(hgrp,hgrp,sep,gene_prog,gene_prog,sep,sample,sample,sep,patient,patient,sep,primary_recurrent,primary_recurrent)
    rownames(annotations)<-c('hcluster','','','gene_prog','','','sample','','','patient','','','primary_recurrent','')
    
    jpeg(paste(out,'final_heatmap/',sname,'_cut',cut,"_Cluster_Heatmaptumoronly_legend.jpeg",sep=""), height=20*250, width=1000, res=100) 
    plot.new()    
    legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
    legend("top",inset=c(2,0.5),legend=r,col=rcol,title='primary/recurrent', pch=15, cex=3, bty='n')
    
    legend("topright",legend=p, col=pcol,title='patient',pch=15, cex=3, bty='n')
    legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
    legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')
    dev.off()
    
  }else
  {
    annotations<-rbind(hgrp,hgrp,sep,gene_prog,gene_prog,sep,sample,sample,sep,patient,patient)
    rownames(annotations)<-c('hcluster','','','gene_prog','','','sample','','','patient','')
    
    jpeg(paste(out,'final_heatmap/',sname,'_cut',cut,"_Cluster_Heatmaptumoronly_legend.jpeg",sep=""), height=20*250, width=1000, res=100)
    plot.new()
    legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
    legend("topright",legend=p, col=pcol,title='patient',pch=15, cex=3, bty='n')
    legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
    legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')
    dev.off()
    
  }
}else
{
  if(length(unique(primary_recurrent_label))>1)
  {
    annotations<-rbind(hgrp,hgrp,sep,gene_prog,gene_prog,sep,sample,sample,sep,primary_recurrent,primary_recurrent)
    rownames(annotations)<-c('hcluster','','','gene_prog','','','sample','','','primary_recurrent','')
    
    jpeg(paste(out,'final_heatmap/',sname,'_cut',cut,"_Cluster_Heatmaptumoronly_legend.jpeg",sep=""), height=20*250, width=2000, res=100)  
    plot.new() 
    legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
    legend("top",inset=c(2,0.5),legend=r,col=rcol,title='primary/recurrent', pch=15, cex=3, bty='n')
    legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
    legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')
    dev.off()
    
  }else
  {
    annotations<-rbind(hgrp,hgrp,sep,gene_prog,gene_prog,sep,sample,sample)
    rownames(annotations)<-c('hcluster','','','gene_prog','','','sample','')
    
    jpeg(paste(out,'final_heatmap/',sname,'_cut',cut,"_Cluster_Heatmaptumoronly_legend.jpeg",sep=""), height=20*250, width=2000, res=100)   
    plot.new()
    legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
    legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
    legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')
    dev.off() 
  }
}
   jpeg(paste(out,'final_heatmap/',sname,'_cut',cut,"_Cluster_Heatmaptumoronly.jpeg",sep=""), height=h*250, width=4000, res=100)

    copykat::heatmap.3(t(tumor.mat),dendrogram="row",
                       ColSideColors=chr1,
                       col=my_palette,RowSideColors = annotations,
                       notecol="black",
                       Colv=NA,Rowv = dend_colrd,
                       breaks=col_breaks, key=TRUE,
                       keysize=1, density.info="histogram",denscol='black',trace="none",
                       cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                       symm=F,symkey=F,symbreaks=F, cex=1,cex.main=4, margins=c(10,10))
 
    dev.off()
    

