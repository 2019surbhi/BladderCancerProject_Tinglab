### This script creates heatmap by calculating mean copykat CNV values per clade ###

# Author: Surbhi Sona #

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(devtools)
library(dplyr)
library(dendextend)
library(copykat)


### 0. Set vars ###

sname<-'primary_recurrent'
cut<-10

sname<-'samples17_CNA_clean'
cut<-17

out<-'/home/sonas/copykat/'
mode<-'saved'

# Default inputs
sfile<-'/home/sonas/data/sample_table.txt'
meta<-'/home/sonas/data/cluster_metadata.rds'

#arrayidx<-as.numeric(unlist(strsplit(arr,split=',')))
#seurat_obj<-'/home/sonas/beegfs/data/v11.pass7p.uroepithelial.clean.v2.rds'


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
  


##  Hierarchical Clustering ##
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



### 1.5 Process hclust and create dendrogram ### 

# Cut the hclust tree
cat('cutting tree using specified k \n')

hc.umap<- cutree(hcc,k=cut)
#names(hc.umap)<-gsub('\\.','-',names(hc.umap))

#Remove non chr annotation entries (if present)
check<-c('chrom','chrompos','abspos')
rem<-match(check,names(hc.umap))
rem<-rem[!is.na(rem)] # Remove NA values

# skip if these strings are not found
if(length(rem)>0) #check if there is at least 1 match
{
  
  hc.umap<-hc.umap[-rem]
}


# dendrogram
cat('creating dendrogram \n')
#dend<-as.dendrogram(hcc)
#rm(hcc)


###  2. Collapse matrix by clade ###

clades<-unique(hc.umap)

clade_bc<-list()
for(i in 1:length(clades))
{
  select<-which(hc.umap==clades[i])
  clade_bc[[i]]<-names(hc.umap)[select]
}
names(clade_bc)<-paste0('clade',clades)

clade_mean<-list()

# copy chrom annotations
mean_mat<-tumor.mat[,1:3]
cname<-colnames(mean_mat)

for(i in 1:length(clades))
{
  sub_mat<-tumor.mat[,clade_bc[[i]]]
  clade_mean<-apply(sub_mat,1,mean)
  mean_mat<-cbind(mean_mat,clade_mean)
}

colnames(mean_mat)<-c(cname,names(clade_bc))


# remove chr annotations from matrix
rem<-vector()
check<-c('chrom','chrompos','abspos')
for(i in 1:length(check))
{
  r<-which(colnames(mean_mat)==check[i])
  rem<-c(rem,r)
}

rem<-rem[!is.na(rem)] # Remove NA values

# skip if these strings are not found
if(length(rem>0)) #check if there is at least 1 match
{
  
  mean_mat2<-mean_mat[,-rem]
}


# Order cell barcodes based on hclust cell order
cat('order clades by recurrent status \n')


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

my_col2<-colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu"))


# Heatmap colors #
my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(n = 999)

# chr annotation color #
chr <- as.numeric(mean_mat$chrom) %% 2+1
rbPal1<-colorRampPalette(c('black','grey'))
CHR<-rbPal1(2)[as.numeric(chr)]
chr1 <- cbind(CHR,CHR)

# color breaks # - Do not use
# col_breaks = c(seq(-1,-0.4,length=50),
#                seq(-0.4,-0.2,length=150),
#                seq(-0.2,0.2,length=600),
#                seq(0.2,0.4,length=150),
#                seq(0.4, 1,length=50))
# 


### 5b. Print Heatmap ### 

h<-7
w<-6000 #sample17
w<-4000 #sample8
#h<-25
#options(bitmapType='cairo')

#annotations<-rbind(hgrp,hgrp,sep,gene_prog,gene_prog,sep,sample,sample,sep,patient,patient,sep,primary_recurrent,primary_recurrent)
#rownames(annotations)<-c('hcluster','','','gene_prog','','','sample','','','patient','','','primary_recurrent','')

#odr<-c('clade1','clade4','clade5','clade6','clade10','clade2','clade3','clade7','clade8','clade9')

### Grouped by R. and NR ###

# sample17
odr<-c(1,4,6,8,9,10,11,12,13,14,16,17,2,3,5,7,15)

# sample8



odr<-paste0('clade',odr)

mean_mat2<-mean_mat2[,odr]

# Clade anno and legend
clade_anno<-my_col(17)[as.numeric(factor(colnames(mean_mat2)))]
names(clade_anno)<-colnames(mean_mat2)

ccol<-clade_anno[!duplicated(clade_anno)]
cl<-colnames(mean_mat2)


annotations<-rbind(clade_anno,clade_anno)
rownames(annotations)<-c('clades','')

jpeg(paste(out,sname,'_cut',cut,"_mean_mat_Cluster_Heatmap_tumoronly_grouped.jpeg",sep=""), height=h*250, width=w, res=100)

copykat::heatmap.3(t(mean_mat2),dendrogram="none",
                   ColSideColors=chr1,
                   col=my_palette,
                   notecol="black",
                   Colv=NA,Rowv=NA,
                   RowSideColors = annotations,
                   key=TRUE,
                   keysize=1,trace="none",
                   cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                   symm=F,symkey=F,symbreaks=F, cex=1,cex.main=4, margins=c(10,10))
# 
#legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
# legend("top",inset=c(2,0.5),legend=r,col=rcol,title='primary/recurrent', pch=15, cex=3, bty='n')
# 

legend("topleft",inset=c(0.2,0.0001),legend=cl, col=ccol,title='clade',pch=15, cex=3, bty='n',horiz = TRUE)
# legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
# legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')

dev.off()


### Ordered manually ##

# sample 17

odr<-17:1
odr<-paste0('clade',odr)

# sample8
odr<-c(5,6,1,3,4,8,2,7,9,10)
odr<-paste0('clade',odr)


mean_mat2<-mean_mat[,odr]

# Clade anno and legend
clade_anno<-my_col(17)[as.numeric(factor(colnames(mean_mat2)))]
names(clade_anno)<-colnames(mean_mat2)

ccol<-clade_anno[!duplicated(clade_anno)]
cl<-colnames(mean_mat2)

annotations<-rbind(clade_anno,clade_anno)
rownames(annotations)<-c('clades','')


jpeg(paste(out,sname,'_cut',cut,"_mean_mat_Cluster_Heatmaptumoronly_ordered_by_main_heatmap.jpeg",sep=""), height=h*250, width=w, res=100)

copykat::heatmap.3(t(mean_mat2),dendrogram="none",
                   ColSideColors=chr1,
                   col=my_palette,
                   notecol="black",
                   Colv=NA,Rowv=NA,
                   RowSideColors = annotations,
                   key=TRUE,
                   keysize=1,trace="none",
                   cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                   symm=F,symkey=F,symbreaks=F, cex=1,cex.main=4, margins=c(10,10))
# 
#legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
# legend("top",inset=c(2,0.5),legend=r,col=rcol,title='primary/recurrent', pch=15, cex=3, bty='n')
# 

legend("topleft",inset=c(0.2,0.0001),legend=cl, col=ccol,title='clade',pch=15, cex=3, bty='n',horiz = TRUE)


# legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
# legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')

dev.off()

## clustered ##

png(paste(out,sname,'_cut',cut,"_mean_mat_Cluster_Heatmaptumoronly_clustered.jpeg",sep=""), height=h*250, width=w, res=100)

copykat::heatmap.3(t(mean_mat2),dendrogram="row",
                   ColSideColors=chr1,
                   col=my_palette,
                   notecol="black",
                   Colv=NA, Rowv=TRUE,
                   RowSideColors = annotations,
                   key=TRUE,
                   keysize=1,trace="none",
                   cexRow=0.1,cexCol=0.1,cex.main=1,cex.lab=0.1,
                   symm=F,symkey=F,symbreaks=F, cex=1,cex.main=4, margins=c(10,10))
# 
#legend("topleft", legend=s,col=scol,title='sample', pch=15, cex=3, bty='n')
# legend("top",inset=c(2,0.5),legend=r,col=rcol,title='primary/recurrent', pch=15, cex=3, bty='n')
# 

legend("topleft",inset=c(0.2,0.0001),legend=cl, col=ccol,title='clade',pch=15, cex=3, bty='n',horiz = TRUE)
# legend("bottomright",legend=hg,col=hcol,title='h-clusters', pch=15, cex=3, bty='n')
# legend("bottomleft",legend=g,col=gcol,title='gene-program',pch=15, cex=3, bty='n')

dev.off()

