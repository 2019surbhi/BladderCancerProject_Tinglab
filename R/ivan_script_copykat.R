library(data.table)
library(phylogram)
library(pheatmap)


##### 1. Load newick tree ######

#source('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/newick_tree_copyKat.txt')
#source('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/newick_tree_copyKat_17samples.txt')
#source('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/newick_tree_copyKat_17samples_filtered.txt')

#source('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/newick_tree_copyKat_17samples_filtered.txt')

# save_pheatmap_png <- function(x, filename, width=2048, height=2048) 
# { stopifnot(!missing(x))
#   stopifnot(!missing(filename))
#   png(filename, width=width, height=height)
#   grid::grid.newpage() 
#   grid::grid.draw(x$gtable) 
#   dev.off()
# }

#H_ = 200

sname<-'samples27'
out<-'/home/sonas/copykat/'

sub<-0.8
sub<-0.7
sub<-0.6
sub<-0.5

#d.1<-fread(paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/newick_tree/table.',sub,'.csv'), sep = '\t')

source(paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/newick_tree_per_sample_data/newick_tree_copyKat.',sub,'.txt'))

# Change index to 1 from 0
x = read.dendrogram(text = newick_tree)
xh = as.hclust(x)
#plot(xh, labels = FALSE)
label_index = as.numeric(xh$labels) + 1

### save hclust as rds if you need it. 
#saveRDS(xh,'/home/sonas/copykat/copyKat_dendogram_full_data_unordered.rds')

###### 2. Load data ########

# Order 
#d.1 = fread('/home/sonas/copyKat/all_CNA_clean.txt', sep = '\t') 
#d.1 = fread('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/all_CNA_clean.txt', sep = '\t')
#d.1 = fread('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/samples17_CNA_clean.txt', sep = ',')

#d.1 = fread('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/samples17_CNA_clean_filtered.txt',sep = ',')

#d.1 = fread('/home/jurici/Ting_SingleCell_Bladder/copyKat/cluster/samples17_CNA_clean_filtered.txt',sep = ',')

#d.1<-fread(paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/newick_tree/table.',sub,'.csv'))

dir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/newick_tree_per_sample_data/'

## Select Naive only ##



## Select Receurrent only ##



fname<-paste0(dir,'table.',sub,'.S1.csv')
tab<-fread(fname)
d.1<-tab[,-1]

for(i in 2:27)
{
 fname<-paste0(dir,'table.',sub,'.S',i,'.csv')
 tab<-fread(fname)
 tab<-tab[,5:ncol(tab)]
 d.1<-cbind(d.1,tab)
}


# png(paste0('/home/sonas/copykat/pheatmap_h.png'))
# PH.1.python<-pheatmap(d.1r, 
#                       cluster_rows = TRUE, 
#                       show_colnames = FALSE,
#                       color = colorRampPalette(c("blue", "white", "red"))(50), 
#                       scale = 'row', border_color = NA)
# 
# dev.off()

#####

# correct barcode

#d.1<-fread(paste0('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/newick_tree/sub',sub,'_bc.csv'), header=TRUE)

#bc<-bc[label_index] # order cells based on index

#d.1r = tumor.mat[,4:length(tumor.mat)]
 
#### OR #######

## Add labels without changing order - doesn't make sense
#xh2<-xh

#d.1r = d.1[,5:length(d.1)]
#bc<-colnames(d.1r)
#bc<-gsub('\\.','-',bc)
#xh$labels<-bc # add barcode back to tree
#write.table(bc,paste0('/home/sonas/copykat/data/bc_',sub,'.txt'),sep='\t',row.names=FALSE,col.names = FALSE)


#saveRDS(xh,paste0('/home/sonas/copykat/data/samples27_sub',sub,'_Hclust_ordered_by_matrix.rds'))

## Add labels by hclust order
d.1r = d.1[,4:ncol(d.1)]
d.1r = setcolorder(d.1r, label_index)
xh$labels = colnames(d.1r)

#plot(xh, labels = FALSE)


### save hclust as rds if you need it. 
#saveRDS(xh,'/home/sonas/copykat/copyKat_dendogram_full_data.rds')
saveRDS(xh,paste0('/home/sonas/copykat/data/samples27_sub',sub,'_CNA_Hclust.rds'))


## save matrix ##

saveRDS(d.1,paste0('/home/sonas/copykat/data/samples27_sub',sub,'_matrix.rds'))

# Add back copykat coord
d.1r2<-cbind(d.1[,1:3],d.1r)

saveRDS(d.1r2,paste0('/home/sonas/copykat/data/samples27_sub',sub,'_ordered_matrix.rds'))


### make clusters
# xh1 = cutree(xh, h = H_)
# my_sample_col = data.frame(cluster = xh1) 
# my_sample_col$cluster = as.factor(my_sample_col$cluster) 
# row.names(my_sample_col) <- names(xh1)

### plot heatmap using pheatmap
# png(paste0('/home/sonas/copykat/pheatmap_h',H_,'.png'))
#  PH.1.python<-pheatmap(d.1r, cluster_cols = xh, cluster_rows = TRUE, show_colnames = FALSE,
#                         color = colorRampPalette(c("blue", "white", "red"))(50), scale = 'row', border_color = NA, annotation_col = my_sample_col)
#  
#  dev.off()

 ### pick location somewhere where you can save file. 
#save_pheatmap_png(PH.1.python, paste0('/home/sonas/copykat/pheatmap_h',H_,'.png'))
                 

#### Plot heatmap ####


