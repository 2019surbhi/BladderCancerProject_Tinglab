library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)

out<-'/home/sonas/cellphonedb/cpdb_output100/'

#dir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/cellphonedb/group_comparison/cpdb_outputs100/'

dir<-'/home/sonas/cellphonedb/cpdb_outputs100/'

prefix_n<-'Naive'
prefix_r<-'Recurrent'

# pval files
pvalues_n<-paste0(dir,prefix_n,'/pvalues.txt')
pvalues_r<-paste0(dir,prefix_r,'/pvalues.txt')

# meta files
#inputdir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/cellphonedb/group_comparison/cpdb_inputs100/'

inputdir<-'/home/sonas/cellphonedb/cpdb_input100/'

meta_n<-paste0(inputdir,prefix_n,'/Naive_metadata.txt')
meta_r<-paste0(inputdir,prefix_r,'/Recurrent_metadata.txt')


# Function from original cellphonedb repo at https://github.com/Teichlab/cellphonedb/edit/master/cellphonedb/src/plotters/R/plot_heatmaps.R

get_count_mat<-function(metafile,meta_sep='\t',pvalues_file,pvalues_sep='\t', pvalue=0.05)
{
  
  # Extract interaction pairs
  
  meta = read.csv(metafile, comment.char = '', sep=meta_sep)
  
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, sep=pvalues_sep, comment.char = '', check.names = F)
  
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)]
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1_all = unique(meta[,2])
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  # count interactions
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  # Create data matrix
  count_matrix = matrix(count1, nrow=length(unique(meta[,2])), ncol=length(unique(meta[,2])))
  rownames(count_matrix)= unique(meta[,2])
  colnames(count_matrix)= unique(meta[,2])
  return(count_matrix)
}

print_cpdb_heatmap<-function(count_matrix,cluster_cols = T,cluster_rows = T,main='',split='none',fname,w=10,h=10,col_breaks='auto')
{
  # Initialize colors
  #col1<-'dodgerblue4'
  #col2<-'peachpuff'
  #col3<-'deeppink4'
  
  col1<-'blue'
  col2<-'white'
  col3<-'red'
  
  hcol<-colorRampPalette(c(col1,col2,col3 ))(1000)
  
  if(col_breaks!='auto')
  {
    col_fun<-colorRamp2(col_breaks,hcol)
  }else{
    col_fun<-hcol
  }
  
  #col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
  
  # Print heatmap
  
  if(split=='none')
    
  {pdf(fname,width=w,height=h)
    
    draw(Heatmap(count_matrix, col= col_fun,
                 show_heatmap_legend = TRUE,
                 name = main ,
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_cols,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE,
                 rect_gp = gpar(col = "white", lwd = 2)))
    dev.off()   
  }else{   
    pdf(fname,width=w,height=h)
    
    draw(Heatmap(count_matrix, col= col_fun,
                 show_heatmap_legend = TRUE,
                 name = main ,
                 cluster_rows = cluster_rows,
                 cluster_columns = cluster_cols,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE,
                 column_split = split,
                 row_split = split,
                 rect_gp = gpar(col = "white", lwd = 2)))
    dev.off()   
  }
  # pheatmap(count_matrix, show_rownames = TRUE, 
  #          show_colnames = TRUE, scale="none", 
  #          cluster_cols = cluster_cols,
  #          border_color='white', cluster_rows = cluster_rows,
  #          fontsize_row = 11, fontsize_col = 11,
  #          main = main, treeheight_row = 0,color = col.heatmap, 
  #          treeheight_col = 0)
  
} 


#### Heatmap ###

col_breaks<-c(
  seq(0,50,length.out=200),
  seq(50,100,length.out=200),
  seq(100,150,length.out=250),
  seq(150,200,length.out=250),
  seq(200,300,length.out=75),
  seq(250,300,length.out=25)
)

out_plot<-'/home/sonas/cellphonedb/cpdb_current/plots/'

# Step 1 get both count matrices

# Get interaction counts
counts_n<-get_count_mat(metafile=meta_n,meta_sep='\t',pvalues_file=pvalues_n,pvalues_sep='\t', pvalue=0.001)

counts_r<-get_count_mat(metafile=meta_r,meta_sep='\t',pvalues_file=pvalues_r,pvalues_sep='\t', pvalue=0.001)

# Since Naive has fewer cell pairs, subset the same for recurrent
counts_r2<-counts_r[rownames(counts_n),colnames(counts_n)]


# Create diff matrix (R-N)

identical(colnames(counts_n),colnames(counts_r2))
identical(rownames(counts_n),rownames(counts_r2))

diff_mat<-counts_r2-counts_n

# Find total interaction (abs colsum)
diff_mat_abs<-abs(diff_mat)
total<-colSums((diff_mat_abs))

#Scale diff matrix
diff_mat_scaled<-diff_mat/total

# Split by rows

odr<-c(1:6,25,14,7:13,15:23,26,24)
diff_mat_scaled_ordered<-diff_mat_scaled[odr,odr]
group<-c(rep('uro',8),rep('immune',17),'stromal')
grp<-factor(group,
            levels = c('uro','immune','stromal'),
            ordered = TRUE)


fname<-paste0(out_plot,'Recurrent_minus_Naive_cpdb_heatmap_scaled_clustered.pdf')

print_cpdb_heatmap(count_matrix = diff_mat_scaled_ordered,cluster_cols = T,cluster_rows = T,main='scaled (R-N)',fname=fname, split = grp)


fname1<-paste0(out_plot,'Recurrent_minus_Naive_cpdb_heatmap.pdf')

print_cpdb_heatmap(count_matrix = diff_mat,cluster_cols = T,cluster_rows = T,main='scaled (R-N)',fname=fname1, split = grp)
