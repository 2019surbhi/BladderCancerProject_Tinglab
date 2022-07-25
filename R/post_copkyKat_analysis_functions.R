# Functions for CNV based analysis using copyKat outputs for scRNA data #
# Author: Surbhi Sona (Ting Lab) #


# Load Dependencies #

library('TxDb.Hsapiens.UCSC.hg38.knownGene') # hg38 gene annotation ref
library('AnnotationDbi')
library('org.Hs.eg.db') # to get gene symbol for Entrez gene.id
library('plyrange') # to apply dplyr functions on GRange obj
library(dplyr)
library(gtools)

## Function to plot genomic coordinate widths ##

# width: genomic coodinate width
# var: gain/loss
# cutoff: gain or loss cutoff
# bw: binwidth for histogram
# out: output dir
# prefix: fileame prefix
# save: if TRUE plot is saved else ggplot obj returned

plot_genomic_coord_width<-function(width,var,cutoff,bw=200,out,prefix,save=FALSE)
{
  
width<-as.data.frame(width)

g<-ggplot(width,aes(x=width))+
  geom_histogram(binwidth = bw) +
  xlab('genomic coodindate size') +
  ggtitle(label = 'Genomic coordinate size distribution',
          subtitle = paste0('[',var, 'cutoff=',cutoff,']'))+
  labs(caption = paste0('Range: ',paste0('[',min(width$width),',',max(width$width),']'))) +
  scale_x_continuous(breaks = round(seq(200000, max(width$width), by = 10000),digits = 0)) +
  theme(plot.caption=element_text(color='navy blue',size=10)) +
  theme_bw()
if(save==FALSE)
{return(g)}else
  {
    ggsave(paste0(out,prefix,'_genomic_coord_width_',
                  var,'_cutoff',cutoff,'.png'),g)

  }
}


## Function to intersect genomic intervals with hg38 (or any other) ref ##

# gr: genomic range obj created using gain/loss dataframe
# ref: Annotation reference database to use [Default: TxDb.Hsapiens.UCSC.hg38.knownGene]

add_gene_anno<-function(gr,ref='TxDb.Hsapiens.UCSC.hg38.knownGene')
{

# Create selected ref database
gr_db<-genes(keepStandardChromosomes(ref))

#gr_db_ts<- transcripts(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene))

# Check overlaps

# returns ranges of gene interval of intersecting ranges
#sub<-subsetByOverlaps(gr_db,gr)

int<-join_overlap_intersect(x=gr_db,y=gr)
int<-sort(int)

#int_ts<-join_overlap_intersect(x=gr_db_ts,y=gr)
#int_ts<-sort(int_ts)

# Some genes span >1 interval and these duplicate entries need to be removed
int<-int[unique(int$gene_id),]

# Get annotations for gene symbol
gene_id<-int$gene_id
#ts_id<-as.character(int_ts$tx_id)

anno<-AnnotationDbi::select(org.Hs.eg.db, keys=gene_id, columns=c('SYMBOL',"GENENAME","GENETYPE"), keytype='ENTREZID')

# keytypes(org.Hs.eg.db)
# anno_ts<-AnnotationDbi::select(org.Hs.eg.db, keys=ts_id, columns='SYMBOL', keytype='ENSEMBLTRANS')


int_df<-as.data.frame(int)
int_anno_df<-merge(int_df,anno,by.x="gene_id",by.y='ENTREZID')

return(int_anno_df)

}


## Function to add ends per chr ###

# df_by_chr: dataframe containing genomic intervals from copykat table subset by chr
# chr_size : chr size table for matching ref

add_ends<-function(df_by_chr,chr_size)
    {
  
    # Generally chromosome column is named as 'chr' but copyKat output uses 'chrom'; GRange uses 'seqname' so need to ensure this function is compatible for all inputs
    start_col<-match(c('chr_pos','chrompos'),colnames(df_by_chr))
    start_col<-start_col[complete.cases(start_col)]
    
    
    # Add chr end (start of next interval-1)
    start<-df_by_chr[,chr_col]
    end<-start[2:length(start)]
    end<-end-1
    
    chr_col<-match(c('chr','chrom','seqname'),colnames(df_by_chr))
    chr_col<-chr_col[complete.cases(chr_col)]
    
    # Add chr last segment based on chr size
    i<-match(unique(df_by_chr[,chr_col]),chr_size$chr)
    end[length(end)+1]<-chr_size$size[i]
    
    # Add chr based on res
    #end[length(end)+1]<-start[length(start)]+res
    
    df_by_chr$end<-end
     if('chrompos' %in% colnames(df_by_chr))
    {
    colnames(df_by_chr)<-gsub('chrompos','start', colnames(df_by_chr))
    }
    
    
    if('chrom' %in% colnames(df_by_chr))
    {
    colnames(df_by_chr)<-gsub('chrom','chr', colnames(df_by_chr))
    }else if('seqname' %in% colnames(df_by_chr))
        {
         colnames(df_by_chr)<-gsub('seqname','chr', colnames(df_by_chr))
        }
    
    # sort columns
    
    odr<-c('chr','start','end','abspos')
    cols<-colnames(df_by_chr)

    #add back any other remaining columns in the dataframe
    rem_col<-setdiff(cols,ord)
    ord<-c(ord,rem_col)
  
    df_by_chr<-df_by_chr[,odr]
    return(df_by_chr)
  }

## Function to add abspos ends per chr ###

# df_by_chr: dataframe containing genomic intervals from copykat table subset by chr
# chr_size : chr size table (must have abspos size too)

add_abspos_ends<-function(df_by_chr,chr_size)
{
  {
    # Generally chromosome column is named as 'chr' but copyKat output uses 'chrom'; GRange uses 'seqname' so need to ensure this function is compatible for all inputs
    chr_col<-match(c('chr','chrom'),colnames(df_by_chr))
    chr_col<-chr_col[complete.cases(chr_col)]
    
    # Add abspos end (start of next interval-1)
    abspos_start<-df_by_chr$abspos
    abspos_end<-abspos_start[2:length(abspos_start)]
    abspos_end<-abspos_end-1
    
    # Add abspos last segment based on abspos size
    i<-match(unique(df_by_chr[,chr_col]),chr_size$chr)
    abspos_end[length(abspos_end)+1]<-chr_size$abs_size[i]
    
    #Add abspos end to the table
    df_by_chr$abspos_end<-abspos_end
    
     colnames(df_by_chr)<-gsub('\\babspos\\b','abspos_start', colnames(df_by_chr))
     
     if('chrom' %in% colnames(df_by_chr))
     {
     colnames(df_by_chr)<-gsub('chrom','chr', colnames(df_by_chr))
     }else if('seqname' %in% colnames(df_by_chr))
         {
          colnames(df_by_chr)<-gsub('seqname','chr', colnames(df_by_chr))
         }
     
    cols<-colnames(df_by_chr)
    ord<-c('chr','start','end','abspos_start','abspos_end')
    
    #add any other remaining columns in the dataframe
    rem_col<-setdiff(cols,ord)
    ord<-c(ord,rem_col)
    
    df_by_chr<-df_by_chr[,ord]
    return(df_by_chr)
  }
}

## Function to process gain-loss outputs ##

# df: dataframe containing gain/loss info (chrom,chrompos, abspos)
# var: gain or loss
# plot: if TRUE genomic widths are plotted
# bw: bindwidth for genomic coordinate range plot
# out: output dir (to save the plot in)
# prefix: filename prefix
# save: if FALSE, the

annotate_genomic_ranges<-function(df,var,cutoff,plot=FALSE,bw=200,out='./',prefix='',save=FALSE,chr_size)
{

# Add range and width as additional columns
width<-df$end-df$start
df$width<-width
df$gi_start<-df$start
df$gi_end<-df$end

if(plot==TRUE)
{
  plot_genomic_coord_width(width,var,cutoff,bw,out,prefix,save=TRUE)
}

# Create GRange obj
rownames(df2)<-NULL


# Rename chr23
if('chr23' %in% df$chr)
{
  idx<-which(df$chr=='chr23')
  df$ch[idx]<-'chrX'
}
gr<-makeGRangesFromDataFrame(df = df,keep.extra.columns = TRUE)

# Add gene symbol
gr_anno<-add_gene_anno(gr)

return(gr_anno)

}


### Functions related to subsetting copyKat genomic intervals based on a cutoff ###

## Function to generate count of cells with CNV values > cutoff ##

# rec_mat: matrix with CNV values for all cells (columns) for all genomic intervals (rows) generated for a given run
# cutoff: gain/loss cutoff for subsetting
# var: specify whether the variation is 'gain' or 'loss'

cell_count_by_cutoff<-function(rec_mat,cutoff,var)
{

# For each row of input matrix, this function calculates number of cells having values> or < cutoff

cell_count<-vector()

if(var=='loss')
{
for(i in 1:nrow(rec_mat))
  {
  len<-length(which(rec_mat[i,]<= (cutoff)))
  cell_count<-append(cell_count,len)
  }
}else if(var=='gain')
{
  for(i in 1:nrow(rec_mat))
  {
    len<-length(which(rec_mat[i,]>= (cutoff)))
    cell_count<-append(cell_count,len)
  }
}else
{
  cat('please specify either loss or gain \n')
}

return(cell_count)
}


## Function to apply a cell count filter for subsets of genomic intervals ##

# cell_count: cell count for cells that passed user defined threshold/cutoff
# tot_cells: total cell count in the given matrix

subset_indx<-function(cell_count,tot_cells)
{


i75<-which(cell_count>=(0.75*tot_cells))
#i70<-which(cell_count>=(0.70*tot_cells))
#i60<-which(cell_count>=(0.60*tot_cells))
#i50<-which(cell_count>=(0.50*tot_cells))
#i25<-which(cell_count>=(0.25*tot_cells))

#idx_list<-list(i75,i70,i60,i50,i25)

#return(idx_list)
return(i75)

}

## Function to generate count of cells with CNV values > cutoff ##

# rec_mat2: matrix with CNV values for all cells (columns) for all genomic intervals (rows) generated for a given run
# usr_cutoff: gain/loss cutoff for subsetting
# var: specify whether the variation is 'gain' or 'loss'
# out: output directory path
# fname: filename prefix

subset_by_cutoff<-function(rec_mat2,usr_cutoff,var,out,fname)
{
rec_mat<-rec_mat2[,-(1:5)]
cols<-ncol(rec_mat)

 if(var=='gain')
{
cc<-cell_count_by_cutoff(rec_mat,cutoff=usr_cutoff,var)
g_idx<-subset_indx(cell_count=cc,tot_cells=cols)
gain<-lapply(1:length(g_idx),function(x){return(rec_mat2[g_idx[[x]],1:5])})
#names(gain)<-c('75g','70g','60g','50g','25g')

saveRDS(object = gain,paste0(out,sname,'_cutoff',usr_cutoff,'_gain.rds'))
}else if(var=='loss')
{
 cc<-cell_count_by_cutoff(rec_mat,cutoff=usr_cutoff,var)
 l_idx<-subset_indx(cell_count=cc,tot_cells=cols)
 loss<-lapply(1:length(l_idx),function(x){return(rec_mat2[l_idx[[x]],1:5])})
 #names(loss)<-c('75l','70l','60l','50l','25l')
 saveRDS(object = loss,paste0(out,fname,'_cutoff',usr_cutoff,'_loss.rds'))
}else
 {
  cat('you need to specify variation as either gain or loss \n')
 }
}


