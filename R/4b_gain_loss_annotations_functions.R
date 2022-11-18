# Functions for CNV based analysis using copyKat outputs for scRNA data #
# Author: Surbhi Sona (Ting Lab) #


# Load Dependencies #

library(GenomicRanges)
library(plyranges)
#library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gtools)

library(openxlsx)
library(readxl)
library(data.table)
library(readr)
library(ComplexHeatmap)

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

add_gene_anno<-function(gr)
{
  # Following libraries must be loaded
  #p<-c('TxDb.Hsapiens.UCSC.hg38.knownGene',
  #     'AnnotationDbi','org.Hs.eg.db',plyrange)
  
  # Create ref
  gr_db<-genes(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  #gr_db_ts<- transcripts(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  # Check overlaps
  
  # returns ranges of gene interval of intersecting ranges
  #sub<-subsetByOverlaps(gr_db,gr) 
  
  int<-join_overlap_intersect(x=gr_db,y=gr)
  int<-sort(int)
  
  #Save a copy of intersect table
  bc<-int %>% as.data.frame(row.names = NULL)
  
  #int_ts<-join_overlap_intersect(x=gr_db_ts,y=gr)
  #int_ts<-sort(int_ts)
  
  # Some genes span >1 interval and these duplicate entries need to be removed else genes() yields error
  
  int<-int[unique(int$gene_id),]
  
  # Get annotations for gene symbol
  gene_id<-int$gene_id
  #ts_id<-as.character(int_ts$tx_id)
  
  anno<-AnnotationDbi::select(org.Hs.eg.db, keys=gene_id, columns=c('SYMBOL',"GENENAME"), keytype='ENTREZID')
  
  # keytypes(org.Hs.eg.db)
  # anno_ts<-AnnotationDbi::select(org.Hs.eg.db, keys=ts_id, columns='SYMBOL', keytype='ENSEMBLTRANS')
  
  int_anno_df<-merge(bc,anno,by.x="gene_id",by.y='ENTREZID')
  
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

annotate_genomic_ranges<-function(df,var,cutoff,plot=FALSE,bw=200,out='./',prefix='',save=TRUE,chr_size)
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
  rownames(df)<-NULL
  
  # Rename chr23
  if('chr23' %in% df$chr)
  {
    idx<-which(df$chr=='chr23')
    df$chr[idx]<-'chrX'
  }
  gr<-makeGRangesFromDataFrame(df = df,keep.extra.columns = TRUE)
  
  # Add gene symbol
  gr_anno<-add_gene_anno(gr)
  
  return(gr_anno)
  
}

## Function to process gain-loss outputs per clade ##

# thresh_dir: path to dir where threholds are stored
# clade: vector with clades eg. calde<-c(1,2,3,4,5,6,7,8,9,10)
# cutoff: CNV cutoff (used to retrieve threshold table/list from threhold directory)
# var: specify variation is 'gain' or 'loss'
# p: cell percent cutoff (e.g. 70)
# chr_sizes: table with chr size where chr annotation matches the chr annotation in other tables i.e. '1' vs 'chr1'
# gene_tab: copykay gene table to check the overlap with
# rclade: vector specifying recurrent clades
# nrclade: vector specifying non-recurrent clades
# fprefix: file prefix to use for saved file
# out: path to output dir

get_per_clade_annotated_table<-function(thresh_dir,clade,cutoff,var,p,chr_sizes,gene_tab,rclade,nrclade,fprefix,out)
{
  df_lst<-lapply(1:length(clade),function(x){readRDS(paste0(thresh_dir,fprefix,'_',clade[x],'_cutoff',cutoff,'_',var,'.rds'))})
  
  
  if(var=='gain')
  {p_name<-paste0(p,'g')}else{
    p_name<-paste0(p,'l')
  }
  
  p_idx<-grep(p_name,names(df_lst[[1]]))
  
  df_lst2<-lapply(df_lst,'[[',p_idx)
  
  names(df_lst2)<-clade
  
  write.xlsx(df_lst2,paste0(out,fprefix,'_per_clade_genomic_intervals_',var,'_',p,'p_cutoff',cutoff,'.xlsx'))
  
  # Filter clades that have 0 rows
  
  rem<-which((unlist(lapply(df_lst2,nrow)))==0)
  df_lst2_fil<-df_lst2[-rem]
  #df_l_lst80_fil<-df_l_lst80
  
  # Merge all tables
  merged_all1<-lapply(df_lst2_fil,rbind)
  
  # remove row where start>end - why do we have these?
  merged_all<-lapply(1:length(merged_all1),function(x){rem<-which((merged_all1[[x]]$end-merged_all1[[x]]$start)<0);
  if(length(rem)>0)
  {merged_all1[[x]]<-merged_all1[[x]][-rem,]};
  return(merged_all1[[x]])})
  
  names(merged_all)<-names(merged_all1)
  
  #merged_tab<-full_join(x = merged_all[[1]],y = merged_all[[2]],by=c("chr","start","end","abspos_start","abspos_end"),suffix=c('',''))
  
  merged_tab<-full_join(x = merged_all[[1]],y = merged_all[[2]],by=c("chr","start","end","abspos_start","abspos_end","cell_cnt","cell_pct"),suffix=c('',''))
  
  j<-2
  for(i in 1:(length(merged_all)-2))
  {
    j<-j+1
    #merged_tab<-full_join(x=merged_tab,y = merged_all[[j]],by=c("chr","start","end","abspos_start","abspos_end"),suffix=c('',''))
    
  merged_tab<-full_join(x=merged_tab,y = merged_all[[j]],by=c("chr","start","end","abspos_start","abspos_end","cell_cnt","cell_pct"),suffix=c('',''))
    
  }
  
  ## Annotate ##
  
  # Remove rows where start>end
  
  # rem<-which((merged_tab$end-merged_tab$start)<0)
  # merged_tab<-merged_tab[-rem,]
  
  # annotate main table
  anno_df<-annotate_genomic_ranges(df = merged_tab,var='loss',
                                   cutoff=cutoff_l,plot=FALSE,
                                   out=out,
                                   prefix=paste0(fprefix,'_loss_cutoff',cutoff_l),
                                   save=FALSE,chr_size=chr_size)
  
  
  # annotate per clade tables
  
  anno_df_lst<-list()
  for(i in 1:length(merged_all))
  {
    anno_df_lst[[i]]<-annotate_genomic_ranges(df = merged_all[[i]],var=var,
                                              cutoff=cutoff_l,plot=FALSE,out=out,
                                              prefix=paste0(fprefix,'_',var,'_cutoff',
                                                            cutoff_l),
                                              save=FALSE,chr_size=chr_size)
    
  }
  
  names(anno_df_lst)<-names(merged_all)
  
  
  write.xlsx(anno_df_lst,paste0(out,fprefix,'_per_clade_annotated_',var,'_',p,'p_cutoff',cutoff,'.xlsx'))
  
  # Mark/count clades that contain a given segement from the merged table
  
  clades<-list()
  clades<-lapply(1:length(anno_df_lst), function(x){cd<-rep(x=0,nrow(anno_df));
  return(cd)})
  
  for(i in 1:length(anno_df_lst))
  {
    cd<-vector()
    cd<-clades[[i]]
    comm<-which(anno_df$SYMBOL %in% anno_df_lst[[i]]$SYMBOL)
    cd[comm]<-1
    clades[[i]]<-cd
  }
  
  
  names(clades)<-names(anno_df_lst)
  
  clade_df<-clades[[1]] %>% as.data.frame()
  for(i in 2:(length(clades)))
  {
    clade_df<-cbind(clade_df,clades[[i]])
  }
  
  colnames(clade_df)<-names(anno_df_lst)
  
# Cell count
  cell_pct<-list()
  cell_pct<-lapply(1:length(anno_df_lst), function(x){cp<-rep(x=0,nrow(anno_df));
  return(cp)})
  
  for(i in 1:length(anno_df_lst))
  {
    cp<-vector()
    cp<-cell_pct[[i]]
    df_idx<-which(anno_df$SYMBOL %in% anno_df_lst[[i]]$SYMBOL)
    #df_lst_idx<-which(anno_df_lst[[i]]$SYMBOL %in% anno_df$SYMBOL)
    cp[df_idx]<-anno_df_lst[[i]]$cell_pct
    cell_pct[[i]]<-cp
  }
  
  
  names(cell_pct)<-names(anno_df_lst)
  
  cell_pct_df<-cell_pct[[1]] %>% as.data.frame()
  for(i in 2:(length(cell_pct)))
  {
   cell_pct_df<-cbind(cell_pct_df,cell_pct[[i]])
  }
  
  colnames(cell_pct_df)<-paste0('cell_pct_',names(anno_df_lst))
  
  
 # anno_df2<-cbind(anno_df,clade_df)
  
 anno_df2<-cbind(anno_df,clade_df,cell_pct_df)

  
  
  # Remove ambiguous gene entry (has geneID but not gene name/symbol)
  rem<-which(is.na(anno_df2$SYMBOL)==TRUE)
  anno_df2<-anno_df2[-rem,]
  
  ## Add R and NR count and % ##
  
  rcol<-match(rclade,colnames(anno_df2))
  R<-rowSums(anno_df2[,rcol])
  nrcol<-match(nrclade,colnames(anno_df2))
  NR<-rowSums(anno_df2[,nrcol])
  pct.R<-(R/length(rclade)) *100
  pct.NR<-(NR/length(rclade))*100
  
  
   rpcol<-match(paste0('cell_pct_',rclade),colnames(anno_df2))
   R_cp<-rowSums(anno_df2[,rpcol])
   nrpcol<-match(paste0('cell_pct_',nrclade),colnames(anno_df2))
   NR_cp<-rowSums(anno_df2[,nrpcol])
   pct.R_cp<-(R_cp/length(rclade)) *100
   pct.NR_cp<-(NR_cp/length(rclade))*100


  #count_df<-cbind(R,NR,pct.R,pct.NR) %>% as.data.frame()
  
  count_df<-cbind(R,NR,pct.R,pct.NR, R_cp,NR_cp,pct.R_cp,pct.NR_cp) %>% as.data.frame()
   
  count_df$pct.R<-paste0(count_df$pct.R,'%')
  count_df$pct.NR<-paste0(count_df$pct.NR,'%')
  count_df$pct.R_cp<-paste0(count_df$pct.R_cp,'%')
  count_df$pct.NR_cp<-paste0(count_df$pct.NR_cp,'%')
 
  anno_df2<-cbind(anno_df2,count_df)
  
  ## Add overlap column - overlap with old table ## - no longer need it!
  
  
  # if(anno_df_old!='')
  # {
  #   # Intersect to get common genes #
  #   comm<-intersect(anno_df_old$SYMBOL,anno_df2$SYMBOL)
  #   
  #   anno_df2$overlap_w_old<-rep(0,nrow(anno_df2))
  #   idx<-which(anno_df2$SYMBOL %in% comm)
  #   anno_df2$overlap_w_old[idx]<-1
  # }
  # 
  
  # overlap with copykat table
  
  ov<-rep(0,nrow(anno_df2))
  int<-intersect(gene_tab$hgnc_symbol,anno_df2$SYMBOL)
  ov[which(anno_df2$SYMBOL %in% int)]<-1
  anno_df2$overlap_w_ck_gene<-ov
  
  
  write.xlsx(anno_df2,paste0(out,fprefix,'_per_clade_merged_',var,'_',p,'p_cutoff',cutoff,'_final.xlsx'))
  
  
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
p<-(cell_count/tot_cells)*100
p<-round(p,digits=0)

pvec<-c(90,80,75,70)

tab_lst<-list()
for(i in 1:length(pvec))
{
idx<-which(p>=pvec[i])
tab_lst[[i]]<-data.frame("seg_indx"=idx,"cell_cnt"=cell_count[idx],"cell_pct"=p[idx])
  
}
return(tab_lst)

#i90<-which(p>=90)
#i80<-which(p>=80)
#i75<-which(p>=75)
#i70<-which(p>=70)

#i90<-which(cell_count>=(0.90*tot_cells))
#i80<-which(cell_count>=(0.80*tot_cells))
#i75<-which(cell_count>=(0.75*tot_cells))
#i70<-which(cell_count>=(0.70*tot_cells))

#idx_list<-list(i90,i80,i75,i70)

#return(idx_list)

}


## Function to generate count of cells with CNV values > cutoff ##

# rec_mat2: matrix with CNV values for all cells (columns) for all genomic intervals (rows) generated for a given run
# usr_cutoff: gain/loss cutoff for subsetting
# var: specify whether the variation is 'gain' or 'loss'
# out: output directory path
# fname: filename prefix

subset_by_cutoff<-function(rec_mat,usr_cutoff,var,out,sname)
{
#rec_mat<-rec_mat2[,-(1:3)]
cols<-ncol(rec_mat)

 if(var=='gain')
{
cc<-cell_count_by_cutoff(rec_mat,cutoff=usr_cutoff,var)
g_tab<-subset_indx(cell_count=cc,tot_cells=cols)
g_idx<-sapply(g_tab,'[[',1)
cell_cnt<-sapply(g_tab,'[[',2)
cell_pct<-sapply(g_tab,'[[',3)
#gain<-lapply(1:length(g_idx),function(x){return(rec_mat2[g_idx[[x]],1:3])})
gain<-lapply(1:length(g_idx), function(x){final_tab<-(gi_ref[g_idx[[x]],]);
                                          final_tab$cell_cnt<-cell_cnt[[x]];
                                          final_tab$cell_pct<-cell_pct[[x]];
                                          return(final_tab)})
names(gain)<-c('90g','80g','75g','70g')
saveRDS(object = gain,paste0(out,sname,'_cutoff',usr_cutoff,'_gain.rds'))
}else if(var=='loss')
{
 cc<-cell_count_by_cutoff(rec_mat,cutoff=usr_cutoff,var)
 l_tab<-subset_indx(cell_count=cc,tot_cells=cols)
 l_idx<-sapply(l_tab,'[[',1)
 cell_cnt<-sapply(l_tab,'[[',2)
 cell_pct<-sapply(l_tab,'[[',3)
 #loss<-lapply(1:length(l_idx),function(x){return(rec_mat2[l_idx[[x]],1:3])})
 loss<-lapply(1:length(l_idx), function(x){final_tab<-(gi_ref[l_idx[[x]],]);
                                          final_tab$cell_cnt<-cell_cnt[[x]];
                                          final_tab$cell_pct<-cell_pct[[x]];
                                          return(final_tab)})
 names(loss)<-c('90l','80l','75l','70l')
 saveRDS(object = loss,paste0(out,sname,'_cutoff',usr_cutoff,'_loss.rds'))
  }else
 {
  cat('you need to specify variation as either gain or loss \n')
 }
}


### Archived Functions ###

intersect_by_abs<-function(df,ref_tab,out,var,cell_cutoff,cnv_cutoff,save=FALSE)
{
df_anno<-data.frame()
not_found<-data.frame()

for(i in 1:nrow(df))
{

# Get genes within i(th) interval
sub<-ref_tab[dplyr::between(ref_tab$gene_abspos,df$abspos_start[i],df$abspos_end[i]),]

if(nrow(sub)>0)
{
#repeat gi rows
n<-nrow(sub)
if(n>1)
{
gi<-df[rep(i, times=n), ]
}else{
  gi<-df[i,]
}
temp<-cbind(gi,sub)
df_anno<-rbind(df_anno,temp)
}else{
 gi<-df[i,]
 not_found<-rbind(not_found,gi)
}

}

gi_count<-df_anno %>% group_by(chr,start,end) %>% tally
gi_count<-gi_count %>% arrange(.,chr,start,end)

nf_count<-not_found %>% group_by(chr,start,end) %>% tally()
nf_count$n<-0
nf_count<-nf_count %>% arrange(.,chr,start,end)

# Merge the tables
count<-rbind(gi_count,nf_count)
count<-count %>% arrange(.,chr,start,end)

# Create concise anno table
concise_anno<-get_concise_table(df=df,df_anno=df_anno,
                                ref_tab,out,var,cell_cutoff,cnv_cutoff)


concise_anno<- concise_anno %>% arrange(.,chr,start,end)
concise_anno2<-merge(concise_anno,count,by.x=c('chr','start','end'),by.y=c('chr','start','end'))

concise_anno2<-concise_anno2 %>% arrange(chr,start,end)

if(save==TRUE)
{

write.xlsx(concise_anno2,paste0(out,'concise_annotated_genomic_ranges_',var,'_',cell_cutoff,'_p_cutoff',cnv_cutoff,'_filtered_by_copyKat_gene_tab_by_abspos.xlsx'))

}else{
  return(concise_anno2)
}

}


## Function to create concise annotated table ##

# df: dataframe containing genomic interval info (chr, start, end,abspos_start, abspos_end)
# df_anno: dataframe containing annotated genomic interval info (chr, start, end, abspos_start, abspos_end,gene_abspos,gene_chr,gene_start, gene_end, ensembl_gene_id, hgnc_symbol,band)
# ref_tab: copyKat gene tab (create non-redundant union of genes if > 1 sample by merging individual tables)
# var: specify 'gain' or 'loss'
# cell_cutoff: specify the cell cutoff used for create df
# cnv_cutoff: specify copyKat CNV value cutoff used for creating df

get_concise_table<-function(df,df_anno,ref_tab,out,var,cell_cutoff,cnv_cutoff)
{

gene_id<-list()
gene<-list()

for(i in 1:nrow(df))
{
j1<-which(df_anno$abspos_start %in%df$abspos_start[i])
if(length(j1)>0)
{
gene_id[[i]]<-df_anno$ensembl_gene_id[j1]
gene[[i]]<-df_anno$hgnc_symbol[j1]
}else{
  gene_id[[i]]<-'NA'
  gene[[i]]<-'NA'
}

}

anno_concise<-df
anno_concise$gene<-gene
anno_concise$gene_id<-gene_id

return(anno_concise)
}

