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


## Functions ##

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
  
  rem<-which(colnames(merged_tab) %in% c('cell_cnt','cell_pct'))
  merged_tab<-merged_tab[,-rem]
  ## Annotate ##
  
  # Remove rows where start>end
  
  # rem<-which((merged_tab$end-merged_tab$start)<0)
  # merged_tab<-merged_tab[-rem,]
  
  # Remove redundant rows
  uniq<-unique(merged_tab$abspos_start)
  merged_tab<-merged_tab[match(uniq,merged_tab$abspos_start),]
  
  # annotate main table
  anno_df<-annotate_genomic_ranges(df = merged_tab,var=var,
                                   cutoff=cutoff_l,plot=FALSE,
                                   out=out,
                                   prefix=paste0(fprefix,'_',var,'_cutoff',cutoff),
                                   save=FALSE,chr_size=chr_size)
  
  
  # annotate per clade tables
  
  anno_df_lst<-list()
  for(i in 1:length(merged_all))
  {
    anno_df_lst[[i]]<-annotate_genomic_ranges(df = merged_all[[i]],var=var,
                                              cutoff=cutoff,plot=FALSE,out=out,
                                              prefix=paste0(fprefix,'_',var,'_cutoff',
                                                            cutoff),
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
  pct.R_cp<-(R_cp/length(rclade))
  pct.NR_cp<-(NR_cp/length(rclade))
  
  
  #count_df<-cbind(R,NR,pct.R,pct.NR) %>% as.data.frame()
  
  count_df<-cbind(R,NR,pct.R,pct.NR, R_cp,NR_cp,pct.R_cp,pct.NR_cp) %>% as.data.frame()
  
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



get_concise_table<-function(df,df_anno,out,var)
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
