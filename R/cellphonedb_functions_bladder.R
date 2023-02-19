library(Seurat)
library(SCENIC)
library(AUCell)

library(GenomicRanges)
library(plyranges)
]library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gtools)

library(ggplot2)
library(scales)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(gridExtra)
library(cowplot)
library(ggrepel)
library(lattice)
library(ComplexHeatmap)
library(corrplot)
library(stringr)
library(viridis)
library(data.table)

library(clusterProfiler)
library(DOSE)
library(enrichplot)

library(karyoploteR)

library(openxlsx)
library(readxl)
library(readr)

library(Hmisc)
library(dplyr)
library(tidyr)
#library(Xmisc)
library(purrr)
library(forcats)
library(iTALK)

#library(tidyverse)



## This function reads the cellphonedb means output to filter for uro->imm or imm-> uro cell-pairs (thus fitlering out imm->imm and uro->uro cell pairs)

# mtab: cellphonedb significant means table; dir and fname can be passed to this function instead
# dir: path to respective cellphonedb output dir
# fname: table name
# col: which column # defines the beginning of cell-pairs columns [Default: 12]

prep_table<-function(dir=NULL,fname=NULL,mtab,col='',cols_to_select)
{
if((!is.null(dir)) | !(is.null(fname)))
{
#mtab<-read.table(paste0(dir,fname),header=TRUE,stringsAsFactors = FALSE,sep = '\t')
 mtab<-fread(paste0(dir,fname)) %>% as.data.frame()
}
  

# Since rank is the last attribute column, we can select for it using which()
if(col=='')
{
col<-which(colnames(mtab)=='rank')
if(length(col)==0)
{
      col<-11
   }else{
      col<-12
    }
}
  
ftab<-mtab[,1:col]

# Since we are interested in interactions between R and NR with immune cells, only select those interaction pairs

select_r<-grep('recurrent_associated',colnames(mtab))
select_nr<-grep('non_recurrent_associated',colnames(mtab))

select<-c(select_r,select_nr)

select<-unique(select)

mtab<-mtab[,select]

# Now filter the uro->uro pairs

 str<-c('recurrent_associated|recurrent_associated',
        'non_recurrent_associated|non_recurrent_associated',
        'recurrent_associated|non_recurrent_associated',
        'non_recurrent_associated|recurrent_associated')

 rem<-vector()
 for(i in 1:length(str))
 {
 r<-which(colnames(mtab) %in% str[i])
 rem<-c(rem,r)
 
 }

 mtab<-mtab[,-rem]

ftab<-cbind(ftab,mtab)

# Reduce the sample name
colnames(ftab)<-gsub('non_recurrent_associated','NR',colnames(ftab))
colnames(ftab)<-gsub('recurrent_associated','R',colnames(ftab))

# Reduce 1 imm cell name
colnames(ftab)<-gsub('Macrophages that engulf tumor cells','MacrophageTETC',colnames(ftab))

# Remove + and spaces
colnames(ftab)<-gsub('\\+','',colnames(ftab))
colnames(ftab)<-gsub(' ','_',colnames(ftab))

#colnames(ftab)<-gsub('\\|B\\b', 'B cells',colnames(ftab))
#colnames(ftab)<-gsub('\\B|\\b', 'B cells',colnames(ftab))

return(ftab)

}


add_sample_id_suffix<-function(tab,sample,cell_types)
{
  fix_col<-colnames(tab)
  exact_sname<-paste0('\\',sample,'\\b')
  
  m<-lapply(1:length(exact_sname),function(x){
    grep(pattern =exact_sname[x] ,
         x = fix_col)})
  
  # select which samples were found
  s_found<-which(lapply(m,length)>0)
  
  for(i in 1:length(s_found))
  {
    fixed_col<-fix_col[m[[s_found[i]]]]
    
    for(j in 1:length(cell_types))
    {
    repl<-paste0(cell_types[j],'_',sample[s_found[i]])
    k<-grep(paste0('\\b',cell_types[j],'\\b'),fixed_col)
    if(length(k)>0)
      {fixed_col[k]<-gsub(cell_types[j],repl,fixed_col[k])}
    }
  # Replace actual cols
  fix_col[m[[s_found[i]]]] <-fixed_col
  }
  
  colnames(tab)<-fix_col
  
  return(tab)
  
}

merge_cols<-function(tab,colname_to_merge,cell_pair_col_start=13)
{
  # select columns to be merged
  select<-grep(colname_to_merge,colnames(tab))
  cols_to_merge<-colnames(tab)[select]
  
  # Collapse columns
  tab_selected<-tab[,select]

  new_col<-tab_selected[,1]
  for(i in 1:(ncol(tab_selected)))
  {
  j<-i+1
  repl<-which(is.na(new_col)==TRUE)
  if(length(repl)>0) # if there are NA rows, replace them
    {
     new_col[repl]<-tab_selected[repl,j]
    }
  }
  
  # Replace selected columns with collapsed single column
  new_tab<-tab[,-select]
  new_tab[,colname_to_merge]<-new_col
  
  return(new_tab)
  
  }

## This function merges cpdb output table from all samples ##
# tab_lst: list of cpdb table for all patients
# prefix:
# colnames_to_merge:
# cell_types:

merge_tables_for_all_samples<-function(tab_lst,prefix='',colnames_to_merge,cell_types)
{
  # Strip prefix
  if(prefix=='')
    {
     prefix<-'primary_recurrent_k10_'
  }
  
  names(tab_lst)<-gsub(prefix,'',names(tab_lst))
  sname<-names(tab_lst)
  suff<-paste0('_',sname)
  
  merged_tab<-full_join(x=tab_lst[[1]],tab_lst[[2]],
                          by='id_cp_interaction',suffix=c(suff[1],suff[2]))
    
  for(i in 1:((length(suff))-1))
  {
    j<-i+1
    merged_tab<-full_join(x=merged_tab,tab_lst[[j]],
                          by='id_cp_interaction',suffix=c('',suff[j]))
    
  }
  
  # Remove extra col as a result of merge (except for com_cols)
  rem_cols<-lapply(tab_lst,colnames) %>%
    do.call(c,.) %>% unique()
  
  r<-which(rem_cols %in% colnames_to_merge)
  rem_cols<-rem_cols[-r]
  
  rem<-which(colnames(merged_tab) %in% rem_cols)
  merged_tab<-merged_tab[,-rem]
  
  
  # Add suffix
  
  merged_tab<-add_sample_id_suffix(merged_tab,sample=sname,cell_types=cell_types)

  
  # Retain id_cp_interaction column
  r<-which(colnames_to_merge=='id_cp_interaction')
  colnames_to_merge2<-colnames_to_merge[-r]
  
   # Now merge redundant common cols
  for(i in 1:length(colnames_to_merge2))
  {
   merged_tab<- merge_cols(merged_tab,colname_to_merge = colnames_to_merge2[i])
  }
  
  
    
  # Sort
  mcol<-c('id_cp_interaction',colnames_to_merge2)
  others<-setdiff(colnames(merged_tab),mcol)
  
  merged_tab<-merged_tab[,c(mcol,others)]
 
  
  return(merged_tab)
  
}

## This function further filters for the specified immune clusters ##

# tab: filtered table containing imm->uro or uro->imm interactions only
# clus_name: immune cluster name to select for
# start_col: if complete table is passed, then specify the column number that marks the start of cell-pairs columns

select_imm_clusters<-function(tab,clus_name,add_cols)
{
  select<-grep(clus_name,colnames(tab),fixed = TRUE)
  
  selected_tab<-tab[,select]
  tab_add<-tab[,which(colnames(tab) %in% add_cols)]
  
  tab_add<-tab_add[,mixedsort(colnames(tab_add))]
  selected_tab<-selected_tab[,mixedsort(colnames(selected_tab))]
  
  final_tab<-cbind(tab_add,selected_tab)
  
  return(final_tab)
}


## This function splits table filtered table as uro-to and uro-from tables ##

# tab: cellphonedb table filtered to include only uro interactions to and from specific immune cluster
# to_cols:
# from_cols:

split_cpdb_filtered_tab<-function(tab,to_cols,from_cols,add_cols='')
{
  if(add_cols=='')
  {
    add_cols<-colnames(tab)[1:12]
  }
  tab_add<-tab[,which(colnames(tab) %in% add_cols)]
  tab_add<-tab_add[,mixedsort(colnames(tab_add))]
  
  to<-vector()
  
  for(i in 1:length(to_cols))
  {
    
  t<-grep(to_cols[i],colnames(tab),value = TRUE)
  
  #t1<-grep('\\|R',colnames(tab),value = TRUE)
  #t2<-grep('\\|NR',colnames(tab),value=TRUE)
  
  to<-c(to,t)
  }
  
  
  rem<-match(to,colnames(tab))
  from<-colnames(tab)[-rem]
  
  tab_to<-cbind(tab_add,tab[,to])
  tab_from<-tab[,from]
  
  tab_lst<-list(tab_to,tab_from)
  names(tab_lst)<-c('to','from')
  
    
  return(tab_lst)

}
