library(dplyr)
library(data.table)
library(tidyverse)



## Function #####

select_cpdb_rows_cols<-function(sig_mean,cells,fprefix,mean_cutoff=1.5)
{
  sig<-fread(sig_mean) %>% as.data.frame()
  
  for(i in 1:length(cells))
  {
    s<-grep(cells[i],colnames(sig))
    sig_sub<-sig[,s]
    
    select_cols<-colnames(sig_sub)
        
    # Create row max
    sig_sub<-replace(sig_sub, is.na(sig_sub), 0)
    sig_sub$max<-apply(sig_sub,1,max)
    
    sig_sub<-cbind(sig[,1:12],sig_sub)
    rows<-which(sig_sub$max>=mean_cutoff)
    
    select_rows<-sig_sub$interacting_pair[rows]
  }
  
  return(list(select_rows,select_cols))
}


select_cpdb_common_rows_cols<-function(sig_mean_r,sig_mean_n,cells,diff_cutoff=1,sub=NULL)
{
  sign<-fread(sig_mean_n) %>% as.data.frame()
  sigr<-fread(sig_mean_r) %>% as.data.frame()
  
  for(i in 1:length(cells))
  {
    sn<-grep(cells[i],colnames(sign))
    sign_sub<-sign[,sn]
    
    sr<-grep(cells[i],colnames(sigr))
    sigr_sub<-sigr[,sr]
    
    # subset common cols
    
    comm<-intersect(colnames(sign_sub),colnames(sigr_sub))
    comm<-sort(comm)
    from<-comm[grep(paste0(cells[i],'|'),comm,fixed = TRUE)]
    to<-comm[grep(paste0('|',cells[i]),comm,fixed = TRUE)]
    
    if(is.null(sub)==FALSE)
    {
      f<-strsplit(from,split='\\|') %>% sapply(.,'[[',2)
      from<-from[which(f %in% sub)]
      
      t<-strsplit(to,split='\\|') %>% sapply(.,'[[',1)
      to<-to[which(t %in% sub)]
      
    }
    
    comm2<-c(from,to)
    
    sign_sub<-sign_sub[,comm2]
    sigr_sub<-sigr_sub[,comm2]
    
      
    # Create row max
    sign_sub<-replace(sign_sub, is.na(sign_sub), 0)
    sigr_sub<-replace(sigr_sub, is.na(sigr_sub), 0)
    
    # Max diff
    sign_sub2<-cbind(sign[,1:12],sign_sub)
    sigr_sub2<-cbind(sigr[,1:12],sigr_sub)
    
    comm_int<-intersect(sign$interacting_pair,sigr$interacting_pair)
    n<-match(comm_int,sign_sub2$interacting_pair)
    r<-match(comm_int,sigr_sub2$interacting_pair)
    
    sig_sub<-sign_sub2[n,13:ncol(sign_sub2)]-sigr_sub2[r,13:ncol(sigr_sub2)]
    sig_sub$max_diff<-apply(sig_sub,1,max)
    sig_sub$min_diff<-apply(sig_sub,1,min)
    sig_sub<-cbind(comm_int,sig_sub)
    
    rowsn<-which(sig_sub$max>=diff_cutoff)
    rowsr<-which(abs(sig_sub$min)>=diff_cutoff)
    rows<-c(rowsn,rowsr)
    
    select_rows<-sig_sub$comm_int[rows]
    }
  return(list(select_rows,comm2))
}


### Output ###
out<-'/home/sonas/cellphonedb/cpdb_current/dotplots/'
# cells<-c('MMDSC','Type 2 conventional dendritic',
#                 'Treg','Macrophages','CD4 T Central Memory 2',
#                 'Proliferating lymphocytes')

cells<-c('CD8 T effector')

naive<-fread('/home/sonas/cellphonedb/cpdb_outputs100/Naive/means.txt')
naive<-naive[,-c(1:11)]
cols_n<-colnames(naive)
all_cells_n<-cols_n %>% strsplit(cols_n,split='\\|') %>% sapply(.,'[[',1) %>% unique()

uro_n<-all_cells_n[c(10:14,22,25)]

# Same set of cells so no need to repeat it for recurrent



## Naive ##
sig_mean_n<-'/home/sonas/cellphonedb/cpdb_outputs100/Naive/significant_means.txt'
fprefix_n<-'naive'
# 
#select_cpdb_rows_cols(sig_mean=sig_mean_n,cells,fprefix=fprefix_n,mean_cutoff=1.5,out)
# 
# ## Recurrent ##
# 
sig_mean_r<-'/home/sonas/cellphonedb/cpdb_outputs100/Recurrent/significant_means.txt'
fprefix_r<-'recurrent'
# 
# select_cpdb_rows_cols(sig_mean=sig_mean_r,cells,fprefix=fprefix_r,mean_cutoff=1.5,out)



## Common cols and rows

out<-'/home/sonas/cellphonedb/cpdb_current/dotplots/'

# Common rows
comm_rc<-select_cpdb_common_rows_cols(sig_mean_r,sig_mean_n,cells,diff_cutoff=0.5,sub = uro_n)

select_comm_rows<-comm_rc[[1]]
select_comm_cols<-comm_rc[[2]]

write.table(select_comm_rows,paste0(out,'common_',cells2,'_rows.txt'),
            row.names = FALSE,col.names = FALSE,quote=FALSE)

write.table(select_comm_cols,paste0(out,'common_',cells2,'_cols.txt'),
            row.names = FALSE,col.names = FALSE,quote=FALSE)
