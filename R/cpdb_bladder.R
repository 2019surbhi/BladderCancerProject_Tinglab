# The testing is done per group (total 3) based on BCG treatment status. 


# functions 

save_norm_exp_mat1<-function(obj,assay,fname)
{
  raw<-obj@assays[[assay]]@counts %>% as.data.frame.matrix()
  expr<-apply(raw,2,function(x) {x/sum(x) *10000})
  expr<-cbind(rownames(expr),expr)
  rownames(expr)<-NULL
  colnames(expr)<-c('Gene',colnames(expr)[-1])
  expr<-as.data.frame(expr)
  write.table(x = expr,
              file =fname, sep = '\t', 
              row.names = F, col.names = T, quote = F)
}

save_norm_exp_mat<-function(obj,assay,fname)
{
  dat<-obj@assays[[assay]]@data %>% as.data.frame.matrix()
  
  rownames(dat)<-NULL
  colnames(dat)<-c('Gene',colnames(dat)[-1])
  dat<-as.data.frame(dat)
  write.table(x = dat,
              file =fname, sep = '\t', 
              row.names = F, col.names = T, quote = F)
}

save_metadata<-function(obj,meta,fname)
{
 metadat<-as.data.frame(obj@meta.data[,meta])
 metadat<-cbind(rownames(obj@meta.data),metadat)
 rownames(metadat)<-NULL
 colnames(metadat)<-c('Cell',meta)
 
# Save meta data and table
 
 write.table(x = metadat,
             file =fname, sep = '\t', row.names = F, 
             col.names = T, quote = F)
}

metadata_cutoff_based_subset<-function(obj,meta='cell_label',cutoff=10)
{
tab<-table(obj@meta.data[,meta]) %>% as.data.frame()
k<-which(tab$Freq>cutoff)
keep<-as.character(tab$Var1[k])
Idents(obj)<-meta
sub_obj<-subset(obj,idents=keep)
return(sub_obj)
}


# read obj and subset by group

obj<-readRDS('dta_cancer_final.rds')
Idents(obj)<-'Groups'
sub_nw<-subset(obj,ident='Naive_w') # 32909 cells
sub_nwo<-subset(obj,ident='Naive_wo') # 103534 cells
sub_rec<-subset(obj,ident='Recurrent') # 113786 cells

# now subset by compartment
Idents(sub_nw)<-'low_res_compartment'
sub_nw_ui<-subset(sub_nw,ident=c('immune','uroepithelial')) # 3379 and 29491 cells respectively
                  
# check cell count
cell_tab<-table(sub_nw_ui@meta.data$cell.names) %>% as.data.frame()
colnames(cell_tab)<-c('cell_type','count')
write.table(cell_tab,'uro_imm_cell_count_Naive_w.txt',sep='\t',row.names=FALSE)

# Now create cpdb inputs

prefix<-'Naive_w'

sub<-metadata_cutoff_based_subset(obj=sub_nw,meta='cell.names',cutoff=10) # 30696 cells remain
### Save metadata ###
save_metadata(sub,meta='cell.names',fname=paste0(out,prefix,'_metadata.txt'))
### Save norm expression matrix ###
save_norm_exp_mat(sub,assay='RNACleaned',fname=paste0(out,prefix,'_counts.txt'))

                  
