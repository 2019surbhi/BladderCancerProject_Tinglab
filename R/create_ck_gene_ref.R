#!/usr/bin/env Rscript

library(dplyr)

fprefix<-'sample17_'

dir='/home/sonas/beegfs/results/scRNA/bladder_copyKat_heatmap/'
gene_tab_dir=paste0(dir,'gene_by_cell_tables/')

cat('Reading files \n')
f<-list.files(gene_tab_dir,full.names = TRUE)

files<-lapply(1:length(f), function(x){tab<-read.table(f[x],header = TRUE);
tab_selected<-tab[,(1:6)];
return(tab_selected)})

#saveRDS(files,paste0(dir,'gene_tab_list.rds'))

cat('Merging files \n')
merged_tab<-do.call(rbind,files)

class(merged_tab)

cat('Saving merged table as compressed file \n')
saveRDS(merged_tab,paste0(dir,'merged_gene_tab.rds'))

uniq<-unique(merged_tab$ensembl_gene_id)

uniq_merged_tab<-merged_tab[match(uniq,merged_tab$ensembl_gene_id),]

cat('Saving merged table \n')
write.table(uniq_merged_tab,paste0(dir,fprefix,'gene_ref.txt'),row.names=FALSE)
