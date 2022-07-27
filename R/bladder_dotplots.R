library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(readxl)
library(gridExtra)

out<-'/home/sonas/copykat/plots/'
prefix<-'primary_recurrent_k10'


### 0. Functions ###

# DotPlots per sample

dotplot_per_sample<-function(obj,genes,out,fname,w=20,h=12,meta)
{
  obj_lst<-SplitObject(obj,split.by = 'orig.ident')
  sname<-names(obj_lst)
  
  fname2<-paste(fname,sname,sep='_')
  plist<-list()
  for(i in 1:length(obj_lst))
  {
    DefaultAssay(obj_lst[[i]])<-'RNACleaned'
    obj_lst[[i]]<-add_cell_count_to_cell_names(obj_lst[[i]],meta = meta)
    Idents(obj_lst[[i]])<-meta
    
    #obj_lst[[i]]@active.ident<-factor(x=obj_lst[[i]]@active.ident, levels=new_levels)
    
    d<-DotPlot(obj_lst[[i]],features=genes) + theme(axis.text.x = element_text(angle = 90))
    plist[[i]]<-d
    ggsave(filename = paste0(out,fname2[i],'_Dotplot.png'),width =w,height=h,units="in",d)
    
  }
  
  m<-marrangeGrob(grobs = plist,ncol = 1,nrow = 1)
  
  ggsave(paste0(out,fname,'_Dotplots.pdf'),width=w,height = h,units="in",m)
  
}

# Function to add cell count to cell_names (in future modify this for any meta col)

add_cell_count_to_cell_names<-function(obj,meta='seurat_clusters')
{
  cells<-table(obj@meta.data[,meta]) %>% as.data.frame()
  new_meta<-rep('none',nrow(obj@meta.data))
  for(i in 1:nrow(cells))
  {
    idx<-which(obj@meta.data[,meta] %in% cells$Var1[i])
    new<-paste0(obj@meta.data[idx,meta],' (',cells$Freq[i],')')
    new_meta[idx]<-new
  }
  
  obj<-AddMetaData(obj,metadata = new_meta,col.name =meta)
  
  return(obj)
}
### 1. Read objects ###

# Immune
imm_all<-readRDS('/home/sonas/data/immune.v11.harmony_umap.pass7p.v2.cn1.rds')

# stromal

endo_all<-readRDS('/home/sonas/data/v11.pass7p.endothelial.clean.v2.rds')
fibro_all<-readRDS('/home/sonas/data/v11.pass7p.fibroblasts.clean.v2.rds')

# Uro
pr_obj<-readRDS('/home/sonas/data/v11.pass7p.uroepithelial.clean.v2_hclust_primary_recurrent.rds')


### 2. obj Processing ###

## IMMUNE - Subset to include only the 8 samples ## 
s<-unique(pr_obj$orig.ident)
Idents(imm_all)<-'orig.ident'
imm<-subset(imm_all,idents=s)

## STROMAL 

Idents(fibro_all)<-'orig.ident'
#remove missing sample
rem<-which((s %in% levels(fibro_all))==FALSE)
s2<-s[-rem]
fibro<-subset(fibro_all,idents=s2)


Idents(endo_all)<-'orig.ident'
#one of the samples is missing
rem<-which((s %in% levels(endo_all))==FALSE)
s2<-s[-rem]
endo<-subset(endo_all,idents=s2)

  
## URO - add recurrent/non-recurrent annotation ##
rec_clus<-c(1,4,5,6,10)
non_rec_clus<-c(2,3,7,8,9)

r<-vector("integer")
nr<-vector("integer")

for(i in 1:length(rec_clus))
{
  idx<-which(pr_obj@meta.data$hclusters==rec_clus[i])
  r<-c(r,idx)
}

for(i in 1:length(non_rec_clus))
{
  idx2<-which(pr_obj@meta.data$hclusters==non_rec_clus[i])
  nr<-c(nr,idx2)
}

meta<-rep('clusters',nrow(pr_obj@meta.data))
meta[r]<-'recurrent_associated'
meta[nr]<-'non_recurrent_associated'

pr_obj<-AddMetaData(pr_obj,metadata = meta,col.name = 'recurrent')



### 3. Read in genes ###

tab_cpi<-read_xlsx('/home/sonas/copykat/misc/CHECKPOINT GENES_updated_7_27_2022.xlsx',sheet = 'cpi')
tab_ck<-read_xlsx('/home/sonas/copykat/misc/CHECKPOINT GENES_updated_7_27_2022.xlsx',sheet = 'cytokines')
tab_cp<-read_xlsx('/home/sonas/copykat/misc/cellphonedb_gene_list.xlsx')

genes1<-tab_cpi$gene %>% unique()
genes2<-tab_ck$gene %>% unique()
genes3<-c((tab_cp$gene1 %>% unique()),tab_cp$gene2 %>% unique())


### 4. Prepare variables for plotting ###

# IMMUNE ##
DefaultAssay(imm)<-'RNACleaned'
imm<-add_cell_count_to_cell_names(imm,meta = 'cell_names')
Idents(imm)<-'cell_names2'
prefix_imm<-paste0(prefix,'_immune')


## STROMAL ##


## URO ##
DefaultAssay(pr_obj)<-'RNACleaned'
pr_obj<-add_cell_count_to_cell_names(pr_obj,meta='recurrent')
Idents(pr_obj)<-'recurrent'
prefix_uro<-paste0(prefix,'_uro')


### 5. Dotplots ###

## (A) Checkpoint inhibitors ##

out_cpi<-paste0(out,'checkpoint_inhibitors_plots/')

# Immune
di<-DotPlot(imm,features=genes1) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_cpi,'immune/',
              prefix_imm,'_chkpoint_inhibitor_Dotplot.png'), width = 20,height=12,units="in",di)

dotplot_per_sample(obj=imm,genes=genes1,out=paste0(out_cpi,'immune/'),
                   fname=paste0(prefix_imm,'_chkpoint_inhibitors'),meta = 'cell_names')
  

# Stromal

# Uro
du<-DotPlot(pr_obj,features=genes1) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_cpi,'uro/',prefix_uro,
              '_chkpoint_inhibitor_Dotplot.png'),width = 20,height=12,units="in",du)

dotplot_per_sample(obj=pr_obj,genes=genes1,out=paste0(out_cpi,'uro/'),
                   fname=paste0(prefix_uro,'_chkpoint_inhibitors'),meta = 'recurrent')


## (B) Cytokines/Chemokines ##

out_ck<-paste0(out,'cytokines_plots/')

# Immune
di<-DotPlot(imm,features=genes2) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_ck,'immune/',prefix_imm,
              '_cytokines_Dotplot.png'),width = 20,height=12,units="in",di)

dotplot_per_sample(obj=imm,genes=genes2,out=paste0(out_ck,'immune/'),
                   fname=paste0(prefix_imm,'_cytokines'),meta = 'cell_names')

# stromal

# uro
du<-DotPlot(pr_obj,features=genes2) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_ck,'uro/',prefix_uro,
              '_cytokines_Dotplot.png'),width = 20,height=12,units="in",du)

dotplot_per_sample(obj=pr_obj,genes=genes2,out=paste0(out_ck,'uro/'),
                   fname=paste0(prefix_uro,'_cytokines'),meta = 'recurrent')

## (C) Cellphonedb genes ## - incomplete


# immune
di<-DotPlot(imm,features=genes3) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_imm,prefix_imm,
              '_cpdb_Dotplot.png'),width = 20,height=12,units="in",di)

# stromal


# uro
du<-DotPlot(pr_obj,features=genes3) + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(out_uro,prefix_uro,
              '_cpdb_Dotplot.png'),width = 20,height=12,units="in",du)


