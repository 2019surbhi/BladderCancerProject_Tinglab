library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)

##### INPUTS ##### 

dir<-'/home/sonas/cellphonedb/cpdb_outputs100/'
prefix_n<-'Naive'
prefix_r<-'Recurrent'
file<-'count_network.txt'
ifile<-'significant_means.txt'       



###### FUNCTIONS #####

### Function to extract significant interactions ###

# dir: input dir 
# prefix: prefix that defines subdirectory
# ifile: interaction file name
# cell_pairs: vector of cell pairs to subset for (name and format must match those in the significant mean table)


get_sig_int2<-function(dir,prefix,ifile,cell_pairs)
{
  itab<-fread(paste0(dir,prefix,'/',ifile))
  reg_cols<-colnames(itab)[1:12]
  col<-colnames(itab)[-c(1:12)]
  col<-gsub('|','_',col,fixed = TRUE)
  colnames(itab)<-c(reg_cols,col)
  
  sub<-list()
  
  for(i in 1:length(cell_pairs))
  {
    col<-match(cell_pairs[i],colnames(itab))
    if(is.na(col)==TRUE)
    {
      sub_tab<-as.data.frame(cbind('none',0))
      colnames(sub_tab)<-c('interacting_pair',cell_pairs[i])
    }else{
      selected_rows<-which(is.na(itab[,..col])==FALSE)
      selected_cols<-c(grep('interacting_pair',colnames(itab)),col)
      sub_tab<-itab[selected_rows,..selected_cols]
      sub_tab<-sub_tab[order(sub_tab[,2],decreasing = TRUE),]
    }  
    sub[[i]]<-sub_tab
  }
  
  names(sub)<-new_pairs
  
  return(sub)
}

### Function to merge significant interactions ###

#sub: list of significant tables (1 per ipair) 
#all_int: non-redundant list of ligand-receptor interactions

merge_sig_int<-function(sub,all_int)
{
  
  final_int<-list()
  
  for(i in 1:length(sub))
  {
    int<-sub[[i]][(match(all_int,(sub[[i]]$interacting_pair))),2]
    if(length(int)==length(is.na(int)==TRUE))
    {
      int<-rep(0,length(int)) %>% as.data.frame()
      colnames(int)<-names(sub[i])
      
    }else{
      int<-int %>% replace(is.na(.),0)
      int<-round(int,digits = 2)
    }
    
    final_int[[i]]<-int
  }
  
  #names(final_int)<-names(sub)
  final_int<-do.call(cbind,final_int)
  
  return(final_int)
}

### Function to filter rows ###

# tab: table of significant interactions
# cutoff: combined mean value cutoff for filtering

filter_rows<-function(tab,cutoff=0.5)
{
  max<-apply(tab,1,max)
  rem<-which(max<=cutoff)
  tab2<-tab[-rem,]
  return(tab2)
  
}

filter_rows_rev<-function(tab,cutoff=0.5)
{
  max<-apply(tab,1,max)
  rem<-which(max>=cutoff)
  tab2<-tab[-rem,]
  return(tab2)
  
}


### Function to select cell pairs ###

select_cell_pairs<-function(cell_group)
{

# Rename cell pairs
select_pairs<-vector()

for(i in 1:length(cell_group))
{
  cells<-paste0(cell_group[i],'_',cell_group)
  select_pairs<-c(select_pairs,cells)
}

return(select_pairs)

}

### Function to create barplots ###

# cell_pairs: vector of cell pairs in correct format: cell1_cell2 instead of cell1|cell2
# n: number of top and bottom interactions to subset for (e.g. n=20 for top20 and bottom 20 interactions)
# out: output dir path
# fprefix: file name prefix (e.g. uro-uro)
# w: plot width in inches (default is 8)
# h: plot height in inches (default is 11)

count_diff_barplot<-function(cell_pairs,n='all',out,fprefix,w=8,h=11)
{
    # Get significant means for selected cell pairs in Naive and Recurrent groups
  s
  sub_n<-get_sig_int2(dir,prefix_n,ifile,cell_pairs)
  sub_r<-get_sig_int2(dir,prefix_r,ifile,cell_pairs)
  
  # Get non-redundant union list of ligands and receptors
  
  all_int_n<-sapply(sub_n,'[[',1) %>% unlist() %>% unique()
  all_int_r<-sapply(sub_r,'[[',1) %>% unlist() %>% unique()
  all_int<-c(all_int_n,all_int_r) %>% unique()
  
  # Get merged significant mean table
  naive<-merge_sig_int(sub_n,all_int) %>% as.data.frame()
  rownames(naive)<-all_int
  rec<-merge_sig_int(sub_r,all_int) %>% as.data.frame()
  rownames(rec)<-all_int
  
 
  # Get interaction count
  ncount<-rowSums(naive!=0)
  rcount<-rowSums(rec!=0)
  
  #sanity check
  identical(names(ncount),names(rcount))
  
  # Get count difference
  diff<-rcount-ncount
  names(diff)<-names(ncount)
  
  # Sort
  diff<-sort(diff,decreasing = TRUE)
  
  # extract top and bottom diff (optional)
  if(n!='all')
  {
    t1<-diff[1:n]
    t2<-diff[(length(diff)-n):(length(diff))]
    diff<-c(t1,t2)
  }
  
  diff_tab<-data.frame('LR'=names(diff),'count_diff'=diff)
  diff_tab$value<-ifelse(diff_tab$count_diff>0,'R-high','N-high')
  
  # Remove diff count that is 0
  rem<-which(diff_tab$count_diff==0)
  if(length(rem)!=0)
  {
  diff_tab<-diff_tab[-rem,]
  }
  # print barplot
  
  b<-ggplot(diff_tab,aes_string(x='reorder(LR,-count_diff)',
                                y='count_diff',
                                fill='value'),
            position='stat',
            colour='black',width = 0.9) +
    geom_col() +
    scale_fill_manual(
      values = c("R-high" = "red",
                 "N-high" = "blue")) +
    xlab('Ligand-Receptor') +
    ylab('Count diff (R-N)') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line('black'),
          axis.text.x = element_text(angle=90)) +
    coord_flip()
  
  ggsave(paste0(out,fprefix,'_R-N_count_diff_barplot.png'),width = w,height = h,units = "in",dpi = 600,limitsize = FALSE,b)
  
}


##### Read INPUTS ##### 

tn<-fread(paste0(dir,prefix_n,'/',file),header=TRUE,sep='\t')
pn<-paste(tn$SOURCE,tn$TARGET,sep='_')
tn$pairs<-pn

tr<-fread(paste0(dir,prefix_r,'/',file),header=TRUE,sep='\t')
pr<-paste(tr$SOURCE,tr$TARGET,sep='_')
tr$pairs<-pr

pairs<-c(pn,pr) %>% unique()


##### Now read interactions #####

# select cell pair groups for uro-uro, imm-imm and imm-uro interactions

all_cells<-c(tn$SOURCE,tn$TARGET,tr$SOURCE,tr$TARGET) %>% unique()

uro<-all_cells[c(1:6,25)] %>% sort()
imm<-all_cells[c(7:23,26:28)] %>% sort()


######## Print barplots ##########

out<-'/home/sonas/cellphonedb/cpdb_outputs100/plots/'

# all pairs #

all_pairs<-select_cell_pairs(c(imm,uro))

fprefix<-'all_imm_uro_top20'
count_diff_barplot(cell_pairs=all_pairs,n=20,out,fprefix,w=8,h=11)

fprefix<-'all_imm_uro_all_sig'
count_diff_barplot(cell_pairs =all_pairs,n='all',out,fprefix,w=8,h=55)

# uro-uro interactions #

uro_pairs<-select_cell_pairs(uro)

fprefix<-'uro-uro_top20'
count_diff_barplot(cell_pairs=uro_pairs,n=20,out,fprefix,w=8,h=11)

fprefix<-'uro-uro_all_sig'
count_diff_barplot(cell_pairs =uro_pairs,n='all',out,fprefix,w=8,h=22)

# imm-imm interactions #

imm_pairs<-select_cell_pairs(imm)

fprefix<-'imm-imm_top20'
count_diff_barplot(cell_pairs=imm_pairs,n=20,out,fprefix,w=8,h=11)

fprefix<-'imm-imm_all_sig'
count_diff_barplot(cell_pairs=imm_pairs,n = 'all',out,fprefix,w=8,h=41)


# imm-uro interactions #

imm_uro_pairs<-vector()

for(i in 1:length(imm))
  {
    cells<-paste0(imm[i],'_',uro)
    imm_uro_pairs<-c(imm_uro_pairs,cells)
  }

fprefix<-'imm-uro_top20'
count_diff_barplot(cell_pairs=imm_uro_pairs,n=20,out,fprefix,w=8,h=11)

fprefix<-'imm-uro_all_sig'
count_diff_barplot(cell_pairs =imm_uro_pairs,n='all',out,fprefix,w=8,h=55)
