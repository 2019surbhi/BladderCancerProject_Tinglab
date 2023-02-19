## Author: Surbhi Sona ##

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)


##### INPUTS ##### 

dir<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/cellphonedb/group_comparison/group_comparison/'
prefix_n<-'Naive'
prefix_r<-'Recurrent'
file<-'count_network.txt'
ifile<-'significant_means.txt'       



###### FUNCTIONS #####

### Function to extract significant interactions ###

# dir: input dir 
# prefix: prefix that defines subdirectory
# ifile: interaction file name
# top: top interaction table with ipair and interaction diff as columns

get_sig_int<-function(dir,prefix,ifile,top)
{
  itab<-fread(paste0(dir,prefix,'/',ifile))
  reg_cols<-colnames(itab)[1:12]
  col<-colnames(itab)[-c(1:12)]
  col<-gsub('|','_',col,fixed = TRUE)
  colnames(itab)<-c(reg_cols,col)
  
  sub<-list()
  
  for(i in 1:nrow(top))
  {
    col<-match(top$ipairs[i],colnames(itab))
    if(is.na(col)==TRUE)
    {
      sub_tab<-as.data.frame(cbind('none',0))
      colnames(sub_tab)<-c('interacting_pair',top$ipairs[i])
    }else{
      selected_rows<-which(is.na(itab[,..col])==FALSE)
      selected_cols<-c(grep('interacting_pair',colnames(itab)),col)
      sub_tab<-itab[selected_rows,..selected_cols]
      sub_tab<-sub_tab[order(sub_tab[,2],decreasing = TRUE),]
    }  
    sub[[i]]<-sub_tab
  }
  
  names(sub)<-top$ipairs
  
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
  max<-apply(final_tab,1,max)
  rem<-which(max<=cutoff)
  tab<-tab[-rem,]
  return(tab)
  
}


### Function to prepare heatmap interactions ###

# tab: final merged interaction table 

prep_heatmap_input<-function(tab)
{
  mat<-as.matrix(tab)
  rownames(mat)<-rownames(tab)
  
  sorted_cols<-sort(colnames(mat))
  mat<-mat[,sorted_cols]
  
  cell_pair<-colnames(mat)
  
  # Group annotation
  Group<-vector()
  n<-grep('(N)',colnames(mat),fixed = TRUE)
  r<-grep('(R)',colnames(mat),fixed = TRUE)
  Group[n]<-rep('Naive',length(n))
  Group[r]<-rep('Recurrent',length(r))
  
  mycol<-c(brewer.pal(10,'Paired'),brewer.pal(10,'Set3'))
  cell_col<-rep(mycol,each=2)
  cell_col<-cell_col[1:length(cell_pair)]
  names(cell_col)<-cell_pair
  
  grp_col<-c('green','purple')
  names(grp_col)<-c('Naive','Recurrent')
  
  
  df<-data.frame(Group=Group,cell_pair=cell_pair)
  
  top<-HeatmapAnnotation(df=df,
                         col=list(Group=grp_col,
                                  cell_pair=cell_col),
                         show_annotation_name = TRUE)
  
  return(list(mat,top))
}

##### Read INPUTS ##### 

tn<-read.table(paste0(dir,prefix_n,'/',file),header=TRUE,sep='\t')
pn<-paste(tn$SOURCE,tn$TARGET,sep='_')
tn$pairs<-pn

tr<-read.table(paste0(dir,prefix_r,'/',file),header=TRUE,sep='\t')
pr<-paste(tr$SOURCE,tr$TARGET,sep='_')
tr$pairs<-pr

pairs<-c(pn,pr) %>% unique()


##### Calculate top interaction difference ##### 

idiff2<-vector()

for(i in 1:length(pairs))
{
  indx_n<-match(pairs[i],tn$pairs)
  indx_r<-match(pairs[i],tr$pairs)
  
  # ipair exclusive to Recurrent
  if(is.na(indx_n)==TRUE)
  {
    idiff2[i]<-tr$count[indx_r]
    idiff2[i]<-idiff2[i] *(-1)
  }else if(is.na(indx_r)==TRUE) # ipair exclusive to Naive
    {
     idiff2[i]<-tn$count[indx_n]
     }else{ 
            idiff2[i]<-tn$count[indx_n]-tr$count[indx_r]
           }
}

diff_tab2<-data.frame("ipairs"=pairs,"Naive-Recurrent"=idiff2)

# Remove rows exclusive to 1 group
#rem<-which(is.na(diff_tab2$Naive.Recurrent)==TRUE)
#diff_tab2<-diff_tab2[-rem,]

# Order
diff_tab2<-diff_tab2[order(diff_tab2$Naive.Recurrent,decreasing = TRUE),]

#top10<-diff_tab2[1:10,]

# Top20 stronger or exclusive to Naive
top20n<-diff_tab2[1:20,]


# Top20 stronger or exclusive to Recurrent

#diff_tab2$Naive.Recurrent<-abs(diff_tab2$Naive.Recurrent)
#diff_tab2<-diff_tab2[order(diff_tab2$Naive.Recurrent,decreasing = TRUE),]

n<-nrow(diff_tab2)
top20r<-diff_tab2[(n-19):n,]


##### Now read interactions #####

# (A) Stronger and exclusive to Naive
sub_n<-get_sig_int(dir,prefix_n,ifile,top20n)
sub_r<-get_sig_int(dir,prefix_r,ifile,top20n)


# (B) Stronger and exclusive to Recurrent
sub_n<-get_sig_int(dir,prefix_n,ifile,top20r)
sub_r<-get_sig_int(dir,prefix_r,ifile,top20r)


##### Common code #####

all_int_n<-sapply(sub_n,'[[',1) %>% unlist() %>% unique()
all_int_r<-sapply(sub_r,'[[',1) %>% unlist() %>% unique()
all_int<-c(all_int_n,all_int_r) %>% unique()

naive<-merge_sig_int(sub_n,all_int)
colnames(naive)<-paste0(colnames(naive),'(N)')
rec<-merge_sig_int(sub_r,all_int)
colnames(rec)<-paste0(colnames(rec),'(R)')


final_tab<-cbind(naive,rec) %>% as.data.frame()
rownames(final_tab)<-all_int


# (1) No split - filter rows # 

final_tab2<-filter_rows(final_tab,0.5)
lst<-prep_heatmap_input(final_tab2)
mat<-lst[[1]]
top<-lst[[2]]
fname<-'cpdb_naive_top20_filt0.5.png'
fname<-'cpdb_rec_clus_top20_filt0.5.png'

#fname<-'cpdb_naive_top20.png'
#fname<-'cpdb_rec_top20.png'

# (2) Split #

## Split into top immune and top non-immune cell ipairs ##

# needs to be done manually so inspect ipairs
colnames(final_tab)

i<-c(7,8,9,10,13,14,15,16,17,18,27,28,29,30,33,34,35,36,37,38)
ni<-c(1,2,3,4,5,6,11,12,19,20,21,22,23,24,25,26,31,32,39,40)

# 392 ligand-receptor interactions
final_tab_i<-final_tab[,i]
final_tab_ni<-final_tab[,ni]

# Filter rows 

# 171 ligand-receptor interactions
final_tab_i2<-filter_rows(final_tab_i,0.5)
final_tab_ni2<-filter_rows(final_tab_ni,0.5)

# (a) Immune
lst<-prep_heatmap_input(final_tab_i2)
mat<-lst[[1]]
top<-lst[[2]]
fname<-'cpdb_naive_top_filt0.5_imm.png'


# (a) Non-Immune

lst<-prep_heatmap_input(final_tab_ni2)
mat<-lst[[1]]
top<-lst[[2]]
fname<-'cpdb_naive_top_filt0.5_non_imm.png'


##### PLOTS #####

out<-'/home/sonas/cellphonedb/plots/'
#fname<-'cpdb_naive_top20.png'
#mycol<-c('#2C7BB6','#FFFFBF','#D7191C')
mycol<-c('white','red')
#mycol<-c('red','b','green')

hmcols<-colorRampPalette(mycol)(256)
h<-50
w<-15


png(filename = paste0(out,fname),width = w,height=h,units = "in",res=300)

Heatmap(mat, col=hmcols,column_split = Group,
        show_heatmap_legend = TRUE,
        top_annotation = top,name = 'mean_exp',
        column_labels = rep('',ncol(mat)),
        cluster_column_slices = FALSE,
        cluster_rows = TRUE, row_dend_side = 'right',
        cluster_columns = TRUE,row_names_side = 'left')
dev.off()


