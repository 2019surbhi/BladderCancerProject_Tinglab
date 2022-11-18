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

# Function to apply cell% cutoff

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

