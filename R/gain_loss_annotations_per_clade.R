# Author: Surbhi Sona

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
  
  #int_ts<-join_overlap_intersect(x=gr_db_ts,y=gr)
  #int_ts<-sort(int_ts)
  
  # Some genes span >1 interval and these duplicate entries need to be removed
  int<-int[unique(int$gene_id),]
  
  # Get annotations for gene symbol
  gene_id<-int$gene_id
  #ts_id<-as.character(int_ts$tx_id)
  
  anno<-AnnotationDbi::select(org.Hs.eg.db, keys=gene_id, columns=c('SYMBOL',"GENENAME"), keytype='ENTREZID')
  
  # keytypes(org.Hs.eg.db)
  # anno_ts<-AnnotationDbi::select(org.Hs.eg.db, keys=ts_id, columns='SYMBOL', keytype='ENSEMBLTRANS')
  
  
  int_df<-as.data.frame(int)
  int_anno_df<-merge(int_df,anno,by.x="gene_id",by.y='ENTREZID')
  
  return(int_anno_df)
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

get_per_clade_annotated_table<-function(thresh_dir,clade,cutoff,var,p,chr_sizes,anno_df_old='',fprefix,out)
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

    # remove row where start>end
    merged_all<-lapply(1:length(merged_all1),function(x){rem<-which((merged_all1[[x]]$end-merged_all1[[x]]$start)<0);
    if(length(rem)>0)
    {merged_all1[[x]]<-merged_all1[[x]][-rem,]};
    return(merged_all1[[x]])})

    names(merged_all)<-names(merged_all1)

    merged_tab<-full_join(x = merged_all[[1]],y = merged_all[[2]],by=c("chr","start","end","abspos_start","abspos_end"),suffix=c('',''))
    j<-2
    for(i in 1:(length(merged_all)-2))
    {
      j<-j+1
      merged_tab<-full_join(x=merged_tab,y = merged_all[[j]],by=c("chr","start","end","abspos_start","abspos_end"),suffix=c('',''))
      
    }

    ## Annotate ##

    # Remove rows where start>end

    # rem<-which((merged_tab$end-merged_tab$start)<0)
    # merged_tab<-merged_tab[-rem,]

    # annotate main table
    anno_df<-annotate_genomic_ranges(df = merged_tab,var='loss',
                                     cutoff=cutoff_l,plot=FALSE,
                                     bw=200,out=out,
                                     prefix=paste0(fprefix,'_loss_cutoff',
                                                   cutoff_l),
                                     save=FALSE,chr_size=chr_size)


    # annotate per clade tables

    anno_df_lst<-list()
    for(i in 1:length(merged_all))
    {
    anno_df_lst[[i]]<-annotate_genomic_ranges(df = merged_all[[i]],var='loss',
                                     cutoff=cutoff_l,plot=FALSE,
                                     bw=200,out=out,
                                     prefix=paste0(fprefix,'_loss_cutoff',
                                                   cutoff_l),
                                     save=FALSE,chr_size=chr_size)

    }

    names(anno_df_lst)<-names(merged_all)


    write.xlsx(anno_df_lst,paste0(out,fprefix,'_per_clade_annotated_',var,'_',p,'p_cutoff',cutoff,'.xlsx'))


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
    anno_df2<-cbind(anno_df,clade_df)


    ## Add R and NR count and % ##

    rclade<-paste0('clade',c(1,4,5,6,10))
    rcol<-match(rclade,colnames(anno_df2))
    R<-rowSums(anno_df2[,rcol])

    nrclade<-paste0('clade',c(3,7,8,9))
    nrcol<-match(nrclade,colnames(anno_df2))
    NR<-rowSums(anno_df2[,nrcol])

    pct.R<-(R/5) *100
    pct.NR<-(NR/4)*100

    count_df<-cbind(R,NR,pct.R,pct.NR) %>% as.data.frame()

    #count_df$pct.R<-paste0(count_df$pct.R,'%')
    #count_df$pct.NR<-paste0(count_df$pct.NR,'%')

    anno_df2<-cbind(anno_df2,count_df)


    ## Add overlap column - overlap with old table ##
    
    
    if(anno_df_old!='')
    {
    # Intersect to get common genes #
    comm<-intersect(anno_df_old$SYMBOL,anno_df2$SYMBOL)

    anno_df2$overlap_w_old<-rep(0,nrow(anno_df2))
    idx<-which(anno_df2$SYMBOL %in% comm)
    anno_df2$overlap_w_old[idx]<-1
    }

    write.xlsx(anno_df2,paste0(out,fprefix,'_per_clade_merged_',var,'_',p,'p_cutoff',cutoff,'_final.xlsx'))


}



chr_sizes<-fread('//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/Misc/hg38.chrom.sizes2.txt')



### N=8 ####
out<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/annotated_table/'

thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/primary_recurrent_per_clade/thresholds/'

n8_fprefix<-'primary_recurrent'
clade<-paste0('clade',1:10)

path<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_r/annotated_table/genes_for_oncoplots.xlsx'


## LOSS ##

#loss
cutoff_l<-(-0.03)
pl<-70
var<-'loss'

anno_df_rec<-read.xlsx(xlsxFile = path,sheet = var)


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_l,var='loss',p=pl,chr_sizes=chr_sizes,anno_df_old=anno_df_rec,fprefix=n8_fprefix,out=out)


## GAINS ##

cutoff_g<-(0.03)
pg<-70
var<-'gain'

anno_df_rec<-read.xlsx(xlsxFile = path,sheet = var)


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_g,var='gain',p=pg,chr_sizes=chr_sizes,anno_df_old=anno_df_rec,fprefix=n8_fprefix,out=out)



### N=17 ###

out<-'/home/sonas/copykat/genes/gain_loss/sample17_per_clade/annotated_table/'

thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17_per_clade/thresholds/'

n17_fprefix<-'samples17'
clade<-paste0('clade',1:17)

## LOSS ##
#cutoff_l<-(-0.035)
cutoff_l<-(-0.030)
pl<-70
var<-'loss'

get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_l,var='loss',p=pl,chr_sizes=chr_sizes,fprefix=n17_fprefix,out=out)


## GAINS ##

cutoff_g<-(0.03)
pg<-70
var<-'gain'


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_g,var='gain',p=pg,chr_sizes=chr_sizes,fprefix=n17_fprefix,out=out)




### UPSET plot ###

loss<-read.xlsx(paste0(out,fprefix,'_all_coordinates_',pl,'p_across_clades_cutoff<',cutoff_l,'_final.xlsx'))
gain<-read.xlsx(paste0(out,fprefix,'_all_coordinates_',pg,'p_across_clades_cutoff>',cutoff_g,'_final.xlsx'))

var<-'loss'
tab<-loss

var<-'gain'
tab<-gain

rem<-which(is.na(tab$SYMBOL)==TRUE)
tab<-tab[-rem,]

clades<-paste0('clade',c(1,4,5,6,10,3,7,8,9))

cl_lst<-list()
for(i in 1:length(clades))
{
idx<-which(tab[,match(clades[i],colnames(tab))]==1)
cl_lst[[i]]<-tab$SYMBOL[idx]
}

names(cl_lst)<-clades

out<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/'
fprefix<-'primary_recurrent_per_clade_overlap'

odr<-names(cl_lst)
m<-'distinct'
#m<-'intersect'
w<-25
h<-7

cmat<-make_comb_mat(cl_lst, mode =m)
cs<-comb_size(cmat)

fname<-paste(fprefix,var,m,sep='_')
ht=UpSet(cmat, set_order = odr,
         comb_order = order(comb_size(cmat),decreasing = TRUE),
         top_annotation = HeatmapAnnotation(
           "Intersection size" = anno_barplot(cs,
                                              ylim = c(0, max(cs)*1.1),
                                              border = FALSE,
                                              gp = gpar(fill = "black"),
                                              height = unit(4, "cm")),
           annotation_name_side = "left",
           annotation_name_rot = 90),
         right_annotation = upset_right_annotation(cmat, add_numbers = TRUE))

pdf(paste0(out,fname,'_upset_plot.pdf'),width=w,height = h)
ht=draw(ht)
od=column_order(ht)
ht=ht+decorate_annotation("Intersection size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "mm"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 18, col = "#404040"))})
dev.off()


########################################################


# Intersect

# min 1 clade in R group

#rec<-paste0('clade',c(4,5,6,10))
#nrec<-paste0('clade',c(3,7,8,9))
#2,3,5,7,15

# sample17

#sum of reccurent except clade1

# Remove  blank rows

rem<-which(is.na(anno_df2$SYMBOL))
anno_df2<-anno_df2[-rem,]

anno_df2$rec<-rowSums(anno_df2[,c(15,16,17,21)])
anno_df2$nrec<-rowSums(anno_df2[,c(14,18,19,20)])

c1<-which(anno_df2$clade1>=1)
c1_genes<-anno_df2$SYMBOL[c1]

r<-which(anno_df2$rec>=1)
nr<-which(anno_df2$nrec>=1)

rec_genes<-anno_df2$SYMBOL[r]
nrec_genes<-anno_df2$SYMBOL[nr]

common<-intersect(rec_genes,nrec_genes)

rec_only<-setdiff(rec_genes,nrec_genes)
nrec_only<-setdiff(nrec_genes,rec_genes)

c1_int_rec<-intersect(c1_genes,rec_only)
c1_int_nrec<-intersect(c1_genes,nrec_only)
c1_int_comm<-intersect(c1_genes,common)

length(c1_int_rec)
length(c1_int_nrec)
length(c1_int_comm)


length(rec_genes)
length(nrec_genes)
length(common)

length(rec_only)
length(nrec_only)

length(c1)


### Upset Plot ###

out<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/'

fprefix<-'primary_recurrent'
var<-'loss'
var<-'gain'

dat<-list(c1_genes,common,rec_only,nrec_only)
names(dat)<-c('clade1','common','R_only','NR_only')
odr<-names(dat)
m<-'distinct'
#m<-'intersect'
w<-18
h<-7

cmat<-make_comb_mat(dat, mode =m)
cs<-comb_size(cmat)

fname<-paste(fprefix,var,m,sep='_')
ht=UpSet(cmat, set_order = odr,
         comb_order = order(comb_size(cmat),decreasing = TRUE),
         top_annotation = HeatmapAnnotation(
           "Intersection size" = anno_barplot(cs,
                                              ylim = c(0, max(cs)*1.1),
                                              border = FALSE,
                                              gp = gpar(fill = "black"),
                                              height = unit(4, "cm")
           ),
           annotation_name_side = "left",
           annotation_name_rot = 90,gp = gpar(fontsize=18)))

pdf(paste0(out,fname,'_upset_plot.pdf'),width=w,height = h)
ht=draw(ht)
od=column_order(ht)
ht=ht+decorate_annotation("Intersection size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "mm"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 20, col = "#404040"))})
dev.off()




