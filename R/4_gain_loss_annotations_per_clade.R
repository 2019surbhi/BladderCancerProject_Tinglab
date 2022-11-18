# Author: Surbhi Sona

source('4_gain_loss_annotations_functions.R')

chr_sizes<-fread('//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/Misc/hg38.chrom.sizes2.txt')



### N=8 ####

gene_tab<-read.table('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/Misc/primary_recurrent_gene_ref.txt',sep='\t',header=TRUE)

out<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/annotated_table/'

#thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/primary_recurrent_per_clade/thresholds/'

# New threhold directory which now has cell count and cell% 
thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/primary_recurrent_per_clade/new_threshold/'

n8_fprefix<-'primary_recurrent'
clade<-paste0('clade',1:10)

rclade<-paste0('clade',c(1,4,5,6,10))
nrclade<-paste0('clade',c(3,7,8,9))

path<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_r/annotated_table/genes_for_oncoplots.xlsx'


## LOSS ##

#loss
cutoff_l<-(-0.03)
pl<-70
var<-'loss'


#anno_df_rec<-read.xlsx(xlsxFile = path,sheet = var)

get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_l,var='loss',p=pl,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,fprefix=n8_fprefix,out=out,gene_tab = gene_tab)


## GAINS ##

cutoff_g<-(0.03)
pg<-70
var<-'gain'

#anno_df_rec<-read.xlsx(xlsxFile = path,sheet = var)


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_g,var='gain',p=pg,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,fprefix=n8_fprefix,out=out,gene_tab = gene_tab)




### N=17 ###

gene_tab17<-read.table('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/Misc/sample17_gene_ref.txt',sep='\t',header=TRUE)

out<-'/home/sonas/copykat/genes/gain_loss/sample17_per_clade/annotated_table/'

#thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17_per_clade/thresholds/'

thresh_dir<-'//home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17_per_clade/new_threshold/'

path<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17/annotated_table/'

n17_fprefix<-'samples17'
clade<-paste0('clade',1:17)

rclade<-paste0('clade',c(1,4,6,8,9,10,11,12,13,14,16,17))
nrclade<-paste0('clade',c(2,5,7,15))

## LOSS ##
#cutoff_l<-(-0.035)
cutoff_l<-(-0.030)
pl<-70
var<-'loss'


tab_l<-'samples17_annotated_genomic_ranges_loss_70p_cutoff-0.03.xlsx'
anno_df_rec<-read.xlsx(xlsxFile = paste0(path,tab_l))

get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_l,var='loss',p=pl,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,anno_df_old=anno_df_rec,fprefix=n17_fprefix,out=out)


## GAINS ##

cutoff_g<-(0.03)
pg<-70
var<-'gain'


gtab_g<-'samples17_annotated_genomic_ranges_gain_70p_cutoff0.03.xlsx'
anno_df_rec<-read.xlsx(xlsxFile = paste0(path,tab_g))


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_g,var='gain',p=pg,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,anno_df_old=anno_df_rec,fprefix=n17_fprefix,out=out)




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



