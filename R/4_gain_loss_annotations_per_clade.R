source('/home/sonas/scripts/final_scripts/4_gain_loss_annotations_functions.R')

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

#path<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_r/annotated_table/genes_for_oncoplots.xlsx'


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


#tab_l<-'samples17_annotated_genomic_ranges_loss_70p_cutoff-0.03.xlsx'
#anno_df_rec<-read.xlsx(xlsxFile = paste0(path,tab_l))

get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_l,var='loss',p=pl,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,gene_tab = gene_tab17,fprefix=n17_fprefix,out=out)


## GAINS ##

cutoff_g<-(0.03)
pg<-70
var<-'gain'


#gtab_g<-'samples17_annotated_genomic_ranges_gain_70p_cutoff0.03.xlsx'
#anno_df_rec<-read.xlsx(xlsxFile = paste0(path,tab_g))


get_per_clade_annotated_table(thresh_dir=thresh_dir,clade=clade,cutoff=cutoff_g,var='gain',p=pg,chr_sizes=chr_sizes,rclade=rclade,nrclade=nrclade,gene_tab = gene_tab17,fprefix=n17_fprefix,out=out)




