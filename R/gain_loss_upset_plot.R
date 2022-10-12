# Author: Surbhi Sona (Ting Lab) #


print_upset_plot<-function(dat,out,fprefix,var,w=18,h=7,m='distinct')
{
odr<-names(dat)
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

}





out<-'/home/sonas/copykat/genes/gain_loss/upset_plot/'
fprefix<-'gene_overlap'

### LOSS ###

var<-'loss'

## N17 ##

n17_old_dir<-'/home/sonas/copykat/genes/gain_loss/samples17/annotated_table/'
n17_per_clade_dir<-'/home/sonas/copykat/genes/gain_loss/sample17_per_clade/annotated_table/'

# Read old merged table #

n17_old_l<-read.xlsx(paste0(n17_old_dir, 'samples17_annotated_genomic_ranges_loss_70p_cutoff-0.03.xlsx'))

# Read per clade
n17_pc_l<-read.xlsx(paste0(n17_per_clade_dir,'samples17_per_clade_merged_loss_70p_cutoff-0.03_final.xlsx'))

rem<-which(is.na(n17_pc_l$SYMBOL)==TRUE)
n17_pc_l<-n17_pc_l[-rem,]


## N8 ##

n8_old_file<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_r/annotated_table/genes_for_oncoplots.xlsx'

n8_per_clade_dir<-'/home/sonas/copykat/genes/gain_loss/primary_recurrent_per_clade/annotated_table/'

# Read old merged table #

n8_old_l<-read.xlsx(xlsxFile = n8_old_file,sheet = var)

# Read per clade
n8_pc_l<-read.xlsx(paste0(n8_per_clade_dir,'primary_recurrent_per_clade_merged_loss_70p_cutoff-0.03_final.xlsx'))

rem<-which(is.na(n8_pc_l$SYMBOL)==TRUE)
n8_pc_l<-n8_pc_l[-rem,]


# dat<-list(n8_old_l$SYMBOL,n8_pc_l$SYMBOL)
# names(dat)<-c('n8_old','n8_per_clade')
# names(dat)<-paste0(names(dat),'_l')

#
dat<-list(n8_old_l$SYMBOL,n8_pc_l$SYMBOL,n17_old_l$SYMBOL,n17_pc_l$SYMBOL)
names(dat)<-c('n8_old','n8_per_clade','n17_old','n17_per_clade')
names(dat)<-paste0(names(dat),'_l')

print_upset_plot(dat=dat,out,fprefix=fprefix,var=var,w=18,h=7,m='distinct')




### GAIN ###

# Read old merged table #

var<-'gain'

n17_old_g<-read.xlsx(paste0(n17_old_dir,'samples17_annotated_genomic_ranges_gain_70p_cutoff0.03.xlsx'))
    
# Read per clade
n17_pc_g<-read.xlsx(paste0(n17_per_clade_dir,'samples17_per_clade_merged_gain_70p_cutoff0.03_final.xlsx'))

rem<-which(is.na(n17_pc_g$SYMBOL)==TRUE)
n17_pc_g<-n17_pc_g[-rem,]

## N8 ##

n8_old_g<-read.xlsx(xlsxFile = n8_old_file,sheet = var)

# Read per clade
n8_pc_g<-read.xlsx(paste0(n8_per_clade_dir,'primary_recurrent_per_clade_merged_gain_70p_cutoff0.03_final.xlsx'))

rem<-which(is.na(n8_pc_g$SYMBOL)==TRUE)
n8_pc_g<-n8_pc_g[-rem,]

dat<-list(n8_old_g$SYMBOL,n8_pc_g$SYMBOL,n17_old_g$SYMBOL,n17_pc_g$SYMBOL)
names(dat)<-c('n8_old','n8_per_clade','n17_old','n17_per_clade')
names(dat)<-paste0(names(dat),'_g')

print_upset_plot(dat=dat,out,fprefix=fprefix,var=var,w=18,h=7,m='distinct')


### Upset plot ###

# Only single set #


dat8<-dat[1:2]
dat17<-dat[3:4]

fprefix<-'gene_overlap_n8'
print_upset_plot(dat=dat8,out,fprefix=fprefix,var=var,w=18,h=7,m='distinct')
#print_upset_plot(dat=dat8,out,fprefix=fprefix,var=var,w=18,h=7,m='intersect')

fprefix<-'gene_overlap_n17'
print_upset_plot(dat=dat17,out,fprefix=fprefix,var=var,w=18,h=7,m='distinct')
#print_upset_plot(dat=dat17,out,fprefix=fprefix,var=var,w=18,h=7,m='intersect')


