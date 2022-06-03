# CNV based analysis using copyKat outputs for scRNA data #
# Author: Surbhi Sona (Ting Lab) #


## Load Dependencies ##

library('TxDb.Hsapiens.UCSC.hg38.knownGene') # hg38 gene annotation ref
library('AnnotationDbi')
library('org.Hs.eg.db') # to get gene symbol for Entrez gene.id
library('plyrange') # to apply dplyr functions on GRange obj


### PRIMARY-RECURRENT pairs of sampels (n=8)

## Create subset of genomic intervals for gain/loss based on CNV value cutoff ##

# Step1: First subset CNV output matrix to include only cells from recurrence-associated clades
# Step2: Create genomic interval subset for this subset matrix based on CNV value cutoff
# Step3: Annoatate subset intervals to determine the genes that fall within these intervals


## Step 1a: Create abspos size in chr size table

# Read chr size file for hg38
chr_size<-fread(chr_size_file)
colnames(chr_size)<-c('chr','size')
chr_size<-as.data.frame(chr_size)

# Only reatin the regular chromosomes
chr_size<-chr_size[1:24,]

#Since copyKat doesn't have chrY we will remove it
rem<-which(chr_size$chr=='chrY')
chr_size<-chr_size[-rem,]

# Now rename chrX as chr23
ren<-which(chr_size$chr=='chrX')
chr_size$chr[ren]<-'chr23'

# Sort by chr
odr<-mixedsort(chr_size$chr)
chr_size<-chr_size[match(odr,chr_size$chr),]

# Add abspos size
pos<-chr_size$size[1] %>% as.double()
add<-pos %>% as.double()

n<-nrow(chr_size)
for(i in 2:n)
{
  add<-add+chr_size$size[i]
  pos<-c(pos,add)
}
chr_size$abs_size<-pos

write.table(chr_size,'/home/sonas/copykat/misc/hg38.chrom.sizes2.txt')


## Step 1b: Create Master genomic interval ref ##

'''Each copyKat run (depending on chosen resolution) has a defined set of genomic interval comprising of chrom, chrompos (start position) and abspos (absolute start position). I am simply adding the end position for chromosome start and absolute start positions. For any given start position (pos), the corresponding end is  calculated as end[pos]=start[pos+1]-1. If pos is the last fragment then end=chr/abspos size for the correspoinding chr '''

# Read subset matrix (recurrence-associated) if needed
rec_mat2<-readRDS(rec_mat2,paste0(dir,'matrix/',sname,'_k10_recurrent_subset_matrix.rds'))
df<-as.data.frame(rec_mat2[,1:3])

# If hg38 then add 'chr' to chr number
df$chrom<-paste0('chr',df$chr)
  
# Sort by abspos
df<-df[order(df$chrom),]
  
# Calculate the ends per chr
chrs<-unique(df$chr)

df2<-NULL
for(i in 1:length(chrs))
 {
   df_temp<-add_ends(df_by_chr=df[df$chr==chrs[i],],chr_size = chr_size)
   if(is.null(df2)==TRUE)
   {
     df2<-df_temp
   }else
   {
   df2<-rbind(df2,df_temp)
   }
}

## Now add ends to abspos ##
df3<-NULL
for(i in 1:length(chrs))
 {
   df_temp<-add_abspos_ends(df_by_chr=df2[df2$chr==chrs[i],],chr_size = chr_size)
   if(is.null(df3)==TRUE)
   {
     df3<-df_temp
   }else
   {
   df3<-rbind(df3,df_temp)
   }
}

write.table(df3,'/home/sonas/copykat/misc/copyKat_220kb_genomic_fragments_table.txt',sep='\t')

## Step1: Subset matrix ##

dir<-'/home/sonas/copykat/'
sname<-'recurrent_pairs'
mat<-readRDS(paste0(dir,'matrix/',sname,'_matrix.rds'))
chr_size_file<-'/home/sonas/copykat/misc/hg38.chrom.sizes'


# Read hclus obj
hcc<-readRDS(paste0(dir,'hclust/',sname,'_Hclust.rds'))

# Extract barcodes of cells in recurrence-associated clades
hc.umap<-cutree(hcc,k=10)
rclades<-c(1,4,5,6,10)
idx<-which( hc.umap %in% rclades)
bc<-names(hc.umap)[idx]

# Subset matrix to include only recurrence-associated clades cells
cell_idx<-which(colnames(mat) %in% bc)
rec_mat<-mat[,cell_idx]

# Read genomic interval ref with chr and abspos end added
ref_tab<-read.table('/home/sonas/copykat/misc/copyKat_220kb_genomic_fragments_table.txt',sep='\t')

# Add genomic intervals to this subset matrix
rec_mat2<-cbind(ref_tab,rec_mat)

saveRDS(rec_mat2,paste0(dir,'matrix/',sname,'_k10_recurrent_subset_matrix.rds'))

## Step2: Genomic interval subset ##


# Option1: Do this for a range of values and evaluate to choose final cutoff #
# Option2: Do this for specific chosen cutoffs #


##  Option1 ##

thresh_out<-'/home/sonas/copykat/genes/gain_loss/thresholds/'

# Range of cutoff - chosen based on copyKat heatmap colorkey histogram/density plot
cutoff_g<-seq(from=0.015,to=0.055,by=0.005)
cutoff_l<-cutoff_g*(-1)

if(var=='gain')
{#Gain
lapply(1:length(cutoff_g), function(x){subset_by_cutoff(rec_mat2,usr_cutoff=cutoff_g[x], var='gain',out=thresh_out, sname=sname) })
bl}else if(var=='loss')
{
#Loss
lapply(1:length(cutoff_l), function(x){subset_by_cutoff(rec_mat2,usr_cutoff=cutoff_l[x], var='loss', out=thresh_out, sname=sname) })
}else
{cat('Specify correct variation either gain or loss \n')}


## Plot to evaluate thresholds ##

plot_out<-'/home/sonas/copykat/genes/gain_loss/thresh_plots/'

g<-list.files(thresh_out,pattern='gain',full.names = TRUE)
nm<-parse_number(g)
gains<-lapply(g,readRDS)

p75<-lapply(gains,'[[',1)
p70<-lapply(gains,'[[',2)
p60<-lapply(gains,'[[',3)
p50<-lapply(gains,'[[',4)
p25<-lapply(gains,'[[',5)


# Now get genomic coord table for each set of cutoffs
g_p75_coord<-lapply(p75,nrow) %>% unlist()
g_p70_coord<-lapply(p70,nrow) %>% unlist()
g_p60_coord<-lapply(p60,nrow) %>% unlist()
g_p50_coord<-lapply(p50,nrow) %>% unlist()
g_p25_coord<-lapply(p25,nrow) %>% unlist()

df75<-data.frame("genomic_segments"=g_p75_coord,"percentile"=rep('75p',length(g_p75_coord)))
df70<-data.frame("genomic_segments"=g_p70_coord,"percentile"=rep('70p',length(g_p70_coord)))
df60<-data.frame("genomic_segments"=g_p60_coord,"percentile"=rep('60p',length(g_p60_coord)))
df50<-data.frame("genomic_segments"=g_p50_coord,"percentile"=rep('50p',length(g_p50_coord)))
df25<-data.frame("genomic_segments"=g_p25_coord,"percentile"=rep('25p',length(g_p25_coord)))

df_g<-rbind(df75,df70,df60,df50,df25)
label<-seq(0.015,0.055,by=0.005)
df_g$cutoff<-rep(label,5)

col<-c('#7CAE00','dodgerblue4','#F8766D','red','#00BFC4')


g1<-ggplot(df_g,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=percentile)) +
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',subtitle = '[gain]') +
  scale_y_continuous(breaks = round(seq(min(df_g$genomic_segments), max(df_g$genomic_segments), by = 500),1)) +
  scale_x_continuous(breaks = round(seq(min(df_g$cutoff), max(df_g$cutoff), by = 0.005),digits = 3)) +
  scale_color_manual(values=col)+
  theme_bw()

ggsave(paste0(plot_out,'primary_recurrent_genomic_segments_cnt_by_cutoff_gain.png'),width=8,height = 8,units='in',g1)

write.xlsx(x = df_g,file = paste0(plot_out,'primary_recurrent_gains_cell_cnt_by_cutoff.xlsx'))

rem<-grep(25,df_g$percentile)
df_sub_g<-df_g[-rem,]
  
g2<-ggplot(df_sub_g,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=percentile)) +
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',subtitle = '[gain]') +
  scale_y_continuous(breaks = round(seq(min(df_sub_g$genomic_segments), max(df_sub_g$genomic_segments), by = 100),1)) +
  scale_x_continuous(breaks = round(seq(min(df_sub_g$cutoff), max(df_sub_g$cutoff), by = 0.005),digits = 3)) +
   scale_color_manual(values=col[2:5]) +
  theme_bw()

ggsave(paste0(plot_out,'primary_recurrent_gains_genomic_segments_cnt_by_cutoff_25p_removed_gain.png'),width=8,height = 8,units='in',g2)


### Loss ###

l<-list.files(thresh_out,pattern='loss',full.names = TRUE)
nm<-parse_number(l)

loss<-lapply(l,readRDS)

p75<-lapply(loss,'[[',1)
p70<-lapply(loss,'[[',2)
p60<-lapply(loss,'[[',3)
p50<-lapply(loss,'[[',4)
p25<-lapply(loss,'[[',5)


l_p75_coord<-lapply(p75,nrow) %>% unlist()
l_p70_coord<-lapply(p70,nrow) %>% unlist()
l_p60_coord<-lapply(p60,nrow) %>% unlist()
l_p50_coord<-lapply(p50,nrow) %>% unlist()
l_p25_coord<-lapply(p25,nrow) %>% unlist()

df75<-data.frame("genomic_segments"=l_p75_coord,"percentile"=rep("75p",length(l_p75_coord)))
df70<-data.frame("genomic_segments"=l_p70_coord,"percentile"=rep("70p",length(l_p70_coord)))
df60<-data.frame("genomic_segments"=l_p60_coord,"percentile"=rep("60p",length(l_p60_coord)))
df50<-data.frame("genomic_segments"=l_p50_coord,"percentile"=rep("50p",length(l_p50_coord)))
df25<-data.frame("genomic_segments"=l_p25_coord,"percentile"=rep("25p",length(l_p25_coord)))

df_l<-rbind(df75,df70,df60,df50,df25)
label<-seq(0.015,0.055,by=0.005)
label_l<-label*(-1)
df_l$cutoff<-rep(label_l,5)

col<-c('#7CAE00','dodgerblue4','#F8766D','red','#00BFC4')

g3<-ggplot(df_l,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=factor(percentile))) +
  ggtitle(label = 'Genomic_segments count for CNV values <cutoff',subtitle = '[loss]') +
  scale_y_continuous(breaks = round(seq(min(df_l$genomic_segments), max(df_l$genomic_segments), by = 500),1)) +
  scale_x_continuous(breaks = round(seq(min(df_l$cutoff), max(df_l$cutoff), by = 0.005),3)) +
  scale_color_manual(values=col) +
  theme_bw()

ggsave(paste0(plot_out,'primary_recurrent_genomic_segments_cnt_by_cutoff_loss.png'),width=8,height = 8,units='in',g3)

write.xlsx(x = df_l,file = paste0(plot_out,'primary_recurrent_loss_cell_cnt_by_cutoff.xlsx'))

rem<-grep(25,df_l$percentile)
df_sub_l<-df_l[-rem,]
  
g4<-ggplot(df_sub_l,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=factor(percentile))) +
  ggtitle(label = 'Genomic_segments count for CNV values <cutoff',subtitle = '[loss]') +
  scale_y_continuous(breaks = round(seq(min(df_sub_l$genomic_segments), max(df_sub_l$genomic_segments), by = 100),1)) +
  scale_x_continuous(breaks = round(seq(min(df_sub_l$cutoff), max(df_sub_l$cutoff), by = 0.005),3)) +
  scale_color_manual(values=col[2:5]) +
  theme_bw()

ggsave(paste0(plot_out,'primary_recurrent_genomic_segments_cnt_by_cutoff_25p_removed_loss.png'),width=8,height = 8,units='in',g4)


# We choose 0.03 and -0.03 as gain and loss cutoff with cell count cutoff as 75% ##

loss<-readRDS()
cutoff_l<-(-0.03)
df_l<-loss[[1]] #  81 segments

gain<-readRDS()
cutoff_g<-0.03
df_g<-gain[[1]] #  76 segments

# 199 genes for 81 coordinates
rownames(df_l)<-NULL

## Option 2 ##

# Need to add code for this section #


## Step3: Annotate genomic interval subset table ##

anno_out<-'/home/sonas/copykat/genes/gain_loss/annotated_table/'

anno_l<-annotate_genomic_ranges(df_l,var='loss',cutoff=cutoff_l,plot=FALSE,bw=200,out=anno_out,prefix=paste0('primary_recurrent_loss_cutoff',cutoff_l),save=FALSE,chr_size=chr_size)

# Sort
anno_l<-anno_l %>% arrange(.,seqnames,start,end,strand)

write.xlsx(anno_l,paste0(anno_out,'annotated_genomic_ranges_loss_50p_cutoff',cutoff_l,'.xlsx'))

# 604 genes for 76 coordinates

rownames(df_g)<-NULL
anno_g<-annotate_genomic_ranges(df_g,var='gain',cutoff=cutoff_g,plot=FALSE,bw=200,out=anno_out,prefix=paste0('primary_recurrent_gain_cutoff',cutoff_g),save=FALSE,chr_size=chr_size)

anno_g<-anno_g %>% arrange(.,seqnames,start,end,strand)

write.xlsx(anno_g,paste0(anno_out,'annotated_genomic_ranges_gain_75p_cutoff',cutoff_g,'.xlsx'))


