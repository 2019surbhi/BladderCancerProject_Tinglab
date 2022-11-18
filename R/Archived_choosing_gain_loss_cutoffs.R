
# Author: Surbhi Sona #


library(dplyr)
library(readr)
library(ggplot2)
library(openxlsx)
library(readxl)


#path<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17/thresholds1/'

path<-'/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/genes/gain_loss_subset/samples17/thresholds/'

out<-'/home/sonas/copykat/genes/gain_loss/samples17/threshold_plots/'

run_prefix<-'samples17'
var<-'gain'
var<-'loss'

f<-list.files(path,pattern=var,full.names = TRUE)
nm<-parse_number(gsub('samples17','',basename(f)))
thresh<-lapply(f,readRDS)

 p75<-lapply(thresh,'[[',1)
 p70<-lapply(thresh,'[[',2)
 p60<-lapply(thresh,'[[',3)
 p50<-lapply(thresh,'[[',4)
 p25<-lapply(thresh,'[[',5)


# Now get genomic coord table for each set of cutoffs
 p75_coord<-lapply(p75,nrow) %>% unlist()
 p70_coord<-lapply(p70,nrow) %>% unlist()
 p60_coord<-lapply(p60,nrow) %>% unlist()
 p50_coord<-lapply(p50,nrow) %>% unlist()
 p25_coord<-lapply(p25,nrow) %>% unlist()


 df75<-data.frame("genomic_segments"=p75_coord,"percentile"=rep('75p',length(p75_coord)))
 df70<-data.frame("genomic_segments"=p70_coord,"percentile"=rep('70p',length(p70_coord)))
 df60<-data.frame("genomic_segments"=p60_coord,"percentile"=rep('60p',length(p60_coord)))
 df50<-data.frame("genomic_segments"=p50_coord,"percentile"=rep('50p',length(p50_coord)))
 df25<-data.frame("genomic_segments"=p25_coord,"percentile"=rep('25p',length(p25_coord)))
 
 df<-rbind(df75,df70,df60,df50,df25)


basename(f) %>% strsplit(.,split='_')

#label<-seq(0.015,0.055,by=0.005)
#label<-c(0.02,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.070)
#label<-label*(-1)

label<-nm
df$cutoff<-rep(label,5)

col<-c('#7CAE00','dodgerblue4','#F8766D','red','#00BFC4')


g1<-ggplot(df,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=factor(percentile))) + 
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',subtitle = paste0('[',var,']'))+
  scale_y_continuous(breaks = round(seq(min(df$genomic_segments), max(df$genomic_segments), by = 500),1)) +
  #scale_x_continuous(breaks = round(seq(min(df$cutoff), max(df$cutoff), by = 0.005),digits = 3)) +
  scale_color_manual(values=col)+
  theme_bw() 

ggsave(paste0(out,run_prefix,'_genomic_segments_cnt_by_cutoff_',var,'.png'),width=8,height = 8,units='in',g1)

write.xlsx(x = df,file = paste0(out,run_prefix,'_',var,'_cell_cnt_by_cutoff.xlsx'))

rem<-grep(25,df$percentile)
df_sub<-df[-rem,]

g2<-ggplot(df_sub,aes(x=cutoff,y=genomic_segments,group=percentile)) +
  geom_line(aes(color=percentile)) + 
  ggtitle(label = 'Genomic segments count for CNV values >cutoff',
          subtitle = paste0('[',var,']')) +
  scale_y_continuous(breaks = round(seq(min(df_sub$genomic_segments), max(df_sub$genomic_segments), by = 100),1)) +
  scale_x_continuous(breaks = round(seq(min(df_sub$cutoff), max(df_sub$cutoff), by = 0.005),digits = 3)) +
  scale_color_manual(values=col[2:5]) +
  theme_bw()

ggsave(paste0(out,run_prefix,var,'_cell_cnt_by_cutoff_25p_removed.png'),width=8,height = 8,units='in',g2)


