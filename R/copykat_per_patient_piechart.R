library(ggplot2)
library(gridExtra)

# Read in per patient proportion

tab<-read_csv('/home/sonas/tingalab/Surbhi/PROJECTS_tinglab_drive/BLADDER/copyKat/prop_tables/BL_191065T_trio_cut4_sample_prop_by_hclust_group.csv')



### Get Pie chart ###
slist<-list()
prlist<-list()
tab2<-tab[,-1]
clades<-paste0('clade',1:4)

for(i in 1:nrow(tab))
{
df<-tab2[i,] %>% t() %>% as.data.frame()
df<-cbind(rownames(df),df)
rownames(df)<-NULL
colnames(df)<-c('sample','cells')
df$p_r<-c('R','P','R')
bp1<- ggplot(df, aes(x="", y=cells, fill=sample))+
  geom_bar(width = 1, stat = "identity")
slist[[i]]<- bp1 + coord_polar("y", start=0) + ggtitle(label = clades[i])
bp2<- ggplot(df, aes(x="", y=cells, fill=p_r))+
  geom_bar(width = 1, stat = "identity")
prlist[[i]]<-bp2 + coord_polar("y", start=0) +ggtitle(label = clades[i])

}

ms<-marrangeGrob(slist,nrow=2,ncol=2)
mpr<-marrangeGrob(prlist,nrow=2,ncol=2)

out<-'/home/sonas/copykat/plots/'

ggsave(filename = paste0(out,'BL_191065trio_k4_pie_by_sample.png'),width = 10,height = 10,units='in',ms)
ggsave(filename = paste0(out,'BL_191065trio_k4_pie_by_PR.png'),width = 10,height = 10,units='in',mpr)

