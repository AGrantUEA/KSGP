setwd('E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts')
library(stringr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

library(lattice)
#install.packages("installr")
#library(installr)
#updateR()

load("AR_keep_v2_phyloseq.Rdata")
taxa=as.data.frame(physeq@tax_table[,1])

ARGPM=read.table('E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts/ArchaeaGTDB_mar25.uc',header = FALSE,sep="\t")
str(ARGPM)
xmerged=merge(ARGPM,taxa,by.x="V9",by.y="row.names")
tmerged=xmerged[which(xmerged$Domain!="?"),]
GPMmerged=tmerged[which(tmerged$Domain=="Archaea"),]


ARGPMSK=read.table('E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts/ArchaeaKSGP_mar25.uc',header = FALSE,sep="\t")
str(ARGPMSK)

xmerged=merge(ARGPMSK,taxa,by.x="V9",by.y="row.names")
tmerged=xmerged[which(xmerged$Domain!="?"),]
GPMSKmerged=tmerged[which(tmerged$Domain=="Archaea"),]


ARsilva=read.table('E:/Documents/OneDrive - University of East Anglia/Abdullah/Lotus/final_database/final_runs/Archaea_silva.uc',header = FALSE,sep="\t")


xmerged=merge(ARsilva,taxa,by.x="V9",by.y="row.names")
tmerged=xmerged[which(xmerged$Domain!="?"),]


combined1=merge(GPMSKmerged,ARGPM,by="V9",suffixes=c("_k","_g"))
combined2=merge(combined1,ARsilva,by="V9",suffixes=c("","_s"))
str(combined2)

combined2$V4_g=as.numeric(combined2$V4_g)
combined2$V4_k=as.numeric(combined2$V4_k)
combined2$V4=as.numeric(combined2$V4)




p = ggplot(combined2,aes(x=V4_g,y=V4)) +
  theme(axis.title = element_text(size = 20))   +
  xlab("") +
  
  ylab("Similarity to best match in SILVA")+ xlim(73,100)+ylim(73,100) +geom_abline(intercept = 0,slope=1,linewidth=2,color="black")

p1 = p + 
  geom_point(alpha = 0.21, colour="red") + 
  geom_density2d(linewidth=2) + 
  theme_bw()+ theme(legend.position = "none")+
  theme(axis.title = element_text(size = 20)) + 
  theme(axis.text = element_text(size = 15)) +
  geom_smooth(method = "loess", se = FALSE,color="green",linewidth=2) 
#plot(p1)
h1=ggplot(combined2,aes(V4_g))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+xlab("Similarity to best match in GTDB")+ylab(".") +theme_bw()+scale_y_continuous(limits=c(0.,1250),expand = c(0,0))+ theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))+theme(axis.text.y =element_blank())   

h2=ggplot(combined2,aes(V4))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+ylab(".")+xlab(".")+coord_flip()+theme_bw()+scale_y_continuous(limits=c(0.,1500),expand = c(0,0))+ theme(axis.title = element_text(size = 20))  +theme(axis.text = element_text(size = 15)) +theme(axis.text.x =element_blank())
xx1<-grid.arrange(p1,h2,h1,ncol=2,widths=c(4,1),heights=c(4,1))

#dev.off()

#pdf("final/figure_2a.pdf",width=8,height=8) 


p = ggplot(combined2,aes(x=V4_g,y=V4_k)) +
  theme(axis.title = element_text(size = 20))   +
  xlab("") +
  
  ylab("Similarity to best match in KSGP")+ xlim(73,100)+ylim(73,100) +geom_abline(intercept = 0,slope=1,linewidth=2,color="black")

p1 = p + 
  geom_point(alpha = 0.21, colour="red") + 
  geom_density2d(linewidth=2) + 
  theme_bw()+ theme(legend.position = "none")+
  theme(axis.title = element_text(size = 20)) + 
  theme(axis.text = element_text(size = 15)) +
  geom_smooth(method = "loess", se = FALSE,color="green",linewidth=2) 
#plot(p1)
h1=ggplot(combined2,aes(V4_g))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+xlab("Similarity to best match in GTDB")+ylab("") +theme_bw()+scale_y_continuous(limits=c(0.,1250),expand = c(0,0))+ theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))+theme(axis.text.y =element_blank())   

h2=ggplot(combined2,aes(V4_k))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+ylab("")+xlab("")+coord_flip()+theme_bw()+scale_y_continuous(limits=c(0.,1500),expand = c(0,0))+ theme(axis.title = element_text(size = 20))  +theme(axis.text = element_text(size = 15)) +theme(axis.text.x =element_blank())
xx2<-grid.arrange(p1,h2,h1,ncol=2,widths=c(4,1),heights=c(4,1))


pdf("final/figure_2all.pdf",width=16,height=8) 

grid.arrange(xx1,xx2, ncol=2, nrow =1)

dev.off()



