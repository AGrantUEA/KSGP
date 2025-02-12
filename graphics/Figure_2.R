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

ARGPM=read.table('E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts/ArchaeaGTDBplus_v3.uc',header = FALSE,sep="\t")
str(ARGPM)
xmerged=merge(ARGPM,taxa,by.x="V9",by.y="row.names")
tmerged=xmerged[which(xmerged$Domain!="?"),]
GPMmerged=tmerged[which(tmerged$Domain=="Archaea"),]


ARGPMSK=read.table('E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts/ArchaeaKSGP_v2.uc',header = FALSE,sep="\t")
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



pdf("figure_2b.pdf",width=8,height=8) 

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
h1=ggplot(combined2,aes(V4_g))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+xlab("Similarity to best match in GTDB")+ylab("") +theme_bw()+scale_y_continuous(limits=c(0.,1250),expand = c(0,0))+ theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))+theme(axis.text.y =element_blank())   

h2=ggplot(combined2,aes(V4))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+ylab("")+xlab("")+coord_flip()+theme_bw()+scale_y_continuous(limits=c(0.,1500),expand = c(0,0))+ theme(axis.title = element_text(size = 20))  +theme(axis.text = element_text(size = 15)) +theme(axis.text.x =element_blank())
grid.arrange(p1,h2,h1,ncol=2,widths=c(4,1),heights=c(4,1))

dev.off()

pdf("figure_2a.pdf",width=8,height=8) 


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
grid.arrange(p1,h2,h1,ncol=2,widths=c(4,1),heights=c(4,1))

dev.off()






BacteriaGPMSKl=read.table('BacteriaGPMSKl.uc',header = FALSE,sep="\t")
Bacteriasilva=read.table('Bacteriasilva.uc',header = FALSE,sep="\t")
Bacteriakarst=read.table('Bacteriakarst.uc',header = FALSE,sep="\t")



BacteriaGPMSKlpct=sort(as.numeric(BacteriaGPMSKl$V4),decreasing=TRUE,na.last=TRUE)



Bacteriasilvapct=sort(as.numeric(Bacteriasilva$V4),decreasing=TRUE,na.last=TRUE)


Bacteriakarstpct=sort(as.numeric(Bacteriakarst$V4),decreasing=TRUE,na.last=TRUE)

BacteriaGPM=read.table('BacteriaGPM2.uc',header = FALSE,sep="\t")

BacteriaGPMpct=sort(as.numeric(BacteriaGPM$V4),decreasing=TRUE,na.last=TRUE)

nn2=length(BacteriaGPMSKlpct)
nnn2=100*c(1:nn2)/nn2

plot (nnn2,BacteriaGPMpct,type='l',col='orange',lwd=3,ylim=c(75.,100.),ylab="Similarity of best match (%)", xlab="Cumulative percentage of OTUs/Sequences")


legend(-3,85,legend=c("KPSG","Karst","NCBI nt","GTDB+","Silva","Greengenes2","KPSG by abundance"),col=c("blue","black","purple","orange","brown","green","blue"),lwd=4,lty=c(1,1,1,1,1,1,3),cex=1.0,bty="n")


lines(nnn2,BacteriaGPMSKlpct,col="blue",lty=1,lwd=3)
lines(nnn2,Bacteriasilvapct,col="brown",lty=1,lwd=3)
lines(nnn2,Bacteriakarstpct,col="black",lty=1,lwd=3)



BacteriaGG2_usl2=read.table('BacteriaGG2_usl2_2.uc',header = FALSE,sep="\t")
str(BacteriaGG2_usl2)
BacteriaGG2_usl2pct=sort(as.numeric(BacteriaGG2_usl2$V4),decreasing=TRUE,na.last=TRUE)
lines(nnn2,BacteriaGG2_usl2pct,col="green",lty=1,lwd=3)


Bacteriablast=read.table("515_otus_blast_c",header = FALSE,sep="\t")
str(Bacteriablast)


Bacteriablast$V7[Bacteriablast$V8<50]=NA
Bacteriablastpct=sort(Bacteriablast$V7,decreasing=TRUE,na.last=TRUE)

nnn2x=length(Bacteriablastpct)
lines(nnn2[1:nnn2x],Bacteriablastpct,col="purple",lwd=3)

load("GPMSK_loot2_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])
table(physeq@tax_table@.Data[,5])
table(physeq@tax_table@.Data[,6])



load("515_phyloseq.Rdata")
counts=rowSums(physeq@otu_table)



merged=merge(BacteriaGPMSKl,counts,by.x="V9",by.y="row.names")
merged$sim=as.numeric(merged$V4)
sorted=merged[order(-merged$sim, na.last = TRUE),]
sorted$cum=cumsum(sorted$y)
nn4=(100*sorted$cum)/max(sorted$cum)
nn4[1]=0
lines(nn4,sorted$V4,col="blue",lwd=3,lty=3)




table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])
table(physeq@tax_table@.Data[,5])
table(physeq@tax_table@.Data[,6])



load("SLV_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

load("AR_GG2_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])


load("AR_GPM_loot2_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])


load("AR_GPMSK_loot2_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])

load("GPM_phyloseq.Rdata")


table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

load("515_GPM_phyloseq.Rdata")
load ("GPM_loot2_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])[1]
table(physeq@tax_table@.Data[,3])[1]
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

load("515_GPMSK_final_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])[1]
table(physeq@tax_table@.Data[,3])[1]
table(physeq@tax_table@.Data[,4])[1]
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

load("515_SLV_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])[1]
table(physeq@tax_table@.Data[,3])[1]
table(physeq@tax_table@.Data[,4])[1]
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

load("515_GPMSK_old_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])[1]
table(physeq@tax_table@.Data[,3])[1]
table(physeq@tax_table@.Data[,4])[1]


load("515_GG2_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])[1]
table(physeq@tax_table@.Data[,3])[1]
table(physeq@tax_table@.Data[,1])
table(physeq@tax_table@.Data[,2])
table(physeq@tax_table@.Data[,3])
table(physeq@tax_table@.Data[,4])

gg2tab=read.table("515_GG2_hiera_BLAST.txt",header=TRUE,sep="\t")

gpmsk2tab=read.table("GPMSK_loot2_hiera_BLAST.txt",header=TRUE,sep="\t")

slvtab=read.table("515_SLV_hiera_BLAST.txt",header=TRUE,sep="\t")
merged=merge(gpmsk2tab,gg2tab,by="OTU",suffixes=c(".gpmsk",".gg"))
merged2=merge(merged,slvtab,by="OTU",suffixes=c("",".slv"))
GPMSKtop=read.table("515_GPMSK_tax.0.first",sep="\t")
merged3=merge(merged2,GPMSKtop,by.x="OTU",by.y="V1")
GPMtop=read.table("515_GPM_tax.0.first",sep="\t")
merged4=merge(merged3,GPMtop,by.x="OTU",by.y="V1")
nophylum=merged4[which(merged2$Phylum.gpmsk=="?"),]
table(nophylum$Domain.gpmsk)
ggphyfreq=sort(table(nophylum$Phylum.gg),decreasing = TRUE)
slvphyfreq=sort(table(nophylum$Phylum),decreasing = TRUE)
xtabs(~Phylum.gg+Phylum,data=nophylum)
slvphyfreq[1:10]
ggphyfreq[1:10]
commonslv=row.names(slvphyfreq[1:10])
commonnophylum=nophylum[which(nophylum$Phylum %in% commonslv),]
xtabs(~Phylum.gg+Phylum,data=commonnophylum)

ggphyfreq[1:10]
commongg=row.names(ggphyfreq[1:10])
commonnophylum2=nophylum[which(nophylum$Phylum.gg %in% commongg),]
xtabs(~Phylum+Phylum.gg,data=commonnophylum2)
hist(nophylum$V3.x)
hist(nophylum$V3.y)
plot(nophylum$V3.y,nophylum$V3.x)
plot(merged4$V3.y,merged4$V3.x)


xdata=merged4
xdata=nophylum
p = ggplot(xdata,aes(x=V3.y,y=V3.x)) +
  theme(axis.title = element_text(size = 20))   +
  xlab("") +
  
  ylab("Similarity to best match in GPMSK")+ xlim(73,100)+ylim(73,100) +geom_abline(intercept = 0,slope=1,linewidth=2,color="black")

p1 = p + 
  geom_point(alpha = 0.21, colour="red") + 
  geom_density2d(linewidth=2) + 
  theme_bw()+ theme(legend.position = "none")+
  theme(axis.title = element_text(size = 20)) + 
  theme(axis.text = element_text(size = 15)) +
  geom_smooth(method = "loess", se = FALSE,color="green",linewidth=2) 
#plot(p1)
h1=ggplot(xdata,aes(V3.y))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+xlab("Similarity to best match in GPM")+ylab("") +theme_bw()+scale_y_continuous(limits=c(0.,1250),expand = c(0,0))+ theme(axis.title = element_text(size = 20))+theme(axis.text = element_text(size = 15))+theme(axis.text.y =element_blank())   

h2=ggplot(xdata,aes(V3.x))+geom_histogram(binwidth=1,center=0.5)+xlim(73,100)+ylab("")+xlab("")+coord_flip()+theme_bw()+scale_y_continuous(limits=c(0.,1500),expand = c(0,0))+ theme(axis.title = element_text(size = 20))  +theme(axis.text = element_text(size = 15)) +theme(axis.text.x =element_blank())
grid.arrange(p1,h2,h1,ncol=2,widths=c(4,1),heights=c(4,1))
