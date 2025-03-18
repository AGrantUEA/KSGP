setwd("E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts")

library(phyloseq)
library(vegan)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stats)

taxonomy=read.table("KSGP_v3.1_hiera_BLAST.txt",sep="\t",header=TRUE)

#setwd("E:/Documents/OneDrive - University of East Anglia/Abdullah/Lotus/AR/final_blasts")

#nearest=read.table("karst.uc",sep="\t")
nearestgtdb=read.table("ArchaeaGTDBplus_v3.1.uc",sep="\t", header = FALSE)
both=read.table("ArchaeaKSGP_v3.1.uc",sep="\t",header=FALSE)
#silva=read.table("silva2.uc",sep="\t",header=FALSE)
#both2=data.frame(both$V9,both$V4,both$V10)

annotated=merge(taxonomy,nearestgtdb,by.x="OTU",by.y="V9",suffixes=c("","_g"))

annotated2=merge(annotated,both,by.x="OTU",by.y="V9",suffixes=c("_g","_k"))

#annotated3=merge(annotated2,silva,by.x="OTU",by.y="V9")
annotated4=annotated2[which(annotated2$Domain=="Archaea"),]
table(annotated4$Phylum)

toplot1=cbind(annotated4$Phylum,as.numeric(annotated4$V4_g),1)
toplot2=cbind(annotated4$Phylum,as.numeric(annotated4$V4_k),2)
toplot=data.frame(rbind(toplot1,toplot2))
boxplot((as.numeric(toplot$X2))~toplot$X3*toplot$X1,ylim=c(75,100),col=(c("gold","darkgreen")),cex.axis=0.5,
        at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,41,43,44,46,47,49,50,52,53),xlab="",ylab="Percentage similarity")
lines(c(0,39),c(91.,91),lty=3,col="blue")
lines(c(0,39),c(88.,88),lty=2,col="red")

exclude=c("B1Sed10-29","EX4484-52","Hydrothermarchaeota","Methanobacteriota","Methanobacteriota_B","Nanohaloarchaeota","SpSt-1190")
toplotx=toplot[which(!toplot$X1 %in% exclude),]
pdf("final/boxplot.pdf",width=12)
boxplot((as.numeric(toplotx$X2))~toplotx$X3*toplotx$X1,ylim=c(75,100),col=(c("gold","darkgreen")),cex.axis=0.5,
        at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32),xlab="",ylab="Percentage similarity")
lines(c(0,39),c(91.,91),lty=3,col="blue")
lines(c(0,39),c(88.,88),lty=2,col="red")
dev.off()

labs=rownames(table(toplotx$X1))
labs[1]="Unknown"
xpos=c(1.5,5.5,9.5,13.5,17.5,21.5,25.5,29.5,33.5,37.5,41.5)
pdf("final/Figure4b.pdf",width=12)
par(mar=c(10.1, 4.1, 4.1, 2.1))
boxplot((as.numeric(toplotx$X2))~toplotx$X3*toplotx$X1,ylim=c(75,100),col=(c("gold","darkgreen")),cex.axis=1.0,
        at=c(1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38,41,42),xlab="",ylab="Percentage similarity", xaxt="n")
text(xpos,70,labels=labs, xpd = NA,srt = 55)
lines(c(0,43),c(91.,91),lty=3,col="blue",lwd=2)
lines(c(0,43),c(88.,88),lty=2,col="red",lwd=2)
dev.off()




sp99=read.table("OBEP01_archaea.uc", sep="\t")
str(sp99)
tfamily=read.table("OBEP01_archaea_fam.uc", sep="\t")
tfamily=tfamily[which(tfamily$V1!="C"),]
str(tfamily)
torder=read.table("OBEP01_archaea_order.uc", sep="\t")
torder=torder[which(torder$V1!="C"),]

str(torder)
tclass=read.table("OBEP01_archaea_class.uc", sep="\t")
tclass=tclass[which(tclass$V1!="C"),]


str(tclass)
phylum=read.table("OBEP01_archaea_phylum2.uc", sep="\t")
phylum=phylum[which(phylum$V1!="C"),]
str(phylum)
sintax=read.table("OBEP01_vs_GTDB_archaea2",sep="\t")
combined=as.data.frame(cbind(phylum$V9,phylum$V2,tclass$V2,torder$V2,tfamily$V2))
annotated=merge(nearest,combined,by.x="V10",by.y="V1")
annotated2=annotated[order(annotated$V2.y, annotated$V3.y, annotated$V4.y, annotated$V5.y, annotated$V9),]
annotated2[c(10,5,11:14)]
annotated3=merge(annotated2,sintax,by.x="V10",by.y="V1" )
annotated3=annotated3[order(annotated3$V2.y, annotated3$V3.y, annotated3$V4.y, annotated3$V5.y, annotated3$V9),]
annotated3[c(10,5,11:21)]
annotated4=merge(annotated3,nearestgtdb,by.x="V9",by.y="V9",suffixes = c("","z"))
table(annotated4$V3)
table(annotated4$V12)
table(annotated4$V4)
table(annotated4$V13)
table(annotated4$V5)
table(annotated4$V14)
annotated5=merge(annotated4,both2,by.x="V9",by.y="both.V9",suffixes = c("","b"))
write.csv(annotated4,"annotated4.csv")
boxplot(as.numeric(annotated4$V4.x))
boxplot(as.numeric(annotated4$V4z))
boxplot((as.numeric(annotated5$both.V4)-as.numeric(annotated5$V4z))~annotated5$V3,ylim=c(-10,25),cex.axis=0.7)
boxplot((as.numeric(annotated5$both.V4)-as.numeric(annotated5$V4z))~annotated5$V12,ylim=c(-10,25),cex.axis=0.5)
boxplot((as.numeric(annotated5$both.V4))~annotated5$V3,cex.axis=0.5)
boxplot((as.numeric(annotated5$V4z))~annotated5$V3,ylim=c(75,100),cex.axis=0.5)
toplot1=cbind(annotated5$V3,as.numeric(annotated5$V4z),1)
toplot2=cbind(annotated5$V3,as.numeric(annotated5$both.V4),2)
toplot=data.frame(rbind(toplot1,toplot2))
boxplot((as.numeric(toplot$X2))~toplot$X3*toplot$X1,ylim=c(75,100),col=(c("gold","darkgreen")),cex.axis=0.5,
         at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38),xlab="",ylab="Percentage similarity")
lines(c(0,39),c(91.,91),lty=3,col="blue")
lines(c(0,39),c(88.,88),lty=2,col="red")
histogram(as.numeric(annotated4$V4.x)-as.numeric(annotated4$V4z))

xtabs(~annotated4$V3+annotated4$V12)
xtabs(~annotated4$V4+annotated4$V13)
write.csv(annotated3[c(10,5,11:21)],"annotated.csv")
xtabs(~annotated3$V5+annotated3$V2.y)

length(table(annotated$V2.y))
sort(table(annotated$V2.y),decreasing=TRUE)
length(table(annotated$V3.y))
sort(table(annotated$V3.y),decreasing=TRUE)
sort(table(annotated$V4.y),decreasing=TRUE)
sort(table(combined$V2),decreasing=TRUE)
sort(table(combined$V3),decreasing=TRUE)

histogram(as.numeric(combined$V7))
sort(table(combined$V6),decreasing=TRUE)





nearest=read.table("silva138_karst_archaea_centroids_vs_pr2_GTDB.uc",sep="\t")
sp99=read.table("silva138_karst_archaea.uc", sep="\t")
str(sp99)
tfamily=read.table("silva138_karst_archaea_fam2.uc", sep="\t")
tfamily=tfamily[which(tfamily$V1!="C"),]
str(tfamily)
torder=read.table("silva138_karst_archaea_order2.uc", sep="\t")
torder=torder[which(torder$V1!="C"),]

str(torder)
tclass=read.table("silva138_karst_archaea_class2.uc", sep="\t")
tclass=tclass[which(tclass$V1!="C"),]


str(tclass)
phylum=read.table("silva138_karst_archaea_phylum2.uc", sep="\t")
phylum=phylum[which(phylum$V1!="C"),]
str(phylum)
combined=as.data.frame(cbind(phylum$V9,phylum$V2,tclass$V2,torder$V2,tfamily$V2,nearest$V10,nearest$V4))
sort(table(combined$V2),decreasing=TRUE)
sort(table(combined$V3),decreasing=TRUE)

histogram(as.numeric(combined$V7))
sort(table(combined$V6),decreasing=TRUE)
