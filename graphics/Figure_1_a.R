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

toplot=read.csv("cumulative_plots_a.csv")
toplot
n=length(toplot$file)
load("AR_keep_v2_phyloseq.Rdata")

counts=rowSums(physeq@otu_table)


pdf("final/figure_1a.pdf",width=8,height=8) 

for (j in 1:n){
 
  if (toplot$type[j]==1){
    data=read.table(toplot$file[j],header = FALSE,sep="\t")
    pct=sort(as.numeric(data$V4),decreasing=TRUE,na.last=TRUE)

    nn=length(pct)
    nnn=100*c(1:nn)/nn
    if (j==1){
        cat("\n",j," ",toplot$file[1],"\n") 
  
        plot (nnn,pct,type="l",lty=toplot$lty[j],col=toplot$colour[j],lwd=3,ylim=c(75.,100.),ylab="Similarity of best match (%)", xlab="Cumulative percentage of OTUs/Sequences")
        }else {
        cat("\n",j," just line ",toplot$file[j],toplot$colour[j])
        lines(nnn,pct,lty=toplot$lty[j],col=toplot$colour[j],lwd=3)

        }
  }
  if (toplot$type[j]==2){
    cat("\n",j," ",toplot$file[1],"\n") 
  
    data=read.table(toplot$file[j],header = FALSE,sep="\t")
  
    data$V7[data$V8<50]=NA
    pct=sort(data$V7,decreasing=TRUE,na.last=TRUE)
    blastlen=length(pct)
    lines(nnn[1:blastlen],pct,lty=toplot$lty[j],col=toplot$colour[j],lwd=3)
  
    }
  if (toplot$type[j]==3){
    cat(j," ",toplot$file[1],"\n") 
  
    data=read.table(toplot$file[j],header = FALSE,sep="\t")
    
    merged=merge(data,counts,by.x="V9",by.y="row.names")
    merged$sim=as.numeric(merged$V4)
    sorted=merged[order(-merged$sim, na.last = TRUE),]
    sorted$cum=cumsum(sorted$y)
    nn4=(100*sorted$cum)/max(sorted$cum)
    nn4[1]=0
    lines(nn4,sorted$V4,lty=toplot$lty[j],col=toplot$colour[j],lwd=3)
    
  }

  }



legend(-3,82,legend=toplot$caption,col=toplot$colour,lwd=3,lty=toplot$lty,cex=1.0,bty="n")
dev.off()


