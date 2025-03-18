setwd("E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts")



data=read.csv("Archaea_classifications.csv")

pdf("final/Figure4ax.pdf")
par(mar=c(5.1, 6.1, 4.1, 2.1))
plot(data$GTDB,type="l",col="orange",lwd=5,lty=3,ylab="Percentage Classified at Taxonomic Level",cex.lab=2.0, xaxt = "n",xlab="Taxonomic Level",ylim=c(0,100),cex.axis=1.5)
axis(1, at = seq(1,6, by = 1),
     labels = data$classified,cex.axis=1.3)
lines(data$GTDBplus,col="orange",lwd=5)
lines(data$Silva,col="red",lwd=5)
lines(data$Sintax,col="blue",lwd=5)
lines(data$LCA,col="blue",lwd=5,lty=3)
lines(data$KSGPplus,col="blue",lwd=5,lty=4)
lines(data$GG2,col="green",lwd=5)


legend(4.3,97,c("raw GTDB","GTDB+  ","KSPG Sintax","KSGP LCA","KSGP+","SILVA","GreenGenes2"),lwd=5,col=c("orange","orange","blue","blue","blue","red","green"),lty=c(3,1,1,3,4,1,1),bty="n")
dev.off()


