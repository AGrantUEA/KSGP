#library(devtools)
#install_github("Russel88/MicEco")

setwd("E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts")

library(phyloseq)
library(vegan)
library(ggtree)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
load("KSGP_v3_keep_phyloseq.Rdata")
str(physeq)

table(physeq@tax_table@.Data[,1])


tipcol=ifelse(physeq@tax_table@.Data[,1]=="Eukaryota","blue",NA)
tipcol=ifelse(physeq@tax_table@.Data[,1]=="Unknown","yellow",tipcol)
tipcol=ifelse(physeq@tax_table@.Data[,1]=="Archaea","red",tipcol)
tipcol=ifelse(physeq@tax_table@.Data[,1]=="?","purple",tipcol)
tipcol=ifelse(physeq@tax_table@.Data[,1]=="","cadetblue1",tipcol)
treex=physeq@phy_tree

pp=ggtree(treex)+geom_tippoint(color=tipcol,size = 1)


ncols=c("purple","cadetblue1","red","blue")
labtext=c("No hits","Unresolved by LCA","Archaea","Eukaryota")
yymin=13000
yinc=500

for (ii in 1:4) {
  pp=pp+annotation_custom(textGrob(labtext[ii],gp=gpar(col=ncols[ii]),hjust=0),xmin=4.0,xmax=6.,  ymin=yymin-yinc,ymax=yymin)
  yymin=yymin-yinc
}
pp=pp+annotation_custom(textGrob("a",gp=gpar(col="black"),hjust=0),xmin=0.01,xmax=0.1,ymin=14000,ymax=14500)
pp
pdf("final/figure_3a.pdf",width=8,height=8) 

pp

dev.off()


load("KSGP_v3_phyloseq.Rdata")

table(physeq@tax_table@.Data[,1])
physeq=subset_taxa(physeq, Domain == "Archaea")

treex=physeq@phy_tree



tipcol4=ifelse(physeq@tax_table@.Data[,2]=="Nanoarchaeota","green","pink")
tipcol4=ifelse(physeq@tax_table@.Data[,2]=="Thermoplasmatota","brown",tipcol4)
tipcol4=ifelse(physeq@tax_table@.Data[,2]=="Thermoproteota","red",tipcol4)
tipcol4=ifelse(physeq@tax_table@.Data[,2]=="Aenigmatarchaeota","yellow",tipcol4)
tipcol4=ifelse(physeq@tax_table@.Data[,2]=="Asgardarchaeota","purple",tipcol4)
tipcol4=ifelse(physeq@tax_table@.Data[,2]=="?","cadetblue1",tipcol4)



library(grid)

ncols=c("red","purple","brown","yellow","green","pink","cadetblue1")
labtext=c("Thermoproteota","Asgardarchaeota","Thermoplasmatota","Aenigmatarchaeota","Nanoarchaeota","Other","Unresolved by LCA")
yymin=5000
yinc=500
pp2=ggtree(treex)+geom_tippoint(color=tipcol4,size = 1)

for (ii in 1:7) {
  pp2=pp2+annotation_custom(textGrob(labtext[ii],gp=gpar(col=ncols[ii]),hjust=0),xmin=0.5,xmax=1.5,  ymin=yymin-yinc,ymax=yymin)
yymin=yymin-yinc
}
pp2=pp2+annotation_custom(textGrob("b",gp=gpar(col="black"),hjust=0),xmin=0.0,xmax=0.05,ymin=11000,ymax=11500)
pp2

pdf("final/figure_3b.pdf",width=8,height=8) 

pp2

dev.off()


pdf("final/figure_3_all.pdf",width=16,height=8) 

grid.arrange(pp,pp2,nrow=1)

dev.off()

