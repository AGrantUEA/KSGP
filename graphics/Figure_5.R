setwd("E:/Documents/OneDrive - University of East Anglia/Abdullah/KSGP scripts")
#setwd("E:/OneDrive - University of East Anglia/Abdullah/Lotus/AR")



library(phyloseq)
library(vegan)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stats)
library(reshape)
library(scales)
library(ggsci)
library(viridis)
#install.packages("dendextend")

#load("phyloseq_KSG_trimmed2.Rdata")
load("KSGP_v3_keep_phyloseq.Rdata")
table(physeq@tax_table@.Data[,1])
physeq2=prune_taxa(physeq@tax_table@.Data[,1]=="Archaea",physeq)
table(physeq2@tax_table@.Data[,2])
phylumK=physeq2@tax_table@.Data[,2]
classK=physeq2@tax_table@.Data[,3]
taxonomyK=data.frame(physeq2@tax_table@.Data)
taxonomyK$OTU=rownames(taxonomyK)
taxonomyK$counts=rowSums(physeq2@otu_table)

#load("GTDB_trimmed_phyloseq.Rdata")
#table(physeq@tax_table@.Data[,2])
#phylumG=physeq@tax_table@.Data[,2]
#classG=physeq@tax_table@.Data[,3]
#taxonomyG=data.frame(physeq@tax_table@.Data)
#taxonomyG$OTU=rownames(taxonomyG)
#taxonomyG$counts=rowSums(physeq@otu_table)

#taxonomyG=taxonomyG %>%
#  mutate(Phylum=case_when(
#    Phylum=="?" ~"Unknown",
#    Phylum=="Aenigmatarchaeota" ~ "Aenigmatarchaeota",
#    Phylum=="Asgardarchaeota" ~ "Asgardarchaeota",
#    Phylum=="Nanoarchaeota" ~ "Nanoarchaeota",
#    Phylum=="Thermoplasmatota" ~ "Thermoplasmatota",
#    Phylum=="Thermoproteota" ~ "Thermoproteota",
#    Phylum=="Halobacteriota" ~ "Halobacteriota",
#    Phylum=="Methanobacteriota" ~ "Methanobacteriota",
#    TRUE~" other"
#  ))

taxonomyK=taxonomyK %>%
  mutate(Phylum=case_when(
    Phylum=="?" ~"Unknown",
    Phylum=="Aenigmatarchaeota" ~ "Aenigmatarchaeota",
    Phylum=="Asgardarchaeota" ~ "Asgardarchaeota",
    Phylum=="Nanoarchaeota" ~ "Nanoarchaeota",
    Phylum=="Thermoplasmatota" ~ "Thermoplasmatota",
    Phylum=="Thermoproteota" ~ "Thermoproteota",
        TRUE~" other"
  ))

#nseqsG=taxonomyG %>% 
#  group_by(Phylum) %>%
#  summarise(value=n())
#nseqsG$variable="GTDB"
#nseqsG2=nseqsG %>% relocate(value,.after=variable)

nseqs=taxonomyK %>% 
  group_by(Phylum) %>%
  summarise(KARST=n(),KARSTN=sum(counts))

stacked=reshape2::melt(nseqs,measure.vars = c("KARST","KARSTN"),id.vars="Phylum")


reads=read.table("phylumlist_v2",header=FALSE,sep="\t")

#reads=read.csv("total_matches_r207_good_archaea.uc.csv",header=FALSE)
str(reads)
table(reads$V1)
reads=reads %>%
  mutate(V1=case_when(
    V1=="?" ~"Unknown",
    V1=="Aenigmatarchaeota" ~ "Aenigmatarchaeota",
    V1=="Asgardarchaeota" ~ "Asgardarchaeota",
    V1=="Nanoarchaeota" ~ "Nanoarchaeota",
    V1=="Thermoplasmatota" ~ "Thermoplasmatota",
    V1=="Thermoproteota" ~ "Thermoproteota",
        TRUE~" other"
  ))

nseqsR=reads %>% 
  group_by(V1) %>%
  summarise(value=n())
nseqsR$variable="RNASeq"
nseqsR2=nseqsR %>% relocate(Phylum=V1) %>% 
   relocate(value,.after=variable)

stacked=rbind(stacked,nseqsR2)

p <- ggplot(stacked, aes(variable, value, fill = Phylum)) +
  geom_bar(position = "fill", stat = "identity") + ylab("Percentage abundance")+xlab("")+
  scale_y_continuous(labels = percent,expand=c(0,0))+theme_classic()+ theme(panel.grid.major = element_blank())+
  scale_fill_brewer(palette="Spectral")+theme(axis.title=element_text(size=20))+
  scale_x_discrete(labels=c("Metabarcoding OTUs","Metabarcoding reads","RNASeq reads"))+
  aes(ymin=0,ymax=100)+theme(axis.text=element_text(size=15))+theme(legend.text=element_text(size=13))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))
   
p


pdf("figure_5.pdf",width=16,height=6.2) 
p
dev.off()
