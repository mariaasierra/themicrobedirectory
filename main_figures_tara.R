library(ggpie)
library(ggplot2)
library(dplyr)
library(ape)
library(ggtree)
library(eulerr)
library(ggpubr)
library(reshape2)
library(tidyr)
library(tidyverse)

#Import external datasets
kraken2_datasets=read.csv("~/External_datasets/data_examples_tmd.tsv", header = F,sep='\t') #Taxonomy file from external datasets
colnames(kraken2_datasets)=c("Sample_ID","taxid","taxa", "species","relative_abundance","clade_reads", "taxon_reads", "rank_code") #Change column names
metadata_datasets=read.csv("~/External_datasets/metadata_all_datasets.csv", header = T)

#Import TMD annotations

bac_tmd=read.csv("~/Data/Bacteria_data.csv", header = T)
fungi_tmd=read.csv("~/Data/Fungi_data.csv", header = T)
virus_tmd=read.csv("~/Data/Virus_data.csv", header = T)
algae_tmd=read.csv("~/Data/Algae_data.csv", header = T)


### TARA Oceans ###

#Extremophiles
tara_kraken_sp=tara_kraken_sp[-grep(("c__Mammalia"), tara_kraken$taxa),] #Remove any mammalian assigned reads
tara_kraken_metadata=merge(tara_kraken_sp,tara_metadata) #Merge with samples metadata 
bac_tmd_in_tara=bac_tmd[bac_tmd$species %in%  unique(tara_kraken$species),] #Get bacteria present TARA species list annotated in TMD

xm_tara=select(bac_tmd_in_tara,species,extremophile_type)#Select columns of dataframe
xm_tara=na.omit(xm_tara)
xm_tara <- xm_tara[-which(xm_tara$extremophile_type == ""), ] #Get list of extremophile types present

lev=levels(factor(xm_tara$extremophile_type))#Split extremophile type column
lev=unique(unlist(strsplit(lev, " "))) #Unique types levels of extremophiles
result_extreme <- matrix(data = "", nrow = length(xm_tara$extremophile_type), ncol = length(lev)) #create empty matrix with df dimentions
char.var = xm_tara$extremophile_type

for (i in 1:length(lev)) {
  result_extreme[grep(lev[i], char.var, fixed = F), i] <- "Present"
} #Presence matrix for each extremophile type
result_extreme <- data.frame(result_extreme, stringsAsFactors = TRUE)
colnames(result_extreme) <- lev
xm_tara$extremophile_type=NULL
xm_tara_total = cbind(xm_tara,result_extreme) #Bind new table to principal TARA tmd-xmp table

xm_type=dplyr::select(xm_tara_total,species, lev) #type of extremophile
xm_type_sum=xm_tara_total %>% group_by(species) %>% summarise_all(funs(trimws(paste(unique(.), collapse = " ")))) 
xm_type_sum=as.data.frame(xm_type_sum)
xm_type_sum=na.omit(xm_type_sum)
rownames(xm_type_sum)=xm_type_sum$species
xm_type_sum=xm_type_sum[,-1]

w <- which(xm_type_sum=="Present",arr.ind=TRUE) #Convert "present" to the type of extremophile per column for easier plotting
xm_type_sum[w] <- names(xm_type_sum)[w[,"col"]]

#write.table(row.names(xm_type_sum), "xmp_tara_taxa.txt")
#cladogram_metasub.py #do not run
xm_tree=read.tree("tree_extremophiles.txt") #Import tree generated in phython
xm_tree$tip.label<-gsub("_"," ",xm_tree$tip.label) #convert underscore to space in tree names

xp_tree_plot=ggtree(xm_tree,layout = "fan",size=0.2, branch.length = "none", open.angle =19) + 
  geom_tiplab(size=1.8,  offset = 7)


type_extremophile_heatmap=gheatmap(xp_tree_plot, xm_type_sum, width=0.65, font.size=3,hjust = 1,colnames_angle=90)+
  scale_fill_brewer(palette = "Paired", na.value="gray95") + theme(legend.position = "bottom")


type_extremophile_heatmap







