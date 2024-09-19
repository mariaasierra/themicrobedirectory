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


#### MetaSUB ####

#Spore forming base on coastal city

metasub_metadata=metadata_datasets[metadata_datasets$Project=="Metasub",] #Get only Metasub data
metasub_kraken=kraken2_datasets[kraken2_datasets$Sample_ID %in% metasub_metadata$Sample_ID ,]

metasub_kraken_bac=metasub_kraken[grep(c("Bacteria"), metasub_kraken$taxa,ignore.case=T),]#Get only Bacteria
metasub_kraken_bac_sp=metasub_kraken_bac[metasub_kraken_bac$rank_code=="S",] #Get only species of Bacteria
metasub_kraken_bac_sp=merge(metasub_kraken_bac_sp,metasub_metadata)#Merge with sample metadata


bac_tmd_in_metasub=bac_tmd[bac_tmd$species %in% unique(metasub_kraken_bac_sp$species),] #Get bacterial species present in metasub that are also annotated in TMD
bac_tmd_trim=select(bac_tmd_in_metasub,species,extremophile,extremophile_type, spore_forming) #Trim dataframe
bac_metasub_tmd_metadata=merge(metasub_kraken_bac_sp,bac_tmd_trim) #Merge with sample metadata

bac_metasub_tmd_metadata_spores=as.data.frame(bac_metasub_tmd_metadata[bac_metasub_tmd_metadata$spore_forming=="Yes",]) #Get spore forming bacteria based on TMD
bac_metasub_tmd_metadata_spores=bac_metasub_tmd_metadata_spores[!is.na(bac_metasub_tmd_metadata_spores$species),] #Filter any missing data

bac_metasub_tmd_metadata_spores$clade_reads_relative=100*(bac_metasub_tmd_metadata_spores$clade_reads/
                                                            sum(bac_metasub_tmd_metadata_spores$clade_reads)) #Get relative abundances

p <- ggplot(bac_metasub_tmd_metadata_spores, aes(Beach, clade_reads_relative, fill=Beach))
spore_forming_city=p + geom_boxplot(outlier.shape = NA)+
  geom_point(alpha = 0.4, position=position_jitter( width=.2))+ 
  ylab("Relative abundance of spore forming bacteria")+ xlab("Beach city")+
  stat_compare_means(method = "t.test", size=8) + ylim(0,0.5)+ 
  theme_pubr()+ scale_fill_manual(values = c("tan1", "steelblue2"))

spore_forming_city


#Halophiles in coastal cities

unique(bac_metasub_tmd_metadata$extremophile_type) #Verify extremophile types in TMD
bac_metasub_tmd_metadata_halo=bac_metasub_tmd_metadata[grep("Halophile",bac_metasub_tmd_metadata$extremophile_type),] #Extract halophiles only
bac_metasub_tmd_metadata_halo$clade_reads_relative=100*(bac_metasub_tmd_metadata_halo$clade_reads/sum(bac_metasub_tmd_metadata_halo$clade_reads))#Get relative abundances

g <- ggplot(bac_metasub_tmd_metadata_halo, aes(Beach, clade_reads_relative, fill=Beach))
halophiles_city=g + geom_boxplot(outlier.shape = NA)+
  geom_point(alpha = 0.4, position=position_jitter( width=.2))+ 
  ylab("Relative abundance of Halophiles")+ xlab("Beach city")+
  stat_compare_means(method = "t.test", size=8) + ylim(0,0.5)+ 
  theme_pubr()+ scale_fill_manual(values = c("tan1", "steelblue2"))

halophiles_city


#Heatmap of bacterial host according to TMD

lev=levels(factor(bac_tmd_in_metasub$microbiome_host))#Split microbiome host column
lev=unique(unlist(strsplit(lev, " "))) #unique types levels of microbiome host
result_bac_metasub<- matrix(data = "", nrow = length(bac_tmd_in_metasub$microbiome_host), ncol = length(lev)) #create empty matrix with dataframe dimentions
char.var = bac_tmd_in_metasub$microbiome_host

for (i in 1:length(lev)) {
  matched_indices <- grep(lev[i], char.var, fixed = FALSE)# Find the indices matching the pattern
  result_bac_metasub[matched_indices, i] <- "Present" # Update values to "Present"
  empty_indices <- which(result_bac_metasub[, i] == "") # Replace empty strings with "Absent"
  result_bac_metasub[empty_indices, i] <- "Absent"
}

result_bac_metasub <- data.frame(result_bac_metasub, stringsAsFactors = TRUE)
colnames(result_bac_metasub) <- lev #Replace column names in dataframe for microbiome hosts
bac_tmd_in_metasub_total_host = cbind(bac_tmd_in_metasub,result_bac_metasub) #bind new table to principal metasub abundance table

data_heatmap=select(bac_tmd_in_metasub_total_host, species, Animal,
                    Human, Coral, Fungi, Plant, Sponge, Unknown)

melt_data <- melt(data_heatmap, id = c("species")) #Transform table for plotting
melt_data=melt_data[melt_data$value=="Present",]
rownames(data_heatmap)=data_heatmap$species
data_heatmap=data_heatmap[,-1]
data_heatmap=na.omit(data_heatmap)

heatmap_hosts=ggplot(melt_data, aes(variable,species,fill= value)) + 
  facet_grid(~variable, scales = "free")+
  geom_tile()+scale_fill_manual(values=c("black"))+theme_bw()+theme(axis.text.y = element_text(size=3, color="black"))

#do not run
#To create cladrogram for heatmap 
#write.table(row.names(data_heatmap), "metasub_hosts_taxa.txt")
#cladogram_metasub.py
#metasub_tree=read.tree("tree_hosts_bac_metasub.txt") #Import tree generated in phython
#metasub_tree$tip.label<-gsub("_"," ",metasub_tree$tip.label) #convert underscore to space in tree names
#gg_metasub_tree=ggtree(metasub_tree,layout = "rectangular",size=0.2, branch.length = "none") +  geom_tiplab(size=0,  offset = 0)+xlim_tree(50)+ylim(-50,650)
#gheatmap(gg_metasub_tree, data_heatmap, width = 1.5,high = 1,hjust = 1.1,colnames_angle=90, font.size = 6)+ scale_fill_manual(values = c("gray98","black"), na.value="gray95") + theme(legend.position = "top")





