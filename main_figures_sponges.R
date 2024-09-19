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


### Sponge microbiome ###

#Pathogenic Fungi
sponge_metadata=metadata_datasets[metadata_datasets$Project=="Sponge",] #Get sponge microbiome samples only
sponge_kraken=kraken2_datasets[kraken2_datasets$Sample_ID %in% sponge_metadata$Sample_ID ,] #Extract sample counts from main matrix
sponge_kraken_fungi=sponge_kraken[grep(c("Fungi"), sponge_kraken$taxa,ignore.case=T),]  #Select only Fungi 
sponge_kraken_fungi_sp=sponge_kraken_fungi[sponge_kraken_fungi$rank_code=="S",]  #Select only species within fingu
sponge_kraken_fungi_sp=merge(sponge_kraken_fungi_sp,sponge_metadata) #Merge with metadata

fungi_tmd_in_sponges=fungi_tmd[fungi_tmd$scientific_name %in% unique(sponge_kraken_fungi_sp$species),] #Get fungi present in sponges that are annotated in TMD
names(fungi_tmd_in_sponges)[1]="species"
fungi_tmd_in_sponges_metadata=fungi_tmd_in_sponges

fungi_patho_sponges=fungi_tmd_in_sponges_metadata[fungi_tmd_in_sponges_metadata$pathogen=="Yes",] #Get species in sponges classified as pathogens by TMD
fungi_patho_sponges <- fungi_patho_sponges %>%
  group_by(Sponge) %>%
  mutate(clade_reads_rel = 100*(clade_reads / sum(clade_reads))) #Create relative abundance tables

ggplot(fungi_patho_sponges, aes(x=reorder(species, -clade_reads_rel), y=clade_reads_rel, fill=Sponge)) +
  geom_bar(stat="identity")+facet_grid(~Sponge, scales = "free")+
  theme_pubr()+ scale_fill_brewer(palette = "Dark2")+
  xlab("Pathogen Fungi Species")+ ylab("Rel. Abundance")


#Microbiome sources for bacteria

sponge_kraken_bact=sponge_kraken[grep(c("Bacteria"), sponge_kraken$taxa,ignore.case=T),] #Get only bacteria in sponge microbiome
sponge_kraken_bact_sp=sponge_kraken_bact[sponge_kraken_bact$rank_code=="S",]  #Get only species of bacteria
sponge_kraken_bact_sp=merge(sponge_kraken_bact,sponge_metadata) #Merge counts with sample metadata
bac_tmd_in_sponges=bac_tmd[bac_tmd$species %in% unique(sponge_kraken_bact_sp$species),] #Bacteria species present in sponges that are annotated in TMD

sponge_kraken_bact_sp$tmd_hmp_status=if_else(sponge_kraken_bact_sp$species %in% unique(bac_tmd_in_sponges$species),"Annotated", "Not_annotated") #Classified bacteria species as annotated or not annotated in TMD
table(sponge_kraken_bact_sp$tmd_hmp_status) #Check annotation species counts
sponge_kraken_bact_sp_annt=sponge_kraken_bact_sp[sponge_kraken_bact_sp$tmd_hmp_status=="Annotated",] #Extract only annotated bacteria

merged_kraken_sponge_metadata=merge(sponge_kraken_bact_sp_annt,bac_tmd, by = "species") #Merge matrix abundance with TMD data.


sponge_bac_host=merged_kraken_sponge_metadata[grep(c("Host"), merged_kraken_sponge_metadata$microbiome),] #Get strings that have host microbiome
sponge_bac_host$custom="Host"
sponge_bac_water=merged_kraken_sponge_metadata[grep(c("Water"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have water microbiome
sponge_bac_water$custom="Water"
sponge_bac_soil=merged_kraken_sponge_metadata[grep(c("Soil"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have soil microbiome
sponge_bac_soil$custom="Soil"
sponge_bac_ext=merged_kraken_sponge_metadata[grep(c("Extreme"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have extreme microbiome
sponge_bac_ext$custom="Extreme"
sponge_bac_urb=merged_kraken_sponge_metadata[grep(c("Urban_Environment"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have urban microbiome
sponge_bac_urb$custom="Urban_Environment"
sponge_bac_food=merged_kraken_sponge_metadata[grep(c("Food"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have food microbiome
sponge_bac_food$custom="Food"
sponge_bac_air=merged_kraken_sponge_metadata[grep(c("Air"), merged_kraken_sponge_metadata$microbiome),]#Get strings that have air microbiome
sponge_bac_air$custom="Air"


df=rbind(sponge_bac_host,sponge_bac_water,sponge_bac_soil,sponge_bac_ext, sponge_bac_urb, sponge_bac_food, sponge_bac_air) #Bind the rows together
df$clade_reads=as.numeric(df$clade_reads)

df <- df %>%
  group_by(Sponge) %>%
  mutate(clade_reads_rel = 100*(clade_reads / sum(clade_reads))) #Normalize clade reads by sponge species

df <- df %>%
  group_by(Sponge,custom) %>%
  mutate(microbiome_percent = sum(clade_reads_rel)) #Group by type of microbiome and sponge species

df_select=unique(select(df,custom, Sponge, microbiome_percent)) #Select specific columns to plot

microbiome_source_sponge=ggplot(df_select, aes(custom, Sponge, size=microbiome_percent, color=Sponge, 
                      label=paste(round(microbiome_percent,1),"%", sep = "")))+
  geom_point(show.legend = F)+geom_text(color="black", size=6,vjust=-1)+
  xlab("Microbiome by TMD")+ylab("Sponge species")+
  labs(size="Percentage of bacteria", color="Microbiome by TMD")+
  scale_size(range = c(1, 15))


microbiome_source_sponge


#Sponge microbiome metabolism

names(merged_kraken_sponge_metadata)[18]="Sponge_sp" #Rename column for easier processing
lev=unique(unlist(strsplit(merged_kraken_sponge_metadata$metabolism, " "))) #Split metabolism column
result_metab <- matrix(data = "", nrow = length(merged_kraken_sponge_metadata$metabolism), ncol = length(lev)) #create empty matrix with data frame dimentions
char.var = merged_kraken_sponge_metadata$metabolism

for (i in 1:length(lev)) {
  matched_indices <- grep(lev[i], char.var, fixed = FALSE) # Find the indices matching the pattern
  result_metab[matched_indices, i] <- "Present" # Update values to "Present"
  empty_indices <- which(result_metab[, i] == "") # Replace empty strings with "Absent"
  result_metab[empty_indices, i] <- "Absent"
}

result_metab <- data.frame(result_metab, stringsAsFactors = TRUE)
colnames(result_metab) <- lev
merged_kraken_sponge_metadata$metabolism=NULL
merged_kraken_sponge_metadata_total = cbind(merged_kraken_sponge_metadata,result_metab) #Bind new metabolism split-table to principal sponge microbiome abundances table
metab_summarized=select(merged_kraken_sponge_metadata_total,Sponge_sp, Chemo, Heterotroph, Organo, Photo, Autotroph, Litho, Unknown) #Select the columns to plot

result <- lapply(metab_summarized[, -1], function(col) table(metab_summarized$Sponge_sp, col)) # Get frequencies of each column based on the first column
y <-melt(result, id = "col", measure = c("Absent", "Present")) #Restructure table
colnames(y)=c("Sponge_sp", "Presence", "Frequency", "Metabolism")

z=y[y$Presence=="Present",]
z <- z %>%
  group_by(Sponge_sp) %>%
  mutate(Relative_Abundance = 100*(Frequency / sum(Frequency))) #Convert frequencies (present) into relative abundances

z$Group=""
z$Group=ifelse(z$Metabolism=="Photo"|z$Metabolism=="Chemo", "Energy source","") #Group metabolism into primary nutritional groups
z$Group=ifelse(z$Metabolism=="Organo"|z$Metabolism=="Litho", "Electron source",z$Group)
z$Group=ifelse(z$Metabolism=="Heterotroph"|z$Metabolism=="Autotroph", "Carbon source",z$Group)
z$Group=ifelse(z$Metabolism=="Unknown", "Unknown",z$Group)


metabolism_microbiome_sponges=ggplot(z, aes(x=reorder(Metabolism,-Relative_Abundance),
                                            Relative_Abundance,fill=Metabolism))+
  geom_bar(stat = "identity")+
  facet_grid(~Group, scales = "free")+
  ylab("Percentage of bacteria")+xlab("")+
  scale_fill_brewer(palette = "Dark2")


metabolism_microbiome_sponges
