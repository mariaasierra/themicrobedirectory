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


### Human Microbiome Project###
#Virus host by TMD

hmp_metadata=metadata_datasets[metadata_datasets$Project=="Human_microbiome",] #Get only HMP data
hmd_kraken=kraken2_datasets[kraken2_datasets$Sample_ID %in% hmp_metadata$Sample_ID ,]

hmd_kraken_virus=hmd_kraken[grep(c("Viruses"), hmd_kraken$taxa,ignore.case=T),]  #Get only viruses
hmd_kraken_virus_sp=hmd_kraken_virus[hmd_kraken_virus$rank_code=="S",] #Get only species of viruses
hmd_kraken_virus_sp_metadata=merge(hmd_kraken_virus_sp,hmp_metadata) #Merge with sample metadata

virus_tmd_in_hmp=virus_tmd[virus_tmd$species %in% unique(hmd_kraken_virus_sp_metadata$species),] #Get virus species present in HMP that are also  annotated in TMD
hmd_kraken_virus_sp_metadata$tmd_hmp_status=if_else(hmd_kraken_virus_sp_metadata$species %in% virus_tmd$species,"Annotated", "Not_annotated")
hmd_kraken_virus_sp_metadata_annt=hmd_kraken_virus_sp_metadata[hmd_kraken_virus_sp_metadata$tmd_hmp_status=="Annotated",]

hmd_kraken_virus_sp_metadata_annt=merge(hmd_kraken_virus_sp_metadata_annt,virus_tmd, by = "species") #Merge with sample metadata

hmp_virus_hosts=select(hmd_kraken_virus_sp_metadata_annt,species, Body_site, pathogen_host) #Subset data table

hmp_virus_hosts_sep=separate_rows(hmp_virus_hosts, pathogen_host) 
hmp_virus_hosts_sep= subset(hmp_virus_hosts_sep, pathogen_host != "") #Rename body sites
site_names <- c(
  `right retroauricular crease` = "Right Retroauricular Crease",
  `gingiva` = "Gingiva",
  `gastrointestinal tract` = "Gastrointestinal Tract",
  `throat` = "Throat"
)


hmp_virus_plot=ggplot(hmp_virus_hosts_sep, aes(x=species, y=pathogen_host, col=pathogen_host))+
  geom_point(size=6)+
  facet_wrap(~Body_site,labeller = as_labeller(site_names,label_wrap_gen(multi_line = TRUE)))+
  scale_color_brewer(palette = "Set1", direction = -1)+
  theme_bw()+ 
  xlab("Viral species")+ labs(col="Host by TMD") +ggtitle("Body sites by HMP")


hmp_virus_plot 


#Proportion of Gram-positive and Gram-negative 
hmp_metadata=metadata_datasets[metadata_datasets$Project=="Human_microbiome",] #Get only HMP data
hmd_kraken=kraken2_datasets[kraken2_datasets$Sample_ID %in% hmp_metadata$Sample_ID ,] 

hmd_kraken_bac=hmd_kraken[grep(c("Bacteria"), hmd_kraken$taxa,ignore.case=T),] #Get only bacteria
hmd_kraken_bac_sp=hmd_kraken_bac[hmd_kraken_bac$rank_code=="S",] #Get only species
hmd_kraken_bac_sp=merge(hmd_kraken_bac_sp,hmp_metadata) #Merge list with metadata


hmd_kraken_bac_sp$tmd_hmp_status=if_else(tolower(hmd_kraken_bac_sp$species) %in% tolower(bac_tmd$scientific_name) ,"Annotated", "Not_annotated") #Classify TMD bacteria in annotated and not annotated species if present in TMD data
hmd_kraken_bac_sp_annt=hmd_kraken_bac_sp[hmd_kraken_bac_sp$tmd_hmp_status=="Annotated",] #Get only annotated bacteria species
names(bac_tmd)[1]="species"

merged_kraken_tmd_metadata=hmd_kraken_bac_sp_annt #already merged
merged_kraken_tmd_metadata$bac_gram_stain=if_else(merged_kraken_tmd_metadata$bac_gram_stain=="", "Unknown",merged_kraken_tmd_metadata$bac_gram_stain) #Fill any gaps of not annotated bacteria with "unknown"

gram_stain_hmp_body=ggpie(data = merged_kraken_tmd_metadata,x=bac_gram_stain, offset=1.1,border.color="white", border.width=0.5,
      by=Body_site, nrow=3,label.size=5,legend=T)+scale_fill_brewer(palette = "Dark2")

gram_stain_hmp_body

#Oxygen requirement of bacteria present per body site
lev=unique(unlist(strsplit(merged_kraken_tmd_metadata$bac_oxygen_use, " "))) #Split oxygen-use column into separate columns
result_oxy_use <- matrix(data = "", nrow = length(merged_kraken_tmd_metadata$bac_oxygen_use), ncol = length(lev)) #Create empty matrix with dataframe dimentions
char.var = merged_kraken_tmd_metadata$bac_oxygen_use

for (i in 1:length(lev)) {
  matched_indices <- grep(lev[i], char.var, fixed = FALSE)# Find the indices matching the pattern
  result_oxy_use[matched_indices, i] <- "Present" # Update values to "Present"
  empty_indices <- which(result_oxy_use[, i] == "")  # Replace empty strings with "Absent"
  result_oxy_use[empty_indices, i] <- "Absent"
}

result_oxy_use <- data.frame(result_oxy_use, stringsAsFactors = TRUE) #New data frame with oxygen use requirements in columns
colnames(result_oxy_use) <- lev
merged_kraken_tmd_metadata$bac_oxygen_use=NULL
merged_kraken_tmd_metadata_total = cbind(merged_kraken_tmd_metadata,result_oxy_use) #Merge new oxygen table to principal HMP bacteria table

oxy_use_summarized=select(merged_kraken_tmd_metadata_total, Body_site, Microaerophile, 
                          Aerobes, Anaerobes, Obligate_anaerobes, Facultative_anaerobes, 
                          Obligate_aerobes, Aerotolerant_anaerobes) #Select only desired columns to plot 

result <- lapply(oxy_use_summarized[, -1], function(col) table(oxy_use_summarized$Body_site, col)) # Get frequencies of each column based on the first column
y <-melt(result, id = "col", measure = c("Absent", "Present"))
colnames(y)=c("Body_site", "Presence", "Frequency", "Oxygen_use")

dataframe_oxy=y[y$Presence=="Present",]
dataframe_oxy <- dataframe_oxy %>%
  group_by(Body_site) %>%
  mutate(Relative_Abundance = 100*(Frequency / sum(Frequency))) #Get frequencies table per body site

dataframe_oxy$Oxygen_use <- sub("_", " ", dataframe_oxy$Oxygen_use) #For aesthetics, replace underscore for spaces
oxygen_names=unique(dataframe_oxy$Oxygen_use)

ggplot(dataframe_oxy, aes(x=Body_site, y=Relative_Abundance, fill=Body_site))+ geom_bar(stat="identity")+ 
  facet_grid(~Oxygen_use, scales = "free")+theme_bw()+ scale_fill_brewer(palette = "Set1", labels=site_names)+
  labs(fill="Body site")+xlab("Body site")+ ylab("Rel. Abundance")

dataframe_oxy




