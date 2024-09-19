setwd("~/NCBI/NCBI_microbes_lists")
library(stringr)

#ARCHAEA
archaea=read.csv("sci_names_archaea.csv", stringsAsFactors =F)
head(archaea)


# Create List of Exclusion words
todelete <- c('Candidatus', 'candidate','environmental', 'unclassified','unidentified', 'group', 'Group','OTU', 'JGI', 'SCGC', 'JdFR', 'archaeon', 'DSM', 'culture', '-','_', 
              'endosymbiont', 'symbiont', 'clone', 'enrichment', '/', 'symbiosum', 'strain','str.', 'cluster', 'Type',
              'sp.', '1', '2','3', '4','5','6', '7', '8', '9', '0')
archaea_wo_uncl=archaea[ grep(paste(todelete, collapse = "|"), archaea$Species, invert = TRUE) , ]
head(archaea_wo_uncl)
archaea_wo_uncl=as.data.frame(archaea_wo_uncl)
archaea_unique=unique(archaea_wo_uncl)


#DIATOMS
diatoms=read.csv("~/NCBI/sci_names_diatoms.csv", stringsAsFactors = F)
# Create List of Exclusion words
todelete.diatoms =c('cf.', 'STRAIN', 'unclassified', 'unidentified','uncultured', 'lineage', 'sp.', 'environmental', 'samples', '-', '_','endosymbiont', 'symbiont',
                    '1', '2','3', '4','5','6', '7', '8', '9', '0')
diatom_wo_uncl=diatoms[ grep(paste(todelete.diatoms, collapse = "|"), diatoms$Species, invert = TRUE) , ]
diatom_wo_uncl=as.data.frame(diatom_wo_uncl)
diatom_unique=unique(diatom_wo_uncl)

#BACTERIA
bacteria=read.csv("~/NCBI/sci_names_bacteria.csv", stringsAsFactors = F)
# Create List of Exclusion words
todelete.bac <- c('Candidatus', 'candidate','environmental', 'unclassified','unidentified', 'group', 'Group','OTU', 'culture', '-','_', 
              'endosymbiont', 'symbiont', 'clone', 'enrichment', '/', 'symbiosum', 'strain','str.', 'cluster', 'Type', 'bacterium',
              'bacteria','eubacteria', 'ensemble',
              'sp.', '1', '2','3', '4','5','6', '7', '8', '9', '0')
bacteria_wo_uncl=bacteria[ grep(paste(todelete.bac, collapse = "|"), bacteria$Species, invert = TRUE) , ]
bacteria_wo_uncl=as.data.frame(bacteria_wo_uncl)
bacteria_unique=unique(bacteria_wo_uncl)


#FUNGI
fungi=read.csv("~/NCBI/sci_names_fungi.csv", stringsAsFactors = F)
# Create List of Exclusion words
todelete.fungi <- c('Candidatus', 'candidate','environmental', 'unclassified','unidentified', 'group', 'Group','OTU', 'culture', '-','_', 
                  'endosymbiont', 'symbiont', 'clone', 'enrichment', '/', 'symbiosum', 'strain','str.', 'cluster', 'Type',
                   '1', '2','3', '4','5','6', '7', '8', '9', '0', 'Fungi','fungi')
fungi_wo_uncl=fungi[ grep(paste(todelete.fungi, collapse = "|"), fungi$Species, invert = TRUE) , ]
macrofungi <- c('Chlorociboria','Paxina', 'Helvella', 'Morchella', 'Kalaharituber', 'Peziza', 'Pseudohelotium', 'Terfezia',
                'Anthracobia', 'Isaria', 'Scutellinia', 'Rhizina', 'Phillipsia', 'Tuber', 'Choiromyces',
                'Cordyceps', 'Daldinia', 'Poronia', 'Xylaria', 'Penzigia', 'Agaricus', 'Arachnion', 
                'Battarrea', 'Battarreoides', 'Bovista', 'Calvatia', 'Chlamydopus', 'Chlorophyllum',
                'Coniolepiota', 'Coprinellus', 'Coprinopsis', 'Coprinus', 'Crucibulum', 'Disciseda',
                'Gyrophragmium', 'Langermannia', 'Lepiota', 'Leucoagaricus', 'Lycoperdon', 'Macrolepiota',
                'Montagnea', 'Montagnites', 'Mycenastrum', 'Parasola', 'Polyplocium', 'Psalliota', 'Secotium', 
                'Tulostoma', 'Xanthagaricus', 'Amanita', 'Amanita', 'Saproamanita', 'Bolbitius', 'Conocybe',
                'Galeropsis', 'Pluteolus', 'Broomeia', 'Clavaria', 'Clavaria', 'Clavulinopsis', 'Mucronella', 
                'Cortinarius', 'Locellina', 'Chondrostereum', 'Cyphella', 'Claudopus', 'Clitopilus', 'Entoloma',
                'Fistulina','Laccaria', 'Hygrocybe', 'Hygrophorus', 'Anellaria', 'Panaeolina', 'Panaeolus',
                'Astrosporina', 'Crepidotus', 'Inocybe', 'Phaeoglabrotricha', 'Phaeosolenia', 'Lyophyllum',
                'Podabrella', 'Termitomyces', 'Calyptella','Marasmius', 'Solenia', 'Cruentomycena',
                'Favolaschia', 'Mycena', 'Flagelloscypha', 'Lachnella', 'Cyathus', 'Anthracophyllum', 'Gymnopus',
                'Marasmiellus', 'Omphalotus', 'Phellorinia', 'Armillaria', 'Armillariella', 'Cyptotrama',
                'Oudemansiella', 'Hymenopellis', 'Physalacria', 'Xerula', 'Pleurotus', 'Pluteus', 'Volvariella',
                'Podaxis', 'Ozonium', 'Psathyrella', 'Pterula', 'Schizophyllum','Sebacina', 'Agrocybe',
                'Deconica', 'Flammula', 'Galera', 'Gymnopilus', 'Hebeloma', 'Hymenogaster', 'Hypholoma', 'Kuehneromyces',
                'Leratiomyces', 'Naucoria', 'Pholiota', 'Psilocybe', 'Stropharia', 'Tubaria', 'Amparoina', 'Cellypha',
                'Clitocybe', 'Collybia', 'Lepista', 'Macrocybe', 'Melanoleuca', 'Omphalia', 'Tricholoma', 'Tricholomopsis',
                'Tricholosporum', 'Trogia', 'Auricularia', 'Eichleriella', 'Exidia', 'Heterochaete', 'Aporpium', 'Aureoboletus',
                'Boletus' ,'Buchwaldoboletus', 'Chalciporus', 'Imleria', 'Leccinum', 'Octaviania','Xerocomellus',
                'Phlebopus', 'Coniophora','Gyrodontium', 'Gyroporus', 'Melanogaster', 'Paxillus', 'Rhizopogon', 'Pisolithus',
                'Scleroderma', 'Serpula', 'Suillus', 'Cantharellus', 'Pellicularia', 'Clavulina', 'Corticium', 'Cytidia',
                'Dendrothele', 'Laetiporus', 'Tretopileus', 'Geasteropsis', 'Geastrum', 'Geastrum', 'Anthurus', 'Aseroe',
                'Blumenavia', 'Clathrella', 'Clathrus', 'Ileodictyon', 'Itajahya', 'Jaczewskia', 'Kalchbrennera', 
                'Kalchbrennera', 'Mutinus', 'Phallus', 'Daedalea', 'Fomitopsis', 'Gloeocystidium', 'Phaeolus',
                'Rhodofomitopsis', 'Amauroderma', 'Ganoderma', 'Gloeophyllum', 'Ramaria', 'Clavariadelphus',
                'Coltricia' ,'Fomitoparia', 'Fuscoporia', 'Hydnum','Hymenochaete','Phellinus',
                'Polystictus','Trichaptum', 'Cotylidia', 'Grandinia', 'Heterochaete','Oxyporus','Crustodontia',
                'Acia','Aegerita','Bjerkandera', 'Cymatoderma', 'Gloeoporus', 'Irpex', 'Laschia', 'Merulius', 'Mycoleptodon',
                'Odontia', 'Phlebia', 'Podoscypha', 'Pseudolagarobasidium', 'Abortiporus', 'Coriolopsis', 'Coriolus',
                'Daedaleopsis', 'Favolus','Fomes', 'Funalia', 'Grammothele', 'Heliocybe','Hexagonia','Lentinus', 
                'Lenzites', 'Lenzites', 'Lenzites', 'Lignosus', 'Microporus', 'Nigroporus', 'Neolentinus', 'Panus',
                'Perenniporia', 'Picipes', 'Phellinus', 'Polyporus', 'Pycnoporus', 'Trametes' ,'Lentinellus', 'Dentipellicula',
                'Laxitextum', 'Asterostroma','Dichostereum', 'Lachnocladium', 'Peniophora', 'Lactarius', 'Lactifluus', 'Russula',
                'Aleurodiscus', 'Stereum', 'Hypochnus', 'Thelephora', 'Arrhytidia', 'Calocera', 'Dacrymyces', 'Dacryopinax',
                'Femsjonia', 'Naematoloma', 'Tremella', 'Phaeotremella', 'Sirobasidium', 'Pilobolus', 'Echinostelium', 'Cribraria',
                'Dictydiaethalium','Licea', 'Lycogala', 'Reticularia', 'Tubifera', 'Diachea', 'Diderma', 'Didymium', 'Mucilago',
                'Badhamia', 'Badhamiopsis', 'Craterium', 'Fuligo', 'Leocarpus', 'Physarella', 'Physarum', 'Willkommlangea',
                'Amaurochaete', 'Comatricha', 'Enerthenema', 'Lamproderma', 'Stemonaria', 'Stemonitis', 'Stemonitopsis', 'Calomyxa',
                'Arcyria', 'Hemitrichia', 'Metatrichia', 'Oligonema', 'Perichaena', 'Trichia','Ceratiomyxa','Ceratium'
                )
fungi_wo_macrofungi=fungi_wo_uncl[ grep(paste(macrofungi, collapse = "|"), fungi_wo_uncl$fungi_wo_uncl, invert = TRUE) , ]
fungi_wo_macrofungi=as.data.frame(fungi_wo_macrofungi)
fungi_wo_macrofungi_uniq=unique(fungi_wo_macrofungi)

#VIRUS
virus=read.csv("~/NCBI/sci_names_virus.csv", stringsAsFactors = F)
# Create List of Exclusion words
todelete.virus <- c('Candidatus', 'candidate','environmental', 'unclassified','unidentified', 'group', 'Group','OTU', 'culture',  
                    'endosymbiont', 'symbiont', 'clone', 'enrichment', '/', 'symbiosum', 'cluster', 'Type', 'phage', 'Phage', '[[]', '[]]]', 'isolate', 'particle', 
                    'ISOLATE', 'Isolate', 'plaque', 'PLAQUE', 'Plaque', 'PREDICT', 'SP.', 'sp.', ':', 'noda-like', 'AFVG')
virus_wo_uncl=virus[ grep(paste(todelete.virus, collapse = "|"), virus$Species, invert = TRUE) , ]
virus_wo_uncl=as.data.frame(virus_wo_uncl)
virus_unique=unique(virus_wo_uncl)

#ALGAE

chloro=read.csv("~/NCBI/sci_names_chlorophyta.txt", stringsAsFactors = F, sep = "\t")
strepto=read.csv("~/NCBI/sci_names_streptophyta.txt", stringsAsFactors = F, sep = "\t")
rodo=read.csv("~/NCBI/sci_names_rhodophyta.txt", stringsAsFactors = F, sep = "\t")
glauco=read.csv("~/NCBI/sci_names_glauco.txt", stringsAsFactors = F, sep = "\t")

unique(chloro$Family) 
chloro_1=subset(chloro, !is.na(Family))
unique(strepto$Family)
# Create List of Exclusion words
strepto_1=strepto[strepto$Family !="Leiosporocerotopsida" &
                  strepto$Family !="Sphagnopsida" &
                  strepto$Family !="Ginkgoopsida" &
                  strepto$Family !="Andreaeobryopsida" &
                  strepto$Family !="Takakiopsida" &
                  strepto$Family !="Gnetopsida" &
                  strepto$Family !="Andreaeopsida" &
                  strepto$Family !="Haplomitriopsida" &
                  strepto$Family !="Pinopsida" &
                  strepto$Family !="Bryopsida" &
                  strepto$Family !="Jungermanniopsida" &
                  strepto$Family !="Magnoliopsida" &
                    strepto$Family !="Charophyceae" &
                    strepto$Family !="Oedipodiopsida" &
                    strepto$Family !="Marchantiopsida" &
                    strepto$Family !="Lycopodiopsida" &
                    strepto$Family !="Polytrichopsida" &
                    strepto$Family !="Polypodiopsida" &
                    strepto$Family !="Anthocerotopsida" &
                    strepto$Family !="Tetraphidopsida" &
                    strepto$Family !="Cycadopsida" ,]
strepto_1=subset(strepto_1, !is.na(Family))
unique(rodo$Family)
rodo_1=rodo[rodo$Family !="Florideophyceae" ,]
rodo_1=subset(rodo_1, !is.na(Family))
unique(glauco$Family)


# Merge all dataframes For full (Outer) join
all_algae=list(chloro_1, strepto_1, rodo_1, glauco)
multi_full <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  all_algae
)

# Create List of Exclusion words
todelete.algae <- c('environmental', 'unclassified','unidentified', 'group', 'Group', 'culture',  
                    'endosymbiont', 'symbiont', 'clone', 'enrichment', '/', 'isolate', 'particle', 
                    'ISOLATE', 'Isolate', 'plaque', 'PLAQUE', 'Plaque', 'PREDICT', 'SP.', 'sp.', ':')
algae_wo_uncl=multi_full[ grep(paste(todelete.algae, collapse = "|"), multi_full$Species, invert = TRUE) , ]
algae_wo_uncl=as.data.frame(algae_wo_uncl)
algae_unique=unique(algae_wo_uncl)
