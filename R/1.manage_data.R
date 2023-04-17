## Script to access and process the data ##
## In this script there are the building of the phyloseq object and the table containing diversity values, time and envrionmental values
## Load libraries
library(phyloseq)   # Manipulation of metabarcoding data  
library(tidyverse)  # Plotting and data manipulation 
library(stringi)

## Load data
## Sequence
bact <- read.csv(unz("data/filtered_samples_tables/data_filtered_Bact02.zip", "data_filtered_Bact02.csv"))
coll <- read.csv(unz("data/filtered_samples_tables/data_filtered_Coll01.zip", "data_filtered_Coll01.csv"))
euka <- read.csv(unz("data/filtered_samples_tables/data_filtered_Euka02.zip", "data_filtered_Euka02.csv"))
fung <- read.csv(unz("data/filtered_samples_tables/data_filtered_Fung02.zip", "data_filtered_Fung02.csv"))
inse <- read.csv(unz("data/filtered_samples_tables/data_filtered_Inse01.zip", "data_filtered_Inse01.csv"))
olig <- read.csv(unz("data/filtered_samples_tables/data_filtered_Olig01.zip", "data_filtered_Olig01.csv"))
sper <- read.csv(unz("data/filtered_samples_tables/data_filtered_Sper01.zip", "data_filtered_Sper01.csv"))

# Samples
labels <- read.csv(file = "data/labels.csv")

## Sequence data ####
## Fungi
## Select taxonomy table
fung.tax <- fung[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
fung.tax <-  fung.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
fung.motu <- fung[,-c(2:25)]
## Manipulation of OTU table
fung.motu <- fung.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
fung.motu.relax <- fung.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
fung.phy.relax = phyloseq(otu_table(fung.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(fung.tax),
                          sample_data(fung.motu.relax[,1:3]))
fung.phy.relax <- prune_taxa(taxa_sums(fung.phy.relax) > 0, fung.phy.relax) 

## Plants
## Select taxonomy table
sper.tax <- sper[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
sper.tax <-  sper.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
sper.motu <- sper[,-c(2:25)]
## Manipulation of OTU table
sper.motu <- sper.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
sper.motu.relax <- sper.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
sper.phy.relax = phyloseq(otu_table(sper.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(sper.tax),
                          sample_data(sper.motu.relax[,1:3]))
sper.phy.relax <- prune_taxa(taxa_sums(sper.phy.relax) > 0, sper.phy.relax) 

## Bacteria
## Select taxonomy table
bact.tax <- bact[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
bact.tax <-  bact.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
bact.motu <- bact[,-c(2:25)]
## Manipulation of OTU table
bact.motu <- bact.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
bact.motu.relax <- bact.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
bact.phy.relax = phyloseq(otu_table(bact.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(bact.tax),
                          sample_data(bact.motu.relax[,1:3]))
bact.phy.relax <- prune_taxa(taxa_sums(bact.phy.relax) > 0, bact.phy.relax) 

## Fungi
## Select taxonomy table
coll.tax <- coll[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
coll.tax <-  coll.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
coll.motu <- coll[,-c(2:25)]
## Manipulation of OTU table
coll.motu <- coll.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
coll.motu.relax <- coll.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
coll.phy.relax = phyloseq(otu_table(coll.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(coll.tax),
                          sample_data(coll.motu.relax[,1:3]))
coll.phy.relax <- prune_taxa(taxa_sums(coll.phy.relax) > 0, coll.phy.relax) 

## Insects
## Select taxonomy table
inse.tax <- inse[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
inse.tax <-  inse.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
inse.motu <- inse[,-c(2:25)]
## Manipulation of OTU table
inse.motu <- inse.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
inse.motu.relax <- inse.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
inse.phy.relax = phyloseq(otu_table(inse.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(inse.tax),
                          sample_data(inse.motu.relax[,1:3]))
inse.phy.relax <- prune_taxa(taxa_sums(inse.phy.relax) > 0, inse.phy.relax) 

## Oligochetes
## Select taxonomy table
olig.tax <- olig[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
olig.tax <-  olig.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
olig.motu <- olig[,-c(2:25)]
## Manipulation of OTU table
olig.motu <- olig.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
olig.motu.relax <- olig.motu %>% 
  group_by(Glacier, Year, Spot) %>% # Group by spot (A, B, C, etc.)
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
olig.phy.relax = phyloseq(otu_table(olig.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(olig.tax),
                          sample_data(olig.motu.relax[,1:3]))
olig.phy.relax <- prune_taxa(taxa_sums(olig.phy.relax) > 0, olig.phy.relax) 

## Eukaryotes
## Select taxonomy table
euka.tax <- euka[c("id", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name", "scientific_name", "sequence")]
## Convert to character matrix
euka.tax <-  euka.tax %>% 
  column_to_rownames("id") %>% 
  as.matrix()

## Select OTU table
euka.motu <- euka[,-c(2:25)]
## Manipulation of OTU table
euka.motu <- euka.motu %>%
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "name") %>% 
  filter(!str_detect(name, "preba_")) %>% 
  filter(!str_detect(name, "miage_")) %>%
  filter(!str_detect(name, "_supr_")) %>%
  filter(!str_detect(name, "frajo_1500_")) %>%
  filter(!str_detect(name, "explo_...._.2")) %>%
  separate(col = name, into = c("Glacier", "Year", "Spot", "Replicate"), sep = "_", remove = F) %>%
  column_to_rownames(var = "name")

## Sum PCR replicates (relaxed stringency method cf. Mächler et al. 2021 10.1111/mec.15725)
euka.motu.relax <- euka.motu %>% 
  group_by(Glacier, Year, Spot) %>%
  dplyr::summarise(across(where(is.numeric), sum))

## Make phyloseq objects
euka.phy.relax = phyloseq(otu_table(euka.motu.relax[,-c(1:3)],
                                    taxa_are_rows = F),
                          tax_table(euka.tax),
                          sample_data(euka.motu.relax[,1:3]))
euka.phy.relax <- prune_taxa(taxa_sums(euka.phy.relax) > 0, euka.phy.relax) 

## Diversity ####
## bact
sample_names(bact) <- paste(sample_data(bact)$Glacier, sample_data(bact)$Year, sample_data(bact)$Spot, sep = "_")
sample_data(bact)$uniqPlot <- substr(sample_names(bact), start = 1, stop = 12) 
bact = merge_samples(bact, "uniqPlot")
bact.div <- estimate_richness(bact, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%  
  mutate(q1 = exp(Shannon))

## coll
sample_names(coll) <- paste(sample_data(coll)$Glacier, sample_data(coll)$Year, sample_data(coll)$Spot, sep = "_")
sample_data(coll)$uniqPlot <- substr(sample_names(coll), start = 1, stop = 12)
coll = merge_samples(coll, "uniqPlot")
coll.div <- estimate_richness(coll, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## euka.uni
sample_names(euka) <- paste(sample_data(euka)$Glacier, sample_data(euka)$Year, sample_data(euka)$Spot, sep = "_")
euka.uni = subset_taxa(euka, phylum_name=="Cercozoa"|phylum_name=="Ciliophora"|phylum_name=="Apicomplexa"|phylum_name=="Hemimastigophora"|phylum_name=="Perkinsozoa"|phylum_name=="Picozoa"|phylum_name=="Bacillariophyta"|phylum_name=="Cryptophyta"|class_name=="Dinophyceae"|class_name=="Eustigmatophyceae"|class_name=="Oomycetes"|class_name=="Xanthophyceae"|class_name=="Variosea"|class_name=="Synurophyceae"|class_name=="Chrysophyceae"|class_name=="Hyphochytriomycetes"|class_name=="Xanthophyceae"|class_name=="Synurophyceae"|class_name=="Glaucocystophyceae"|class_name=="Labyrinthulomycetes"|order_name=="Stemonitida"|order_name=="Plasmodiophorida"|order_name=="Eccrinales"|order_name=="Bicosoecida"|order_name=="Leptomyxida"|order_name=="Protosteliales"|order_name=="Arcellinida"|order_name=="Choanoflagellata"|order_name=="Vampyrellida"|order_name=="Longamoebia"|family_name=="Hartmannellidae"|family_name=="Apusomonadidae"|family_name=="Heterophryidae"|family_name=="Oikomonadaceae"|family_name=="Nucleariidae"|family_name=="Echinamoebidae"|family_name=="Pterocystidae"|family_name=="Raphidiophryidae"|family_name=="Vannellidae") # unicellular 
sample_data(euka.uni)$uniqPlot <- substr(sample_names(euka.uni), start = 1, stop = 12)
euka.uni = merge_samples(euka.uni, "uniqPlot")
euka.uni.div <- estimate_richness(euka.uni, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## euka.ani
sample_names(euka) <- paste(sample_data(euka)$Glacier, sample_data(euka)$Year, sample_data(euka)$Spot, sep = "_")
euka.ani = subset_taxa(euka, phylum_name=="Platyhelminthes"|phylum_name=="Nematoda"|phylum_name=="Tardigrada"|phylum_name=="Rotifera"|phylum_name=="Gastrotricha"|phylum_name=="Mollusca"|phylum_name=="Arthropoda") # animals
sample_data(euka.ani)$uniqPlot <- substr(sample_names(euka.ani), start = 1, stop = 12)
euka.ani = merge_samples(euka.ani, "uniqPlot")
euka.ani.div <- estimate_richness(euka.ani, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## fung
sample_names(fung) <- paste(sample_data(fung)$Glacier, sample_data(fung)$Year, sample_data(fung)$Spot, sep = "_")
sample_data(fung)$uniqPlot <- substr(sample_names(fung), start = 1, stop = 12)
fung = merge_samples(fung, "uniqPlot")
fung.div <- estimate_richness(fung, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## inse
sample_names(inse) <- paste(sample_data(inse)$Glacier, sample_data(inse)$Year, sample_data(inse)$Spot, sep = "_")
sample_data(inse)$uniqPlot <- substr(sample_names(inse), start = 1, stop = 12)
inse = merge_samples(inse, "uniqPlot")
inse.div <- estimate_richness(inse, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## olig
sample_names(olig) <- paste(sample_data(olig)$Glacier, sample_data(olig)$Year, sample_data(olig)$Spot, sep = "_")
sample_data(olig)$uniqPlot <- substr(sample_names(olig), start = 1, stop = 12)
olig = merge_samples(olig, "uniqPlot")
olig.div <- estimate_richness(olig, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>% 
  mutate(q1 = exp(Shannon))

## sper
sample_names(sper) <- paste(sample_data(sper)$Glacier, sample_data(sper)$Year, sample_data(sper)$Spot, sep = "_")
sample_data(sper)$uniqPlot <- substr(sample_names(sper), start = 1, stop = 12)
sper = merge_samples(sper, "uniqPlot")
sper.div <- estimate_richness(sper, measures = c("Observed", "Shannon")) %>% 
  rownames_to_column() %>% 
  separate(col = rowname, into = c("Glacier", "Year", "Plot"), sep = "_", remove = T) %>%  # Split name in different column
  mutate(q1 = exp(Shannon))

## Diversity dataset
div.table <- data.frame(Glacier = sper.div$Glacier, Year = sper.div$Year, Plot = sper.div$Plot, 
                        sper.q1 = sper.div$q1, bact.q1 = bact.div$q1, coll.q1 = coll.div$q1, olig.q1 = olig.div$q1, inse.q1 = inse.div$q1, euka.uni.q1 = euka.uni.div$q1, euka.ani.q1 = euka.ani.div$q1, fung.q1 = fung.div$q1,
                        sper.q0 = sper.div$Observed, bact.q0 = bact.div$Observed, coll.q0 = coll.div$Observed, olig.q0 = olig.div$Observed, inse.q0 = inse.div$Observed, euka.uni.q0 = euka.uni.div$Observed, euka.ani.q0 = euka.ani.div$Observed, fung.q0 = fung.div$Observed)

## Time since retreat ####
samples <- read.csv(file = "data/IC_sampled_points_may2021.csv", sep = ",")

## Remove unknown date
samples <- samples %>% filter(!is.na(date)) %>% filter(date != "") %>% filter(date != " ") %>% 
  filter(!str_detect(name, "explo_...._.2")) %>%
  select(name, glacier, dating, date)

length(unique(samples$glacier))# 50

## Replace uncorrect date
samples$date <- as.factor(samples$date) ## Factorize
levels(samples$date)[match("00/05/2015",levels(samples$date))] <- "01/05/2015"   

## Deal with different format type
a <- as.Date(samples$date, format="%d/%m/%Y") 
b <- as.Date(samples$date, format="%Y-%m-%d") 
a[is.na(a)] <- b[!is.na(b)]
samples$date_uniq <- a

## Extract only the year of sampling
samples$date_year <- as.numeric(format(samples$date_uniq, format = "%Y"))

## Caluclate time since glacier retreat
samples$time = samples$date_year - as.numeric(as.character(samples$dating)) 
samples$time.log <- log(samples$time)
samples$time.log.sc <- scale(samples$time.log)

## Create unique plot name
samples$uniqPlot <- substr(samples$name, start = 1, stop = 12)

samples.sel <- samples %>% 
  select(uniqPlot, glacier, time, time.log, time.log.sc) %>%
  group_by(uniqPlot) %>% 
  slice(1) %>% 
  drop_na()

## Environment ####
chem <- read.table('data/data_compete.txt', sep = '', header = T) 

chem.clean <- chem %>% 
  mutate(uniqPlot = substr(code, start = 1, stop = 12)) %>% 
  select(-lon, -lat, -code, -Glacier, -Year, -Plot, -time_retreat)

## Climate
clim <- read.csv("data/16.environmental variables+twi+slope+tpi_new.csv")

## Keep only the relevant glaciers and columns
clim.clean <- clim %>% 
  dplyr::filter(stringr::str_detect(name, "\\?", negate = TRUE)) %>% 
  dplyr::filter(stringr::str_detect(dating, "^S", negate = TRUE)) %>%
  mutate(uniqPlot = substr(name, start = 1, stop = 12)) %>% 
  select(-date, -glacier, -dating, -name) %>% 
  group_by(uniqPlot) %>% 
  summarise_all(mean) %>%
  drop_na()

## Joining tables ####
## Add time to diversity table
div.table$uniqPlot <- paste(div.table$Glacier, div.table$Year, div.table$Plot, sep = '_')
div.table <- div.table %>% 
  left_join(samples.sel, by = "uniqPlot")

## Add climatic and chemistry to diversity table
div.clim.chem <- div.table %>%
  left_join(chem.clean, by = 'uniqPlot') %>% 
  left_join(clim.clean, by = 'uniqPlot')
length(unique(div.clim.chem$Glacier))# 46

div.clim.chem$glacier[div.clim.chem$glacier == "Glacier de Gébroulaz"] <- "Gébroulaz"
div.clim.chem$glacier[div.clim.chem$glacier == "Glacier Blanc-Noir"] <- "Blanc-Noir"
div.clim.chem$glacier[div.clim.chem$glacier == "Glacier d'Estellette"] <- "Estellette"
div.clim.chem$glacier[div.clim.chem$glacier == "Glacier de Nantillons"] <- "Nantillons"
div.clim.chem$glacier[div.clim.chem$glacier == "Glacier des Pélerins"] <- "Pélerins"

