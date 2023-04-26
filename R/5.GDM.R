## GDM
## In this script there is the code for the Generalized Dissimilarity Modeling

## Load packages and data ####
library(phyloseq)
library(tidyverse)
library(gdm)

beta.tot=read.csv("data/beta.biotic.plot.csv",header=T,stringsAsFactors=F)
all_data=read.delim("data/div.clim.chem.csv", sep=";")
div=all_data[!is.na(all_data$lg_n),] # remove rows with NA

## Format table
env <- env %>% 
  rename(Spot = Plot) %>% 
  mutate(Plot = substr(uniqPlot, start = 1, stop = 12)) %>% 
  mutate(Site = substr(uniqPlot, start = 1, stop = 10)) %>% 
  drop_na()

## Transformation
env$t.s <- scale(env$meanT_new)
env$twi.s <- scale(log(env$twi))
env$ndvi.s <- scale(log(env$ndvi))

env.plot <- env %>%
  group_by(Plot) %>%
  dplyr::select(Plot, lon, lat, time.s, t.s, twi.s, ndvi.s, ph.s, n.s, p.s) %>% 
  summarise_all(mean) %>% 
  mutate(Glacier = substr(Plot, start = 1, stop = 5)) 

## Keep only plot with environmental variables measured and MOTU present (from SM script for SEM 2.ecodist-MRM_permute sem_soil.R)
div$all.ani.q0=div$coll.q0+div$olig.q0+div$inse.q0+div$euka.ani.q0
zero_plus=apply(div[,which(names(div) %in% c("sper.q0","bact.q0","fung.q0","all.ani.q0"))],
                MARGIN=1,FUN=function(x){all(x>0)})
tokeep=div$uniqPlot[which(zero_plus==T)]
env.plot=env.plot[which(env.plot$Plot %in% tokeep),] # 563

## Community data from 4 groups: all the animals, plants, bacteria, fungi
## Plants
sper <- readRDS("data/phyloseq/sper.phy.relax.rds")
sample_names(sper) <- paste(sample_data(sper)$Glacier, sample_data(sper)$Year, sample_data(sper)$Spot, sep = "_")
sample_data(sper)$unikPlot <- substr(sample_names(sper), start = 1, stop = 12) 
sper = merge_samples(sper, "unikPlot")
sper <- subset_samples(sper, sample_names(sper) %in% env.plot$Plot)
sper.mat = as.data.frame(as(otu_table(sper), "matrix"))
sper.mat$Plot <- rownames(sper.mat)
gdmTab.sper <- formatsitepair(bioData=sper.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                              predData=env.plot, siteColumn="Plot", dist = "jaccard") 
gdmTab.sper <- gdmTab.sper[-which(gdmTab.sper$s1.Glacier != gdmTab.sper$s2.Glacier),]
gdmTab.sper = gdmTab.sper[,!(names(gdmTab.sper) %in% c("s1.Glacier","s2.Glacier"))]
gdm.sper <- gdm(data=gdmTab.sper, geo=T)
gdm.sper.signif <- gdm.varImp(spTable = gdmTab.sper, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)
imp.sper <- data.frame(deviance = gdm.sper.signif[[2]][,1], pvalue = gdm.sper.signif[[3]][,1], 
                       marker = rep('Vascular plants',length(gdm.sper.signif[[2]][,1])))
imp.sper$variables <- rownames(imp.sper)

## Plant community
sper.pcoa <- ordinate(sper, method="PCoA", distance = "jaccard")
sper.pcoa.load <- data.frame(compo.axis1 = sper.pcoa$vectors[,1], compo.axis2 = sper.pcoa$vectors[,2])
sper.pcoa.load$Plot <- rownames(sper.pcoa.load)
env.plant.plot <- sper.pcoa.load %>% 
  left_join(env.plot, by = "Plot") %>% 
  drop_na()

## Fungi
fung <- readRDS("data/phyloseq/fung.phy.relax.rds")
sample_names(fung) <- paste(sample_data(fung)$Glacier, sample_data(fung)$Year, sample_data(fung)$Spot, sep = "_")
sample_data(fung)$unikPlot <- substr(sample_names(fung), start = 1, stop = 12)
fung = merge_samples(fung, "unikPlot") 
fung <- subset_samples(fung, sample_names(fung) %in% env.plant.plot$Plot)
fung.mat = as.data.frame(as(otu_table(fung), "matrix"))
fung.mat$Plot <- rownames(fung.mat)
gdmTab.fung <- formatsitepair(bioData=fung.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                              predData=env.plant.plot, siteColumn="Plot", dist = "jaccard", abundance = F)
gdmTab.fung <- gdmTab.fung[-which(gdmTab.fung$s1.Glacier != gdmTab.fung$s2.Glacier),]
gdmTab.fung = gdmTab.fung[,!(names(gdmTab.fung) %in% c("s1.Glacier","s2.Glacier"))]
gdm.fung <- gdm(data=gdmTab.fung, geo=T)
gdm.fung.signif <- gdm.varImp(spTable = gdmTab.fung, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)
imp.fung <- data.frame(deviance = gdm.fung.signif[[2]][,1], pvalue = gdm.fung.signif[[3]][,1], 
                       marker = rep('Fungi',length(gdm.fung.signif[[2]][,1])))
imp.fung$variables <- rownames(imp.fung)

## Bacteria
bact <- readRDS("data/phyloseq/bact.phy.relax.rds")
sample_names(bact) <- paste(sample_data(bact)$Glacier, sample_data(bact)$Year, sample_data(bact)$Spot, sep = "_")
sample_data(bact)$unikPlot <- substr(sample_names(bact), start = 1, stop = 12)
bact = merge_samples(bact, "unikPlot") # samples
bact <- subset_samples(bact, sample_names(bact) %in% env.plant.plot$Plot) 
bact.mat = as.data.frame(as(otu_table(bact), "matrix"))
bact.mat$Plot <- rownames(bact.mat)
gdmTab.bact <- formatsitepair(bioData=bact.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                              predData=env.plant.plot, siteColumn="Plot", dist = "jaccard") 
gdmTab.bact <- gdmTab.bact[-which(gdmTab.bact$s1.Glacier != gdmTab.bact$s2.Glacier),]
gdmTab.bact = gdmTab.bact[,!(names(gdmTab.bact) %in% c("s1.Glacier","s2.Glacier"))]
gdm.bact <- gdm(data=gdmTab.bact, geo=T)
gdm.bact.signif <- gdm.varImp(spTable = gdmTab.bact, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)
imp.bact <- data.frame(deviance = gdm.bact.signif[[2]][,1], pvalue = gdm.bact.signif[[3]][,1], 
                       marker = rep('Bacteria',length(gdm.bact.signif[[2]][,1])))
imp.bact$variables <- rownames(imp.bact)

## All the animals
## coll
coll <- readRDS("data/phyloseq/coll.phy.relax.rds")
sample_names(coll) <- paste(sample_data(coll)$Glacier, sample_data(coll)$Year, sample_data(coll)$Spot, sep = "_")
sample_data(coll)$unikPlot <- substr(sample_names(coll), start = 1, stop = 12)
coll = merge_samples(coll, "unikPlot")
coll <- subset_samples(coll, sample_names(coll) %in% env.plant.plot$Plot) 
coll.mat = as.data.frame(as(otu_table(coll), "matrix"))
## euka.ani
euka <- readRDS("data/phyloseq/euka.phy.relax.rds")
euka.animals = subset_taxa(euka, phylum_name=="Platyhelminthes"|phylum_name=="Nematoda"|phylum_name=="Tardigrada"|phylum_name=="Rotifera"|phylum_name=="Gastrotricha"|phylum_name=="Mollusca"|phylum_name=="Arthropoda") # animals
sample_names(euka.animals) <- paste(sample_data(euka.animals)$Glacier, sample_data(euka.animals)$Year, sample_data(euka.animals)$Spot, sep = "_")
sample_data(euka.animals)$unikPlot <- substr(sample_names(euka.animals), start = 1, stop = 12) 
euka.animals = merge_samples(euka.animals, "unikPlot")
euka.animals <- subset_samples(euka.animals, sample_names(euka.animals) %in% env.plant.plot$Plot) 
euka.animals.mat = as.data.frame(as(otu_table(euka.animals), "matrix"))
## inse
inse <- readRDS("data/phyloseq/inse.phy.relax.rds")
sample_names(inse) <- paste(sample_data(inse)$Glacier, sample_data(inse)$Year, sample_data(inse)$Spot, sep = "_")
sample_data(inse)$unikPlot <- substr(sample_names(inse), start = 1, stop = 12) 
inse = merge_samples(inse, "unikPlot") 
inse <- subset_samples(inse, sample_names(inse) %in% env.plant.plot$Plot) 
inse.mat = as.data.frame(as(otu_table(inse), "matrix"))
## olig
olig <- readRDS("data/phyloseq/olig.phy.relax.rds")
sample_names(olig) <- paste(sample_data(olig)$Glacier, sample_data(olig)$Year, sample_data(olig)$Spot, sep = "_")
sample_data(olig)$unikPlot <- substr(sample_names(olig), start = 1, stop = 12)
olig = merge_samples(olig, "unikPlot")
olig <- subset_samples(olig, sample_names(olig) %in% env.plant.plot$Plot) 
olig.mat = as.data.frame(as(otu_table(olig), "matrix"))

## Merge all animals matrices
all.ani.mat <- cbind(coll.mat, euka.animals.mat, inse.mat, olig.mat)
all.ani.mat <- all.ani.mat[rowSums(all.ani.mat)>0,] #1177
all.ani.mat$Plot <- rownames(all.ani.mat)

## GDM for all animals
gdmTab.all.ani <- formatsitepair(bioData=all.ani.mat, bioFormat=1, XColumn="lon", YColumn="lat",
                                 predData=env.plant.plot, siteColumn="Plot", dist = "jaccard")
gdmTab.all.ani <- gdmTab.all.ani[-which(gdmTab.all.ani$s1.Glacier != gdmTab.all.ani$s2.Glacier),]
gdmTab.all.ani = gdmTab.all.ani[,!(names(gdmTab.all.ani) %in% c("s1.Glacier","s2.Glacier"))]
gdm.all.ani <- gdm(data=gdmTab.all.ani, geo=T)
gdm.all.ani.signif <- gdm.varImp(spTable = gdmTab.all.ani, geo=T, fullModelOnly = TRUE, nPerm = 1000, parallel = TRUE)
imp.all.ani <- data.frame(deviance = gdm.all.ani.signif[[2]][,1], pvalue = gdm.all.ani.signif[[3]][,1], 
                          marker = rep('Animals',length(gdm.all.ani.signif[[2]][,1])))
imp.all.ani$variables <- rownames(imp.all.ani)

## Merge the results for the 4 groups
imp.4groups <- rbind(imp.bact, imp.all.ani, imp.fung, imp.sper)

