## GLMM
## In this script there is the code for the generalized linear mixed-effect models for the biodiversity and soil features

## Load packages and data ####
library(phyloseq)
library(tidyverse)
library(brms); options(mc.cores = parallel::detectCores())
library(ncf)
library(gdm)

div.table <- read.csv("data/div.clim.chem.csv") #from manage.clean.R

## Bact02
bact.q1.logn <- brm(formula = bact.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Coll01
coll.q1.logn <- brm(formula = coll.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.uni
euka.uni.q1.logn <- brm(formula = euka.uni.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                     data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.ani
euka.ani.q1.logn <- brm(formula = euka.ani.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                         data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Fung02
fung.q1.logn <- brm(formula = fung.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Inse01
inse.q1.logn <- brm(formula = inse.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Olig01
olig.q1.logn <- brm(formula = olig.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Sper01
sper.q1.logn <- brm(formula = sper.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## R2
r2_bayes(bact.rich.nb)
r2_bayes(coll.rich.nb)
r2_bayes(euka.uni.rich.nb)
r2_bayes(euka.animals.rich.nb)
r2_bayes(fung.rich.nb)
r2_bayes(inse.rich.nb)
r2_bayes(olig.rich.nb)
r2_bayes(sper.rich.nb)

## Plotting
## Spatial autocorrelation
par(mfrow=c(2,4))
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(bact.q1.logn))
plot(cres, ylim = c(-1,1), main = "Bacteria") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(fung.q1.logn))
plot(cres, ylim = c(-1,1), main = "Fungi") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(euka.uni.q1.logn))
plot(cres, ylim = c(-1,1), main = "Protists") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(sper.q1.logn))
plot(cres, ylim = c(-1,1), main = "Vascular plants") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(inse.q1.logn))
plot(cres, ylim = c(-1,1), main = "Insects") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(coll.q1.logn))
plot(cres, ylim = c(-1,1), main = "Springtails") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(olig.q1.logn))
plot(cres, ylim = c(-1,1), main = "Earthworms") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(euka.ani.q1.logn))
plot(cres, ylim = c(-1,1), main = "Other animals") # No sign spatial autocorrelation

## GLMMs for soil data ####
## Chemistry only
chem.table <- div.table %>% 
  drop_na() %>% 
  mutate(p = p/1000000)

## Nitrogen
n.hg <- brm(formula = n ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table, family=hurdle_gamma, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Carbon
c.logn <- brm(formula = c ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
              data = chem.table[!(chem.table$Glacier=="agola"),], family=lognormal, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Phosphorus
p.logn <- brm(formula = p ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
              data = chem.table, family=lognormal, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## pH
ph.g <- brm(formula = ph ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table, family=gaussian(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Remove young site
ph.g <- brm(formula = ph ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table[!(chem.table$time<2),], family=gaussian(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## NDVI
ndvi.logn <- brm(formula = ndvi ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## R2
r2_bayes(ndvi.logn)
r2_bayes(ph.g)
r2_bayes(p.logn)
r2_bayes(n.hg)
r2_bayes(c.logn)

##Plotting
## Spatial autocorrelation
par(mfrow=c(2,3))
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(p.logn))
plot(cres, ylim = c(-1,1), main = "Phosphorus") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(n.hg))
plot(cres, ylim = c(-1,1), main = "Nitrogen") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(c.logn))
plot(cres, ylim = c(-1,1), main = "Carbon") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(ph.g))
plot(cres, ylim = c(-1,1), main = "pH") # No sign spatial autocorrelation

cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(ndvi.logn))
plot(cres, ylim = c(-1,1), main = "Productivity") # No sign spatial autocorrelation
