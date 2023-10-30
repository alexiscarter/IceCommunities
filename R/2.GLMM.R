## GLMM
## In this script there is the code for the generalized linear mixed-effect models for the biodiversity and soil features

## Load packages and data ####
library(phyloseq)
library(tidyverse)
library(brms); options(mc.cores = parallel::detectCores())
library(performance)
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
r2_bayes(bact.q1.logn)
r2_bayes(coll.q1.logn)
r2_bayes(euka.uni.q1.logn)
r2_bayes(euka.animals.q1.logn)
r2_bayes(fung.q1.logn)
r2_bayes(inse.q1.logn)
r2_bayes(olig.q1.logn)
r2_bayes(sper.q1.logn)

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

## Ratio
## CN
chem.table.ratio.cn <- div.table %>%
  drop_na() %>%
  filter(n>0.01) %>%
  mutate(cn=c/n) %>%
  filter(!Glacier %in% c("agola","rotmo","amola","sorap")) ## remove sites in the Dolomites

cn.logn <- brm(formula = cn ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  data = chem.table.ratio.cn, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## CP
chem.table.ratio.cp <- div.table %>%
  drop_na() %>%
  mutate(p = p/1000000) %>%
  filter(p>0.0005) %>%
  mutate(cp=c/p) %>%  
  filter(!Glacier %in% c("agola","rotmo","amola","sorap"))

cp.logn <- brm(formula = cp ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  data = chem.table.ratio.cp, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## R2
r2_bayes(ndvi.logn)
r2_bayes(ph.g)
r2_bayes(p.logn)
r2_bayes(n.hg)
r2_bayes(c.logn)
r2_bayes(cn.logn)
r2_bayes(cp.logn)

## Plotting
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

## Priors ####
## Priors for biotic factors
prior.bact <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase figure 4 of Pothula and Adams (2022)
prior.coll <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # increase figure 7 of Pothula and Adams (2022)
prior.euka.uni <- c(set_prior("normal(0,2)", class = "b", coef = "time.log.sc"))# no clear prior knowledge
prior.euka.ani <- c(set_prior("normal(0,2)", class = "b", coef = "time.log.sc"))  # depends on the taxa (Pothula and Adams 2022)
prior.fung <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase figure 5 of Pothula and Adams (2022)
prior.inse <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase figure 8 of Pothula and Adams (2022)
prior.olig <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase figure 9 of Pothula and Adams (2022)
prior.sper <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase figure 3 of Pothula and Adams (2022)
## Most taxa are expected to increase after Pothula and Adams so informative prior with a normal distribution with a mean of 1 and sd of 2. Except for euka.uni and euka.ani so informative prior with a normal distribution with a mean of 0 and sd of 2.

## Bact02
bact.q1.logn.prior <- brm(formula = bact.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.bact, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Coll01
coll.q1.logn.prior <- brm(formula = coll.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.coll, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.uni
euka.uni.q1.logn.prior <- brm(formula = euka.uni.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                              prior = prior.euka.uni, sample_prior = "yes",
                              data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.ani
euka.ani.q1.logn.prior <- brm(formula = euka.ani.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                              prior = prior.euka.ani, sample_prior = "yes",
                              data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Fung02
fung.q1.logn.prior <- brm(formula = fung.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.fung, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Inse01
inse.q1.logn.prior <- brm(formula = inse.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.inse, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Olig01
olig.q1.logn.prior <- brm(formula = olig.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.bact, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Sper01
sper.q1.logn.prior <- brm(formula = sper.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                          prior = prior.sper, sample_prior = "yes",
                          data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Priors for abiotic factors
prior.n <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # increase, figure 2b of Khedim et al. (2021)
prior.c <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # increase, figure 2 of Pothula and Adams (2022)
prior.p <- c(set_prior("normal(0,2)", class = "b", coef = "time.log.sc")) # no clear trend, figure 2 of Pothula and Adams (2022)
prior.ph <- c(set_prior("normal(-1,2)", class = "b", coef = "time.log.sc")) # decrease, figure 2 of Pothula and Adams (2022) 
prior.ndvi <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # increase, based on Ren and Gao (2022)
prior.cn <- c(set_prior("normal(1,2)", class = "b", coef = "time.log.sc")) # overall increase, figure 4 of Pothula and Adams (2022)
prior.cp <- c(set_prior("normal(0,2)", class = "b", coef = "time.log.sc")) # no clear trend, figure 2 of Pothula and Adams (2022)

## Nitrogen
n.hg.prior <- brm(formula = n ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  prior = prior.n, sample_prior = "yes",
                  data = chem.table, family=hurdle_gamma, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Carbon
c.logn.prior <- brm(formula = c ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                    prior = prior.c, sample_prior = "yes",
                    data = chem.table[!(chem.table$Glacier=="agola"),], # Agola foreland removed 
                    family=lognormal, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Phosphorus
p.logn.prior <- brm(formula = p ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                    prior = prior.p, sample_prior = "yes",
                    data = chem.table, family=lognormal, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## pH
ph.g.prior <- brm(formula = ph ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  prior = prior.ph, sample_prior = "yes",
                  data = chem.table, family=gaussian(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## NDVI
ndvi.logn.prior <- brm(formula = ndvi ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                       prior = prior.ndvi, sample_prior = "yes",
                       data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)
## Ratio
cn.g.prior <- brm(formula = cn ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  prior = prior.cn, sample_prior = "yes",
                  data = chem.table.ratio.cn, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

cp.g.prior <- brm(formula = cp ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  prior = prior.cp, sample_prior = "yes",
                  data = chem.table.ratio.cp, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

