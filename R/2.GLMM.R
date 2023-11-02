## GLMMs
## This script loads necessary packages, fits Generalized Linear Mixed-Effect Models (GLMMs), calculates R-squared values, and creates spatial correlograms to assess spatial autocorrelation for both biodiversity and soil features data.

## Load packages and data ####
# Load required packages
library(phyloseq)
library(tidyverse)
library(brms); options(mc.cores = parallel::detectCores())
library(performance)
library(ncf)
library(gdm)

# Set the number of CPU cores to use for parallel processing
options(mc.cores = parallel::detectCores())

# Load the dataset
div.table <- read.csv("data/div.clim.chem.csv", sep = ";") #from manage.clean.R
div.table$time.log.sc <- scale(log(div.table$time)) 

## Bact02
# Fit a GLMM for bacterial diversity
bact.q1.logn <- brm(formula = bact.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year), # It specifies the formula for the model. It indicates that you are modeling the variable bact.q1 as a function of time.log.sc and its squared term I(time.log.sc^2) while accounting for nested random factors (1|Glacier/Year). In other words, you are trying to predict bact.q1 based on the linear and quadratic effects of time.log.sc while considering the grouping structure of Glacier and Site (noted as the "Year" of the moraine).
                    data = div.table, # It specifies the data frame from which the variables in the formula are taken. 
                    family=lognormal(), # This sets the family of the response variable to log-normal. It means that you are assuming the distribution of bact.q1 is log-normal, which is often used for positively skewed continuous data.
                    warmup = 500, # This determines the number of warm-up iterations, which are used for tuning the sampler, a warm-up of 5000 iterations is preferred here, 500 for preliminary runs
                    iter = 5000, # This sets the total number of iterations for the Markov Chain Monte Carlo (MCMC) sampling, 50000 iterations here, 5000 for preliminary runs
                    chains = 4, # This specifies the number of parallel chains to run. In this case, you are running 4 chains in parallel, which is a common practice to assess convergence and improve efficiency.
                    thin = 10, # Thinning is a technique used to reduce the autocorrelation between samples in MCMC. It keeps every 10th iteration's sample.
                    refresh = 0, # The refresh interval determines how often progress updates are printed to the console. Setting it to 0 means that no updates will be printed during sampling.
                    silent = TRUE) # This suppresses any output messages or progress updates during the sampling process, making the process run silently.

## Coll01
# Fit a GLMM for springtail diversity
coll.q1.logn <- brm(formula = coll.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.uni
# Fit a GLMM for protist diversity
euka.uni.q1.logn <- brm(formula = euka.uni.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                     data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.ani
# Fit a GLMM for the diversity of other animals
euka.ani.q1.logn <- brm(formula = euka.ani.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                         data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Fung02
# Fit a GLMM for fungal diversity
fung.q1.logn <- brm(formula = fung.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Inse01
# Fit a GLMM for insect diversity
inse.q1.logn <- brm(formula = inse.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Olig01
# Fit a GLMM for earthworm diversity
olig.q1.logn <- brm(formula = olig.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Sper01
# Fit a GLMM for vascular plant diversity
sper.q1.logn <- brm(formula = sper.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

# Calculate R-squared for each of the biodiversity models
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
# Create spatial correlograms for each biodiversity model
par(mfrow=c(2,4))

# Bacteria
cres = spline.correlog(x = div.table$lon, y = div.table$lat, # This specifies the x and y coordinates of the plots 
                       resamp = 100, # This parameter specifies the number of resampling points or distances at which the correlogram will be calculated.
                       z = resid(bact.q1.logn), # This specifies the variable for which the residuals are calculated. It computes the residuals of the GLMM model bact.q1.logn. Residuals are the differences between the observed values and the predicted values from the model.
                       latlon = T) # This specifies that coordinates used are latitude and longitude
plot(cres, ylim = c(-1,1), main = "Bacteria", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Fungi
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(fung.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Fungi", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Protists
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(euka.uni.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Protists", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Vascular plants
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(sper.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Vascular plants", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Insects
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(inse.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Insects", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Springtails
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(coll.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Springtails", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Earthworms
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(olig.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Earthworms", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Other animals
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(euka.ani.q1.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Other animals", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

## GLMMs for soil data ####
# Create a new dataset by removing rows with missing values and scaling the phosphorus
chem.table <- div.table %>% 
  drop_na() %>% 
  mutate(p = p/1000000)

## Nitrogen
n.hg <- brm(formula = n ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table, family=hurdle_gamma, warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Carbon
c.logn <- brm(formula = c ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
              data = chem.table[!(chem.table$Glacier=="agola"),], family=lognormal, warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Phosphorus
p.logn <- brm(formula = p ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
              data = chem.table, family=lognormal, warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## pH
ph.g <- brm(formula = ph ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table, family=gaussian(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Remove young site
ph.g <- brm(formula = ph ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
            data = chem.table[!(chem.table$time<2),], family=gaussian(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## NDVI
ndvi.logn <- brm(formula = ndvi ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Ratio
## CN
chem.table.ratio.cn <- div.table %>%
  drop_na() %>%
  filter(n>0.01) %>%
  mutate(cn=c/n) %>%
  filter(!Glacier %in% c("agola","rotmo","amola","sorap")) ## remove sites in the Dolomites

cn.logn <- brm(formula = cn ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  data = chem.table.ratio.cn, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## CP
chem.table.ratio.cp <- div.table %>%
  drop_na() %>%
  mutate(p = p/1000000) %>%
  filter(p>0.0005) %>%
  mutate(cp=c/p) %>%  
  filter(!Glacier %in% c("agola","rotmo","amola","sorap"))

cp.logn <- brm(formula = cp ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                  data = chem.table.ratio.cp, family=lognormal(), warmup = 500, iter = 5000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## R2
# Calculate R-squared for each of the soil models
r2_bayes(ndvi.logn)
r2_bayes(ph.g)
r2_bayes(p.logn)
r2_bayes(n.hg)
r2_bayes(c.logn)
r2_bayes(cn.logn)
r2_bayes(cp.logn)

## Plotting
## Spatial autocorrelation
# Create spatial correlograms for each soil model
par(mfrow=c(2,4))

# Phosphorus
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(p.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Phosphorus", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Nitrogen
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(n.hg), latlon = T)
plot(cres, ylim = c(-1,1), main = "Nitrogen", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Carbon
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(c.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Carbon", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# pH
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(ph.g), latlon = T)
plot(cres, ylim = c(-1,1), main = "pH", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# Productivity
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(ndvi.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "Productivity", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# CN
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 1000, z = resid(cn.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "C:N ratio", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

# CP
cres = spline.correlog(x = div.table$lon, y = div.table$lat, resamp = 100, z = resid(cp.logn), latlon = T)
plot(cres, ylim = c(-1,1), main = "C:P ratio", xlab = "Distance between plots (km)") # No significant spatial autocorrelation is observed using spatial correlograms.

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
