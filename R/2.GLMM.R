## GLMMs
## This script loads necessary packages, fits Generalized Linear Mixed-Effect Models (GLMMs), calculates R-squared values, and creates spatial correlograms to assess spatial autocorrelation for both biodiversity and soil features data.

## Load packages and data ####
# Load required packages
library(phyloseq)
library(tidyverse)
library(brms)
library(ncf)
library(gdm)

# Set the number of CPU cores to use for parallel processing
options(mc.cores = parallel::detectCores())

# Load the dataset
div.table <- read.csv("data/div.clim.chem.csv") #from manage.clean.R

## Bact02
# Fit a GLMM for bacterial diversity
bact.q1.logn <- brm(formula = bact.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year), # It specifies the formula for the model. It indicates that you are modeling the variable bact.q1 as a function of time.log.sc and its squared term I(time.log.sc^2) while accounting for nested random factors (1|Glacier/Year). In other words, you are trying to predict bact.q1 based on the linear and quadratic effects of time.log.sc while considering the grouping structure of Glacier and Site (noted as the "Year" of the moraine).
                    data = div.table, # It specifies the data frame from which the variables in the formula are taken. 
                    family=lognormal(), # This sets the family of the response variable to log-normal. It means that you are assuming the distribution of bact.q1 is log-normal, which is often used for positively skewed continuous data.
                    warmup = 5000, # This determines the number of warm-up iterations, which are used for tuning the sampler, 5000 warm-up iterations used here.
                    iter = 50000, # This sets the total number of iterations for the Markov Chain Monte Carlo (MCMC) sampling, 50000 iterations here.
                    chains = 4, # This specifies the number of parallel chains to run. In this case, you are running 4 chains in parallel, which is a common practice to assess convergence and improve efficiency.
                    thin = 10, # Thinning is a technique used to reduce the autocorrelation between samples in MCMC. It keeps every 10th iteration's sample.
                    refresh = 0, # The refresh interval determines how often progress updates are printed to the console. Setting it to 0 means that no updates will be printed during sampling.
                    silent = TRUE) # This suppresses any output messages or progress updates during the sampling process, making the process run silently.

## Coll01
# Fit a GLMM for springtail diversity
coll.q1.logn <- brm(formula = coll.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.uni
# Fit a GLMM for protist diversity
euka.uni.q1.logn <- brm(formula = euka.uni.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                     data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## euka.ani
# Fit a GLMM for the diversity of other animals
euka.ani.q1.logn <- brm(formula = euka.ani.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                         data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Fung02
# Fit a GLMM for fungal diversity
fung.q1.logn <- brm(formula = fung.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Inse01
# Fit a GLMM for insect diversity
inse.q1.logn <- brm(formula = inse.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Olig01
# Fit a GLMM for earthworm diversity
olig.q1.logn <- brm(formula = olig.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## Sper01
# Fit a GLMM for vascular plant diversity
sper.q1.logn <- brm(formula = sper.q1 ~ time.log.sc + I(time.log.sc^2) + (1|Glacier/Year),
                 data = div.table, family=lognormal(), warmup = 5000, iter = 50000, chains = 4, thin = 10, refresh = 0, silent = TRUE)

## R2
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
cres = spline.correlog(x = div.table$lat, y = div.table$lon, # This specifies the x and y coordinates of the plots 
                       resamp = 100, # This parameter specifies the number of resampling points or distances at which the correlogram will be calculated.
                       z = resid(bact.q1.logn)) # This specifies the variable for which the residuals are calculated. It computes the residuals of the GLMM model bact.q1.logn. Residuals are the differences between the observed values and the predicted values from the model.
plot(cres, ylim = c(-1,1), main = "Bacteria") # No significant spatial autocorrelation is observed using spatial correlograms.

# Fungi
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(fung.q1.logn))
plot(cres, ylim = c(-1,1), main = "Fungi") # No significant spatial autocorrelation is observed using spatial correlograms.

# Protists
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(euka.uni.q1.logn))
plot(cres, ylim = c(-1,1), main = "Protists") # No significant spatial autocorrelation is observed using spatial correlograms.

# Vascular plants
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(sper.q1.logn))
plot(cres, ylim = c(-1,1), main = "Vascular plants") # No significant spatial autocorrelation is observed using spatial correlograms.

# Insects
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(inse.q1.logn))
plot(cres, ylim = c(-1,1), main = "Insects") # No significant spatial autocorrelation is observed using spatial correlograms.

# Springtails
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(coll.q1.logn))
plot(cres, ylim = c(-1,1), main = "Springtails") # No significant spatial autocorrelation is observed using spatial correlograms.

# Earthworms
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(olig.q1.logn))
plot(cres, ylim = c(-1,1), main = "Earthworms") # No significant spatial autocorrelation is observed using spatial correlograms.

# Other animals
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(euka.ani.q1.logn))
plot(cres, ylim = c(-1,1), main = "Other animals") # No significant spatial autocorrelation is observed using spatial correlograms.

## GLMMs for soil data ####
# Create a new dataset by removing rows with missing values and scaling the phosphorus
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
# Calculate R-squared for each of the soil models
r2_bayes(ndvi.logn)
r2_bayes(ph.g)
r2_bayes(p.logn)
r2_bayes(n.hg)
r2_bayes(c.logn)

##Plotting
## Spatial autocorrelation
# Create spatial correlograms for each soil model
par(mfrow=c(2,3))

# Phosphorus
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(p.logn))
plot(cres, ylim = c(-1,1), main = "Phosphorus") # No significant spatial autocorrelation is observed using spatial correlograms.

# Nitrogen
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(n.hg))
plot(cres, ylim = c(-1,1), main = "Nitrogen") # No significant spatial autocorrelation is observed using spatial correlograms.

# Carbon
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(c.logn))
plot(cres, ylim = c(-1,1), main = "Carbon") # No significant spatial autocorrelation is observed using spatial correlograms.

# pH
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(ph.g))
plot(cres, ylim = c(-1,1), main = "pH") # No significant spatial autocorrelation is observed using spatial correlograms.

# Productivity
cres = spline.correlog(x = div.table$lat, y = div.table$lon, resamp = 100, z = resid(ndvi.logn))
plot(cres, ylim = c(-1,1), main = "Productivity") # No significant spatial autocorrelation is observed using spatial correlograms.
