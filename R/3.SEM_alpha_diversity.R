## SEM alpha diversity
## In this script there is the code for the structural equation modeling for the alpha diversity

library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(MuMIn)
library(MASS)
library(lavaan.survey)
library(ggplot2)
library(survey)
library(gridExtra)
options(max.print = 10000)

## PSEMS were tested with R 4.3.2
## using packages:
## lavaan 0.6-16
## lavaan.survey 1.1.3.1 (note that it cannot be installed directly from CRAN; zip file must be downloaded from  https://cran.r-project.org/src/contrib/Archive/lavaan.survey/)
## piecewiseSEM 2.3.0

## ALSO TESTED WITH R 3.6.2:
## piecewiseSEM 2.1.0
## lavaan 0.6-15
## lavaan.survey 1.1.3.1

setwd("PUT YOUR PATH HERE")

## Data
all_data=read.delim("data/div.clim.chem.csv", sep=",")
full_nona=all_data[!is.na(all_data$n),]      ### for this analysis, only use plots with complete chemical parameters
dim(full_nona)
# 793  33    793 plots; 33 columns

## transform and scale_variables
full_nona$ph.s=as.vector(scale(full_nona$ph))
full_nona$ndvi.s=as.vector(scale(log(full_nona$ndvi)))  ## log transformed to improve normality
full_nona$time.s=as.vector(scale(log(full_nona$time)))  ## log transformed to improve normality
full_nona$twi.s=as.vector(scale(log(full_nona$twi)))  ## log transformed to improve normality
full_nona$n.s=as.vector(scale(log(full_nona$n+0.01)))        ## log transformed to improve normality
full_nona$t.s=as.vector(scale(full_nona$mean.temp))
full_nona$p.s=as.vector(scale(log(full_nona$p)))  ## log transformed to improve normality
full_nona$c.s=as.vector(scale(log(full_nona$c))) ## log transformed to improve normality
full_nona$ph2=full_nona$ph.s^2
full_nona$time2=full_nona$time.s^2

### free variables to allow external correlations with PSEM
n.s=full_nona$n.s
ph.s=full_nona$ph.s
ph2=full_nona$ph2
ndvi.s=full_nona$ndvi.s
p.s=full_nona$p.s
twi.s=full_nona$twi.s

### example of patterns of variation of environmental variables: variation of NDVI
lmer_ndvi=lmer(ndvi.s~time.s+time2+t.s+n.s+twi.s+sper.q1+p.s+(1|Glacier)+(1|site), data=full_nona)
summary(lmer_ndvi)
# Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)  -0.19555    0.11903  35.49819  -1.643   0.1093    
# time.s        0.50703    0.04655 170.65049  10.892  < 2e-16 ***
#   time2         0.10588    0.02509 139.83440   4.221 4.36e-05 ***
#   t.s           0.15766    0.05248 470.81635   3.004   0.0028 ** 
#   n.s           0.04592    0.02003 702.70347   2.293   0.0221 *  
#   twi.s         0.02907    0.01647 714.61492   1.766   0.0779 .  
# sper.q1       0.02754    0.01129 669.72311   2.440   0.0150 *  
#   p.s           0.02300    0.02068 729.47355   1.112   0.2665     
r.squaredGLMM(lmer_ndvi)
# R2m      R2c
# [1,] 0.3218196 0.9216721


### LOG TRANSFORM BIODIVERSITY DATA TO IMPROVE NORMALITY
full_nona$fung.q1.l=log(full_nona$fung.q1)
full_nona$bact.q1.l=log(full_nona$bact.q1)
full_nona$sper.q1.l=log(full_nona$sper.q1)
full_nona$euka.ani.q1.l=log(full_nona$euka.ani.q1)
full_nona$euka.uni.q1.l=log(full_nona$euka.uni.q1)
full_nona$coll.q1.l=log(full_nona$coll.q1)
full_nona$olig.q1.l=log(full_nona$olig.q1)
full_nona$inse.q1.l=log(full_nona$inse.q1)

fung.q1.l=full_nona$fung.q1.l
bact.q1.l=full_nona$bact.q1.l
sper.q1.l=full_nona$sper.q1.l
euka.ani.q1.l=full_nona$euka.ani.q1.l
euka.uni.q1.l=full_nona$euka.uni.q1.l
coll.q1.l=full_nona$coll.q1.l
olig.q1.l=full_nona$olig.q1.l
inse.q1.l=full_nona$inse.q1.l
site=full_nona$site

###PIECEWISE STRUCTURAL EQUATION MODELS

### to take into account non-indipendence of plots within site and within landscape, models representing causal relationships
### included two random factors: Glacier (i.e., name of the proglacial landscape) and site (nested within glacier)

### MODEL 1: NUTRIENT-LED. THE ENVIRONMENTAL FEATURES DIRECTLY AFFECT BIODIVERSITY
psem_all_direct = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_olig=      lmer(olig.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_n= lmer(n.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_p= lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_ph=lmer(ph.s~time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_ndvi=lmer(ndvi.s~time.s+t.s+n.s+twi.s+sper.q1.l+p.s+(1|site)+(1|Glacier), data=full_nona),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  n.s%~~%ph2,
  n.s%~~%p.s,
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  lmer_sper=lmer(sper.q1.l~ph.s+time.s+twi.s+n.s+t.s+p.s+(1|site)+(1|Glacier), data=full_nona),
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
  )

summary(psem_all_direct, standardize = "scale")
# Response method Marginal Conditional
# euka.uni.q1.l   none     0.08        0.29
# euka.ani.q1.l   none     0.21        0.42
# coll.q1.l   none     0.34        0.54
# inse.q1.l   none     0.42        0.58
# olig.q1.l   none     0.22        0.51
# bact.q1.l   none     0.10        0.40
# fung.q1.l   none     0.29        0.59
# n.s   none     0.26        0.74
# p.s   none     0.12        0.72
# ph.s   none     0.05        0.82
# ndvi.s   none     0.31        0.92
# sper.q1.l   none     0.10        0.47

### MODEL 2: BIODIVERSITY-LED  (BIODIVERSITY determines ecosystem attributes)

psem_multifunct = psem(
  lmer_euka.uni=lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona),
  lmer_coll=lmer(coll.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona),
  lmer_inse=lmer(inse.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona),
  lmer_olig=lmer(olig.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona),
  lmer_bact=lmer(bact.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_fung=lmer(fung.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_sper=lmer(sper.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_ph=lmer(ph.s~time.s+twi.s+(1|site)+(1|Glacier), data=full_nona),
  
  lmer_n=lmer(n.s~time.s*t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s+(1|site)+(1|Glacier), data=full_nona),
  lmer_ndvi=lmer(ndvi.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+(1|site)+(1|Glacier), data=full_nona),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  n.s%~~%ph2,
  n.s%~~%p.s,
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

summary(psem_multifunct, standardize = "scale")
# Response method Marginal Conditional
# euka.uni.q1.l   none     0.07        0.30
# euka.ani.q1.l   none     0.19        0.41
# coll.q1.l   none     0.22        0.53
# inse.q1.l   none     0.30        0.58
# olig.q1.l   none     0.12        0.52
# bact.q1.l   none     0.09        0.40
# fung.q1.l   none     0.28        0.59
# sper.q1.l   none     0.10        0.45
# ph.s   none     0.03        0.83
# n.s   none     0.34        0.75
# p.s   none     0.15        0.74
# ndvi.s   none     0.33        0.92

### THIRD APPROACH: CO-VARIATION. ECOSYSTEM ATTRIBUTES COVARY  WITH BIODIVERSITY
### NOTE: I  TRIED TO ADD TIME^2, particularly given its effect on pH, but it increases the BIC of the model so then excluded.
### non-significant quadratic relationships (eg with pH) not included in the model

psem_covary = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),

  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

summary(psem_covary)
# euka.uni.q1.l   none     0.07        0.29
# euka.ani.q1.l   none     0.19        0.40
# coll.q1.l   none     0.22        0.52
# inse.q1.l   none     0.30        0.57
# olig.q1.l   none     0.12        0.51
# bact.q1.l   none     0.09        0.39
# fung.q1.l   none     0.29        0.59
# sper.q1.l   none     0.09        0.45
# ph.s   none     0.05        0.82
# n.s   none     0.26        0.73
# p.s   none     0.13        0.71
# ndvi.s   none     0.30        0.92


## MODEL FIT  
fisherC(psem_all_direct)
# Fisher.C df P.Value
# 1   14.904 10   0.136

fisherC(psem_multifunct)
# Fisher.C df P.Value
# 1     6.55  4   0.162

fisherC(psem_covary)
# Fisher.C df P.Value
# 9.507  8   0.301
dSep(psem_covary)
# Independ.Claim Test.Type       DF Crit.Value   P.Value 
# 1       n.s ~ ph2 + ...      coef 791.4616  1.0472507 0.2953037 
# 2       p.s ~ ph2 + ...      coef 788.2995  0.4053318 0.6853436 
# 3 sper.q1.l ~ ph2 + ...      coef 600.5215 -1.6539979 0.0986505 
# 4 olig.q1.l ~ ph2 + ...      coef 667.0053 -0.7866772 0.4317503 


### TO CALCULATE IMPORTANCE OF DIRECT AND INDIRECT EFFECTS, STORE THE COEFFICIENT IN A TABLE (to be used later)
tab_res_all=summary(psem_covary)$coefficients

#######################################################
### NOW USE LAVAAN.SURVEY TO CALCULATE BIC VALUES. models in Lavaan survey have the same structure of the ones in PSEM
#######################################################

#################### NUTRIENT-LED MODEL
lavaan_direct_full<-'
  euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s
  euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s
  coll.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s
  inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s
  olig.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s
  bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s
  fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s
  sper.q1.l    ~ph.s    +time.s+twi.s+t.s+                 n.s+p.s
  n.s~time.s+t.s+twi.s # lavaan does not support interactions
  p.s~time.s+t.s+twi.s
  ph.s~time.s+twi.s+t.s
  ph2~time.s+twi.s+t.s
  ndvi.s~time.s+t.s+n.s+twi.s+sper.q1.l+p.s
  n.s~~ph.s
  p.s~~ph.s
  n.s~~ph2
  n.s~~p.s
  ph.s~~ph2
  ndvi.s~~ph.s
  euka.ani.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.uni.q1.l
  inse.q1.l ~~ euka.uni.q1.l
  olig.q1.l ~~ euka.uni.q1.l
  bact.q1.l ~~ euka.uni.q1.l
  fung.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ euka.ani.q1.l
  olig.q1.l ~~ euka.ani.q1.l
  bact.q1.l ~~ euka.ani.q1.l
  fung.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ coll.q1.l
  olig.q1.l ~~ coll.q1.l
  bact.q1.l ~~ coll.q1.l
  fung.q1.l ~~ coll.q1.l
  olig.q1.l ~~ inse.q1.l
  bact.q1.l ~~ inse.q1.l
  fung.q1.l ~~ inse.q1.l
  bact.q1.l ~~ olig.q1.l
  fung.q1.l ~~ olig.q1.l
  fung.q1.l ~~ bact.q1.l
  
  bact.q1.l ~~ sper.q1.l
  fung.q1.l ~~ sper.q1.l'

fit_lavaan_direct <- sem(lavaan_direct_full, data=full_nona)
survey.design <- svydesign(ids=~site, prob=~1, data=full_nona) # Adding also Glacier ID to the structure of the model results in models with identical BIC values but with non-identification issues. All conclusions remain identical if Glacier ID is added
survey_direct <- lavaan.survey(lavaan.fit=fit_lavaan_direct, survey.design=survey.design)
fitMeasures(survey_direct, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "SRMR", "BIC"))
# chisq.scaled     df.scaled pvalue.scaled  rmsea.scaled          srmr           bic 
# 29.388         5.000         0.000         0.078         0.028     21160.226 
### BIODIVERSITY-LED MODEL

lavaan_multifunc_full <-'
  euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  coll.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  inse.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  olig.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  bact.q1.l~ph.s+ph2+time.s+twi.s+t.s
  fung.q1.l~ph.s+ph2+time.s+twi.s+t.s
  sper.q1.l~ph.s+ph2+time.s+twi.s+t.s
  ph.s~time.s+twi.s
  ph2~time.s+twi.s
  
  n.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s
  p.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s
  ndvi.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l
  n.s~~ph.s
  p.s~~ph.s
  n.s~~ph2
  n.s~~p.s
  ndvi.s~~ph.s
  ndvi.s~~ph2
  euka.ani.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.uni.q1.l
  inse.q1.l ~~ euka.uni.q1.l
  olig.q1.l ~~ euka.uni.q1.l
  bact.q1.l ~~ euka.uni.q1.l
  fung.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ euka.ani.q1.l
  olig.q1.l ~~ euka.ani.q1.l
  bact.q1.l ~~ euka.ani.q1.l
  fung.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ coll.q1.l
  olig.q1.l ~~ coll.q1.l
  bact.q1.l ~~ coll.q1.l
  fung.q1.l ~~ coll.q1.l
  olig.q1.l ~~ inse.q1.l
  bact.q1.l ~~ inse.q1.l
  fung.q1.l ~~ inse.q1.l
  bact.q1.l ~~ olig.q1.l
  fung.q1.l ~~ olig.q1.l
  fung.q1.l ~~ bact.q1.l
  ## and plants:
  
  bact.q1.l ~~ sper.q1.l
  fung.q1.l ~~ sper.q1.l'

fit_lavaan_multifunc <- sem(lavaan_multifunc_full, data=full_nona)
survey_multifunc <- lavaan.survey(lavaan.fit=fit_lavaan_multifunc, survey.design=survey.design)
fitMeasures(survey_multifunc, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "SRMR", "BIC"))
# chisq.scaled     df.scaled pvalue.scaled  rmsea.scaled          srmr           bic 
# 33.313         4.000         0.000         0.096         0.041     21273.013 

### COVARIATION MODEL

lavaan_covary_full<-'
  euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l
  olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l
  bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s
  fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s
  sper.q1.l    ~ph.s    +time.s+twi.s+t.s
  
  n.s~time.s+t.s+twi.s   # lavaan does not support interactions
  p.s~time.s+t.s+twi.s
  ph.s~time.s+t.s+twi.s 
  ph2~time.s+t.s+twi.s 
  ndvi.s~sper.q1.l+time.s+t.s+twi.s
  n.s~~ph.s
  p.s~~ph.s
  p.s~~euka.uni.q1.l
  p.s~~euka.ani.q1.l
  p.s~~coll.q1.l
  p.s~~inse.q1.l
  p.s~~olig.q1.l
  p.s~~bact.q1.l
  p.s~~fung.q1.l
  p.s~~sper.q1.l
  
  n.s~~p.s
  n.s~~euka.uni.q1.l
  n.s~~euka.ani.q1.l
  n.s~~coll.q1.l
  n.s~~inse.q1.l
  n.s~~olig.q1.l
  n.s~~bact.q1.l
  n.s~~fung.q1.l
  n.s~~sper.q1.l
  
  ndvi.s~~ph.s
  ndvi.s~~ph2
  ndvi.s~~n.s
  ndvi.s~~p.s
  ndvi.s~~euka.uni.q1.l
  ndvi.s~~euka.ani.q1.l
  ndvi.s~~coll.q1.l
  ndvi.s~~inse.q1.l
  ndvi.s~~olig.q1.l
  ndvi.s~~bact.q1.l
  ndvi.s~~fung.q1.l
  
  euka.ani.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.uni.q1.l
  inse.q1.l ~~ euka.uni.q1.l
  olig.q1.l ~~ euka.uni.q1.l
  bact.q1.l ~~ euka.uni.q1.l
  fung.q1.l ~~ euka.uni.q1.l
  coll.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ euka.ani.q1.l
  olig.q1.l ~~ euka.ani.q1.l
  bact.q1.l ~~ euka.ani.q1.l
  fung.q1.l ~~ euka.ani.q1.l
  inse.q1.l ~~ coll.q1.l
  olig.q1.l ~~ coll.q1.l
  bact.q1.l ~~ coll.q1.l
  fung.q1.l ~~ coll.q1.l
  olig.q1.l ~~ inse.q1.l
  bact.q1.l ~~ inse.q1.l
  fung.q1.l ~~ inse.q1.l
  bact.q1.l ~~ olig.q1.l
  fung.q1.l ~~ olig.q1.l
  fung.q1.l ~~ bact.q1.l
  ## and plants:
  
  bact.q1.l ~~ sper.q1.l
  fung.q1.l ~~ sper.q1.l'

fit_lavaan_covary <- sem(lavaan_covary_full, data=full_nona)
survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)
#the fit of this model:
fitMeasures(survey_covary, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "SRMR", "BIC"))
# chisq.scaled     df.scaled pvalue.scaled  rmsea.scaled          srmr           bic 
# 8.622         5.000         0.125         0.030         0.018     21125.878  

###########################################################################################
### ASSESS THE RELATIVE IMPORTANCE OF DIFFERENT PATHS IN THE BEST MODEL
########################################################################################### 

### rationale: starting from the best BIC model (fitted in lavaan.survey), we build nodels iteratively removing all the paths
### BIC is used to assess the fit of each model
### it is thus possible testing the average importance of paths representing the effects of habitat, time and biotic interactions

formlist=list(
  "euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l",
  "euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l",
  "coll.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l",
  "inse.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l",
  "olig.q1.l~ph.s+time.s+twi.s+t.s+sper.q1.l",
  "bact.q1.l~ph.s+ph2+time.s+twi.s+t.s",
  "fung.q1.l~ph.s+ph2+time.s+twi.s+t.s",
  "sper.q1.l~ph.s+time.s+twi.s+t.s",
  "p.s~~euka.uni.q1.l",
  "p.s~~euka.ani.q1.l",
  "p.s~~coll.q1.l",
  "p.s~~inse.q1.l",
  "p.s~~olig.q1.l",
  "p.s~~bact.q1.l",
  "p.s~~fung.q1.l",
  "p.s~~sper.q1.l",
  "n.s~~euka.uni.q1.l",
  "n.s~~euka.ani.q1.l",
  "n.s~~coll.q1.l",
  "n.s~~inse.q1.l",
  "n.s~~olig.q1.l",
  "n.s~~bact.q1.l",
  "n.s~~fung.q1.l",
  "n.s~~sper.q1.l",
  "ndvi.s~~euka.uni.q1.l",
  "ndvi.s~~euka.ani.q1.l",
  "ndvi.s~~coll.q1.l",
  "ndvi.s~~inse.q1.l",
  "ndvi.s~~olig.q1.l",
  "ndvi.s~~bact.q1.l",
  "ndvi.s~~fung.q1.l",
  "euka.ani.q1.l~~euka.uni.q1.l",
  "coll.q1.l~~euka.uni.q1.l",
  "inse.q1.l~~euka.uni.q1.l",
  "olig.q1.l~~euka.uni.q1.l",
  "bact.q1.l~~euka.uni.q1.l",
  "fung.q1.l~~euka.uni.q1.l",
  "coll.q1.l~~euka.ani.q1.l",
  "inse.q1.l~~euka.ani.q1.l",
  "olig.q1.l~~euka.ani.q1.l",
  "bact.q1.l~~euka.ani.q1.l",
  "fung.q1.l~~euka.ani.q1.l",
  "inse.q1.l~~coll.q1.l",
  "olig.q1.l~~coll.q1.l",
  "bact.q1.l~~coll.q1.l",
  "fung.q1.l~~coll.q1.l",
  "olig.q1.l~~inse.q1.l",
  "bact.q1.l~~inse.q1.l",
  "fung.q1.l~~inse.q1.l",
  "bact.q1.l~~olig.q1.l",
  "fung.q1.l~~olig.q1.l",
  "fung.q1.l~~bact.q1.l",
  "bact.q1.l~~sper.q1.l",
  "fung.q1.l~~sper.q1.l",  
  "n.s~time.s+t.s+twi.s",
  "p.s~time.s+t.s+twi.s",
  "ph.s~time.s+t.s+twi.s",
  "ndvi.s~sper.q1.l+time.s+t.s+twi.s",
  "n.s~~ph.s",
  "p.s~~ph.s",
  "n.s~~p.s",
  "ndvi.s~~ph.s",
  "ndvi.s~~ph2",
  "ndvi.s~~n.s",
  "ndvi.s~~p.s",
  "ph2~time.s+t.s+twi.s"
)

out=vector("list",0)
i=1

## remove iteratively all the paths representing potential effects on biodiversity:
for(i in 1:54){
  print(i)
  split=unlist(strsplit(formlist[[i]],split="\\~{1}"))		# this produces three elements when the formula is a co-variation
  
  if(length(split)==3){						# if this line represents a co-variation, completely remove i-th covariation (by using "0*")
    
    term1=split[1]
    term2=split[3]
    
    formula=formlist
    formula[[i]]=paste0(term1,"~~0*",term2)
    
    fit_lavaan_covary <- sem(paste0(formula,collapse="\n "), data=full_nona)
    survey.design <- svydesign(ids=~site, prob=~1, data=full_nona)
    survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)
    bic=BIC(survey_covary)
    aic=AIC(survey_covary)
    logl=fitMeasures(survey_covary, c("logl"))
    out[[length(out)+1]]=c(formlist[[i]],formlist[[i]],bic,aic)
    
  }else{								#  if this line represents a directional effect, iteratively removes each of the independent variables
    
    
    dep=split[1]
    ind=split[2]
    
    var=unlist(strsplit(ind,split="\\+|\\*"))
    
    for(j in 1:length(var)){
      
      jind=paste(var[-j],collapse="+")
      
      formula=formlist
      formula[[i]]=paste(dep,jind,sep="~")
      
      fit_lavaan_covary <- sem(paste0(formula,collapse="\n "), data=full_nona)
      survey.design <- svydesign(ids=~site, prob=~1, data=full_nona)
      survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)
      bic=BIC(survey_covary)
      aic=AIC(survey_covary)
      logl=fitMeasures(survey_covary, c("logl"))
      out[[length(out)+1]]=c(formlist[[i]],var[j],bic,aic)
          }
      }
  }       

res=matrix(unlist(out),length(out),4,byrow=T)
colnames(res)=c("model","removed_variable","bic","aic")
res=data.frame(res)
res$bic=as.numeric(as.character(res$bic))
res$delta_bic=res$bic-21125.878     ### calculate the change in BIC resulting from the removal of each variable. 27832.691is the BIC of the full mdoel

# the summary of delta-bic depending on type of drivers
## for each path, define if it represents biotic effects (b), habitat (h) or time (t)

type=c('Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Biot', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Biot', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Biot', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Biot', 'Hab', 'Time', 'Hab', 'Hab', 'Biot', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Hab', 'Time', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Hab', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot', 'Biot')
typeF=c('Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Biotic', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Biotic', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Biotic', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Biotic', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Biotic', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Habitat', 'Time', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Habitat', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic', 'Biotic')


######################################################
#  evaluate if the importance of processes is similar between boreal, temperate and tropical
######################################################
### CREATE 3 GEOGRAPHICALLY-RESTRICTED SUBSETS: BOREAL

boreal=subset(full_nona, full_nona$lat>60|full_nona$lat< -46)

### SUBSET: TROPICAL

tropical=subset(full_nona, full_nona$lat<33&full_nona$lat> -20)

### SUBSET: temperate

temperate=subset(full_nona, full_nona$lat<60&full_nona$lat> -45)
temperate=subset(temperate, temperate$lat>33|temperate$lat< -20)

######################################################
#  evaluate if the importance of processes is similar between boreal, temperate and tropical areas USING PSEM 
######################################################

######################################################
#  BOREAL
######################################################

### free variables to allow external correlations with PSEM
n.s=boreal$n.s
ph.s=boreal$ph.s
ph2=boreal$ph2
ndvi.s=boreal$ndvi.s
p.s=boreal$p.s
twi.s=boreal$twi.s

psem_covary_boreal = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=boreal, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=boreal, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

tab_res_psem_boreal=summary(psem_covary_boreal)$coefficients

######################################################
#  TEMPERATE
######################################################

### free variables to allow external correlations with PSEM
n.s=temperate$n.s
ph.s=temperate$ph.s
ph2=temperate$ph2
ndvi.s=temperate$ndvi.s
p.s=temperate$p.s
twi.s=temperate$twi.s

psem_covary_temperate = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=temperate, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=temperate, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

tab_res_psem_temperate=summary(psem_covary_temperate)$coefficients

######################################################
#  TROPICAL
######################################################


### free variables to allow external correlations with PSEM
n.s=tropical$n.s
ph.s=tropical$ph.s
ph2=tropical$ph2
ndvi.s=tropical$ndvi.s
p.s=tropical$p.s
twi.s=tropical$twi.s

psem_covary_tropical = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=tropical, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=tropical, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

tab_res_psem_tropical=summary(psem_covary_tropical)$coefficients


### ### ### 
### GO BACK TO THE ORIGINAL VARIABLES (INCLUDING ALL THE FORELANDS)
n.s=full_nona$n.s
ph.s=full_nona$ph.s
ph2=full_nona$ph2
ndvi.s=full_nona$ndvi.s
p.s=full_nona$p.s
twi.s=full_nona$twi.s
### ### ### ### ### ### 
## the typologie of effects
typology_psem=c("h", "h", "t", "h", "h", "b", "h", "h", "t", "h", "h", "b", "h", "h", "t", "h", "h", "b", "h", "h", "t", "h", "h", "b", "h", "t", "h", "h", "b", "h", "h", "t", "h", "h", "h", "h", "t", "h", "h", "h", "t", "h", "h", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "h", "h", "h", "h", "h", "h", "h", "h", "other", "h", "h", "h", "h", "h", "h", "h", "h", "other", "other", "other", "other", "other", "h", "h", "h", "h", "h", "h", "h", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b", "b")
### 

## store all the standardized estimates in a table to bild plots:
tab_res_psem=cbind(abs(tab_res_all$Std.Estimate), abs(tab_res_psem_boreal$Std.Estimate), abs(tab_res_psem_temperate$Std.Estimate), abs(tab_res_psem_tropical$Std.Estimate), typology_psem)
colnames(tab_res_psem)<-c("All", "Boreal", "Temperate", "Tropical", "Type")
tab_res_psem=tab_res_psem[!tab_res_psem[,5]=="other",]
tab_res_psem=data.frame(tab_res_psem)
tab_res_psem$All=as.numeric(as.character(tab_res_psem$All))
tab_res_psem$Boreal=as.numeric(as.character(tab_res_psem$Boreal))
tab_res_psem$Temperate=as.numeric(as.character(tab_res_psem$Temperate))
tab_res_psem$Tropical=as.numeric(as.character(tab_res_psem$Tropical))

### ### ### ### ### ### ### ### ### ### ### 
### PLOTS WITH THE STANDARDIZED EFFECTS
### ### ### ### ### ### ### ### ### ### ### 

## all the landscapes
gg_all=ggplot(tab_res_psem, aes(x = Type, y =  All, fill = type, col=type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.7)+
  ylab("")+
  xlab("")+
  theme_minimal()

## tropical only
gg_trop=ggplot(tab_res_psem, aes(x = Type, y =  Tropical , fill = Type, col=Type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.7)+
  ylab("")+
  xlab("")+
  theme_minimal()

## temperate only
gg_temp=ggplot(tab_res_psem, aes(x = Type, y =  Temperate, fill = type, col=type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.7)+ylab("")+
  xlab("")+
  theme_minimal()

## boreal only
gg_bor=ggplot(tab_res_psem, aes(x = Type, y =  Boreal, fill = type, col=type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.7)+ylab("")+
  xlab("")+
  theme_minimal()

x11(width=5,height=3)
grid.arrange(gg_all, gg_bor, gg_temp, gg_trop, nrow = 1)



### PLOT THE RELATIVE IMPORTANCE ASSESSED BY BIC DROP
x11(width=5,height=3)
ggplot(res, aes(x = typeF, y =  delta_bic, fill = typeF, col=typeF)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, geom = "errorbar", col="black", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  # geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 0), width=0.75, color = "black", show.legend = FALSE)+
  ylab("")+
  xlab("")+
  theme_minimal()

### EXPORT TABLES OF DATA USED TO BILD FIG 3A AND FIG 3B

res=cbind(res, typeF)
write.table(res, "bic_iterative_exclusion_alpha.txt", row.names = F)  # to export and explore the results
write.table(tab_res_psem, "psem_standardized_coefficients_alpha.txt", row.names = F)  # to export and explore the results

## check spatial autocorrelation of the residuals of the best-model (psem_covary)
## this can be done with the library ecogenetics (only works with R 3.6 or previous versions)

library(EcoGenetics)

# protists
res=residuals(psem_covary$lmer_euka.uni)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8], alternative="greater",seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))

# d.mean     obs p.val  size
# d=0-100           45.875 -0.0906     1  1821
# d=100-500        279.349 -0.0201     1  4150
# d=500-1000       702.937  0.0043     1  1978
# d=1000-10000    3333.247 -0.0185     1  3480
# d=10000-50000  34942.905 -0.0016     1 14720
# d=50000-1e+05  73479.503 -0.0014     1 12619
# d=1e+05-2e+05 149708.716  0.0002     1 23754
# d=2e+05-4e+05 243392.167 -0.0003     1 16996
# d=4e+05-1e+06 647560.848  0.0019     1  7540

res=residuals(psem_covary$lmer_euka.multi)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0776     1  1821
# d=100-500        279.349 -0.0161     1  4150
# d=500-1000       702.937 -0.0046     1  1978
# d=1000-10000    3333.072 -0.0285     1  3480
# d=10000-50000  34942.905 -0.0017     1 14720
# d=50000-1e+05  73479.503 -0.0003     1 12619
# d=1e+05-2e+05 149708.716 -0.0003     1 23754
# d=2e+05-4e+05 243392.167 -0.0007     1 16996
# d=4e+05-1e+06 647560.848  0.0000     1  7540


res=residuals(psem_covary$lmer_sper)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0694     1  1821
# d=100-500        279.349 -0.0263     1  4150
# d=500-1000       702.937 -0.0594     1  1978
# d=1000-10000    3333.247 -0.0029     1  3480
# d=10000-50000  34942.905  0.0003     1 14720
# d=50000-1e+05  73479.503  0.0002     1 12619
# d=1e+05-2e+05 149708.716  0.0007     1 23754
# d=2e+05-4e+05 243392.167  0.0013     1 16996
# d=4e+05-1e+06 647560.848  0.0035     1  7540



res=residuals(psem_covary$lmer_coll)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0843     1  1821
# d=100-500        279.349 -0.0116     1  4150
# d=500-1000       702.937 -0.0028     1  1978
# d=1000-10000    3333.247 -0.0404     1  3480
# d=10000-50000  34942.905 -0.0015     1 14720
# d=50000-1e+05  73479.503  0.0008     1 12619
# d=1e+05-2e+05 149708.716  0.0005     1 23754
# d=2e+05-4e+05 243392.167  0.0001     1 16996
# d=4e+05-1e+06 647560.848 -0.0014     1  7540

res=residuals(psem_covary$lmer_inse)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0634     1  1821
# d=100-500        279.349 -0.0230     1  4150
# d=500-1000       702.937 -0.0029     1  1978
# d=1000-10000    3333.247 -0.0397     1  3480
# d=10000-50000  34942.905  0.0000     1 14720
# d=50000-1e+05  73479.503  0.0006     1 12619
# d=1e+05-2e+05 149708.716  0.0020     1 23754
# d=2e+05-4e+05 243392.167  0.0000     1 16996
# d=4e+05-1e+06 647560.848 -0.0021     1  7540

res=residuals(psem_covary$lmer_olig)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],  alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0873     1  1821
# d=100-500        279.349 -0.0108     1  4150
# d=500-1000       702.937 -0.0382     1  1978
# d=1000-10000    3333.247 -0.0248     1  3480
# d=10000-50000  34942.905 -0.0003     1 14720
# d=50000-1e+05  73479.503  0.0011     1 12619
# d=1e+05-2e+05 149708.716  0.0006     1 23754
# d=2e+05-4e+05 243392.167  0.0004     1 16996
# d=4e+05-1e+06 647560.848 -0.0018     1  7540

res=residuals(psem_covary$lmer_bact)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8], alternative="greater",seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0888     1  1821
# d=100-500        279.349 -0.0118     1  4150
# d=500-1000       702.937  0.0033     1  1978
# d=1000-10000    3333.247 -0.0389     1  3480
# d=10000-50000  34942.905  0.0023     1 14720
# d=50000-1e+05  73479.503  0.0013     1 12619
# d=1e+05-2e+05 149708.716  0.0030     1 23754
# d=2e+05-4e+05 243392.167  0.0020     1 16996
# d=4e+05-1e+06 647560.848  0.0016     1  7540

res=residuals(psem_covary$lmer_fung)
hist(res)
(moran=eco.correlog(Z=res, XY=full_nona[,7:8],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
# d.mean     obs p.val  size
# d=0-100           45.875 -0.0975     1  1821
# d=100-500        279.349 -0.0215     1  4150
# d=500-1000       702.937 -0.0189     1  1978
# d=1000-10000    3333.247 -0.0141     1  3480
# d=10000-50000  34942.905 -0.0002     1 14720
# d=50000-1e+05  73479.503  0.0000     1 12619
# d=1e+05-2e+05 149708.716  0.0009     1 23754
# d=2e+05-4e+05 243392.167  0.0011     1 16996
# d=4e+05-1e+06 647560.848  0.0019     1  7540




lmer_p_n=lmer(lg_p~lg_n++(1|site)+(1|Glacier), data=full_nona)
r.squaredGLMM(lmer_p_n)
# R2m       R2c
# [1,] 0.07630883 0.7568362

cor(full_nona[,c(30,31,45)])
#        lg_n      lg_p      lg_c
# lg_n 1.0000000 0.3320786 0.7496999
# lg_p 0.3320786 1.0000000 0.2073620
# lg_c 0.7496999 0.2073620 1.0000000

#############################################################################################
### alternative model: USE CARBON INSTEAD OF NITROGEN
#############################################################################################

### RELATIONSHIP CARBON / NITROGEN

cor(full_nona[,c(38,40,41)])    ##covariation between N, P and C
#        n.s       p.s       c.s
# n.s 1.0000000 0.3320786 0.7496999
# p.s 0.3320786 1.0000000 0.2073620
# c.s 0.7496999 0.2073620 1.0000000

lmer_c_n=lmer(c.s~n.s+(1|site)+(1|Glacier), data=full_nona)
summary(lmer_c_n)
r.squaredGLMM(lmer_c_n)
# R2m       R2c
# [1,] 0.3990582 0.8927772   ## the relationship is too strong to keep C and N toghether in the same model, so we developed separate models

## external variable for co-variances
c.s=full_nona$c.s

psem_carbon = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F), 
  lmer_c=lmer(c.s~time.s*t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  c.s%~~%ph.s,
  c.s%~~%ph2,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  c.s%~~%p.s,
  c.s%~~%euka.uni.q1.l,
  c.s%~~%euka.ani.q1.l,
  c.s%~~%coll.q1.l,
  c.s%~~%inse.q1.l,
  c.s%~~%olig.q1.l,
  c.s%~~%bact.q1.l,
  c.s%~~%fung.q1.l,
  c.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%c.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)

summary(psem_carbon)
# euka.uni.q1.l   none     0.07        0.29
# euka.ani.q1.l   none     0.19        0.40
# coll.q1.l   none     0.22        0.52
# inse.q1.l   none     0.30        0.57
# olig.q1.l   none     0.12        0.51
# bact.q1.l   none     0.09        0.39
# fung.q1.l   none     0.29        0.59
# sper.q1.l   none     0.09        0.45
# ph.s   none     0.05        0.82
# n.s   none     0.26        0.73
# p.s   none     0.13        0.71
# ndvi.s   none     0.30        0.92


## MODEL FIT  
fisherC(psem_carbon)
# Fisher.C df P.Value
# 1    7.068  6   0.315

#############################################################################################
### ADDITIONAL TEST: is there a role for solar radiation and / or aspect (northness?)
#############################################################################################

### LOAD DATA ON TERRAIN: SOLAR RADIATION & ASPECT

terrain=read.delim("data/terrain_swrad.csv", sep=";")

hist(terrain$slope)   ## slope values


solar=c()   ## calculate annual solar radiation as the sum of monthly radiation during snow-free months:

for (i in 1:793){
  print(i)
  s=terrain[i,2:13]
  solar[i]<-sum(s, na.rm = T)
}

# better on untransformed data!!!
full_nona$solar<-solar
full_nona$solar.s=scale(full_nona$so)
cor(full_nona$solar, full_nona$t.s)    #positive relationship between temperature and solar radiation
# 0.5863836

full_nona$aspect<-terrain$aspect
# following Amatulli et al 2018, use cos(aspect) and multiply it by -1 in the S hemisphere
full_nona$aspect_corr<-full_nona$aspect
full_nona$aspect_cos<-cos(full_nona$aspect*pi/180)
full_nona$aspect_cos[full_nona$lat<0]<-full_nona$aspect_cos[full_nona$lat<0]*-1

### THE PSEM INCLUDING SOLAR RADIATION AS ADDITIONAL VARIABLE:
psem_covary_solar = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+solar.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)


summary(psem_covary_solar)  #  below I only pasted the coeffcients of solar radiation
# Response       Predictor Estimate Std.Error       DF Crit.Value P.Value Std.Estimate    
  # euka.uni.q1.l         solar.s   0.0365    0.0696  42.2781     0.5242  0.6028       0.0421    
  # euka.ani.q1.l         solar.s  -0.0032     0.065  42.5198    -0.0488  0.9613      -0.0039    
  # coll.q1.l         solar.s  -0.0450    0.0452  49.5868    -0.9962  0.3240      -0.0840    
  # inse.q1.l         solar.s  -0.0734    0.0577  50.0294    -1.2723  0.2092      -0.1040    
  # olig.q1.l         solar.s  -0.0676    0.0375  50.4454    -1.8028  0.0774      -0.1724    
  # bact.q1.l         solar.s  -0.0977    0.0779  42.6586    -1.2556  0.2161      -0.1116    
  # fung.q1.l         solar.s   0.0563    0.0813  45.4058     0.6923  0.4923       0.0579    
  # sper.q1.l         solar.s  -0.0294    0.0544  49.2933    -0.5401  0.5916      -0.0535   ## NO EFFECT FOR ANY OF THE TAXA

### THE PSEM INCLUDING ASPECT ("northness") AS ADDITIONAL VARIABLE:
psem_covary_aspect = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+aspect_cos+(1|site)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|site)+(1|Glacier), data=full_nona, REML=F),
  # northness%~~%t.s,
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  p.s%~~%euka.uni.q1.l,
  p.s%~~%euka.ani.q1.l,
  p.s%~~%coll.q1.l,
  p.s%~~%inse.q1.l,
  p.s%~~%olig.q1.l,
  p.s%~~%bact.q1.l,
  p.s%~~%fung.q1.l,
  p.s%~~%sper.q1.l,
  
  n.s%~~%p.s,
  n.s%~~%euka.uni.q1.l,
  n.s%~~%euka.ani.q1.l,
  n.s%~~%coll.q1.l,
  n.s%~~%inse.q1.l,
  n.s%~~%olig.q1.l,
  n.s%~~%bact.q1.l,
  n.s%~~%fung.q1.l,
  n.s%~~%sper.q1.l,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  ndvi.s%~~%euka.uni.q1.l,
  ndvi.s%~~%euka.ani.q1.l,
  ndvi.s%~~%coll.q1.l,
  ndvi.s%~~%inse.q1.l,
  ndvi.s%~~%olig.q1.l,
  ndvi.s%~~%bact.q1.l,
  ndvi.s%~~%fung.q1.l,
  
  euka.ani.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.uni.q1.l,
  inse.q1.l %~~% euka.uni.q1.l,
  olig.q1.l %~~% euka.uni.q1.l,
  bact.q1.l %~~% euka.uni.q1.l,
  fung.q1.l %~~% euka.uni.q1.l,
  coll.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% euka.ani.q1.l,
  olig.q1.l %~~% euka.ani.q1.l,
  bact.q1.l %~~% euka.ani.q1.l,
  fung.q1.l %~~% euka.ani.q1.l,
  inse.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% coll.q1.l,
  bact.q1.l %~~% coll.q1.l,
  fung.q1.l %~~% coll.q1.l,
  olig.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% inse.q1.l,
  fung.q1.l %~~% inse.q1.l,
  bact.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% olig.q1.l,
  fung.q1.l %~~% bact.q1.l,
  ## and plants:
  
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
)


summary(psem_covary_aspect)  # [only the effect of aspect is shown]
### ALSO IN THIS CASE, NO EFFECT OF ASPECT:
# Response       Predictor Estimate Std.Error       DF Crit.Value P.Value Std.Estimate    
  # euka.uni.q1.l      aspect_cos   0.0851    0.0535 432.0064     1.5906  0.1124       0.0670    
  # euka.ani.q1.l      aspect_cos  -0.0157    0.0469 459.9296    -0.3356  0.7373      -0.0133    
  # coll.q1.l      aspect_cos  -0.0068      0.03 580.4181    -0.2269  0.8206      -0.0087    
  # inse.q1.l      aspect_cos  -0.0324     0.036 550.1336    -0.9010  0.3680      -0.0314    
  # olig.q1.l      aspect_cos  -0.0381    0.0216 574.2295    -1.7619  0.0786      -0.0664    
  # bact.q1.l      aspect_cos   0.0487    0.0553 534.8191     0.8804  0.3790       0.0380    
  # fung.q1.l      aspect_cos   0.0004    0.0516 618.1803     0.0075  0.9940       0.0003    
  # sper.q1.l      aspect_cos   0.0196    0.0313 518.9681     0.6265  0.5313       0.0244    
