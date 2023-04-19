## SEM alpha diversity
## In this script there is the code for the structural equation modeling for the alpha diversity

### NOTE: ANALYSIS PERFORMED USING  THE VERSION 2.1.2 OF PIECEWISE SEM.
### DIFFERENT VERSIONS OF PIECEWISE SEM YELD SLIGHTLY DIFFERENT VALUES

library(iml)
library(lme4)
library(lmerTest)
library(vegan)
library(piecewiseSEM)
library(MuMIn)
library(nlme)
library(MuMIn)
library(MASS)
library(EcoGenetics)
options(max.print = 10000)

all_data=read.delim("data/div.clim.chem.csv", sep=";")

## remove rows with NA
full_nona=all_data[!is.na(all_data$lg_n),]

## transform and scale_variables
full_nona$ph.s=scale(full_nona$ph)
full_nona$ndvi.s=scale(log(full_nona$ndvi))  ## log transformed to improve normality
full_nona$time.s=scale(full_nona$lg_time)  ## log transformed to improve normality
full_nona$twi.s=scale(log(full_nona$twi))  ## log transformed to improve normality
full_nona$n.s=scale(full_nona$lg_n)        ## log transformed to improve normality
full_nona$t.s=scale(full_nona$meanT_new)
full_nona$p.s=scale(full_nona$lg_p)
full_nona$ph2=full_nona$ph.s^2
full_nona$time2=full_nona$time.s^2
full_nona$lg_c=log(full_nona$c)
full_nona$year_g=paste(full_nona$Glacier, full_nona$Year, sep="_")

### free variables to allow external correlations
n.s=full_nona$n.s
ph.s=full_nona$ph.s
ph2=full_nona$ph2
ndvi.s=full_nona$ndvi.s
p.s=full_nona$p.s
twi.s=full_nona$twi.s
time2=I(full_nona$twi.s^2)

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
year_g=paste(full_nona$Year, "_", full_nona$Glacier, sep="")


### THE SEM ANALYSIS STARTS HERE

###  STRUCTURE: plants co-vary with fungi and bacteria
## in preliminary trials, we tried quadratic effects of time but  never significant so removed

### MODEL 1: THE ENVIRONMENTAL FEATURES DIRECTLY AFFECT BIODIVERSITY?
psem_all_direct = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_olig=      lmer(olig.q1.l    ~ph.s+    time.s+twi.s+t.s+ndvi.s+sper.q1.l+n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+ndvi.s          +n.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_n= lmer(n.s~time.s*t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_p= lmer(p.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_ph=lmer(ph.s~time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_ndvi=lmer(ndvi.s~time.s+t.s+n.s+twi.s+sper.q1.l+p.s+(1|year_g)+(1|Glacier), data=full_nona),
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
  lmer_sper=lmer(sper.q1.l~ph.s+time.s+twi.s+n.s+t.s+p.s+(1|year_g)+(1|Glacier), data=full_nona),
  bact.q1.l %~~% sper.q1.l,
  fung.q1.l %~~% sper.q1.l
  )

BIC(psem_all_direct)
# 876.085
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

### ALTERNATIVE APPROACH: BIODIVERSITY determines ecosystem attributes?
psem_multifunct = psem(
  lmer_euka.uni=lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_coll=lmer(coll.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_inse=lmer(inse.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_olig=lmer(olig.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_bact=lmer(bact.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_fung=lmer(fung.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_sper=lmer(sper.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_ph=lmer(ph.s~time.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona),
  
  lmer_n=lmer(n.s~time.s*t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+ndvi.s+(1|year_g)+(1|Glacier), data=full_nona),
  lmer_ndvi=lmer(ndvi.s~time.s+t.s+twi.s+euka.uni.q1.l+euka.ani.q1.l+coll.q1.l+inse.q1.l+olig.q1.l+fung.q1.l+sper.q1.l+bact.q1.l+(1|year_g)+(1|Glacier), data=full_nona),
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

BIC(psem_multifunct)
# 881.083

### THIRD MODEL: ECOSYSTEM ATTRIBUTES COVARY  WITH BIODIVERSITY?
### NOTE: I  TRIED TO ADD TIME^2, particularly given its effect on pH, but it increases the BIC of the model so then excluded.
### non-significant quadratic relationships (eg with pH) not included in the model

psem_covary = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
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
############# COMPARISON BETWEEN THE 3 MODELS:


BIC(psem_covary)
# 710.468
BIC(psem_all_direct)
# 876.085
BIC(psem_multifunct)
#880.083


dSep(psem_covary)
# Independ.Claim Test.Type       DF Crit.Value   P.Value 
# 1       n.s ~ ph2 + ...      coef 791.4616  1.0472507 0.2953037 
# 2       p.s ~ ph2 + ...      coef 788.2995  0.4053318 0.6853436 
# 3 sper.q1.l ~ ph2 + ...      coef 600.5215 -1.6539979 0.0986505 
# 4 olig.q1.l ~ ph2 + ...      coef 667.0053 -0.7866772 0.4317503 
fisherC(psem_covary)
# Fisher.C df P.Value
# 9.507  8   0.301
summary(psem_covary)


## check spatial autocorrelation of residuals using Moran's I

res=residuals(psem_covary$lmer_euka.uni)
(moran=eco.correlog(Z=res, XY=full_nona[,33:34], alternative="greater",seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))

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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],  alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34], alternative="greater",seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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
(moran=eco.correlog(Z=res, XY=full_nona[,33:34],alternative="greater", seqvec=c(0,100,500,1000,10000,50000,100000,200000,400000,1000000),method="I", latlon=T, nsim=300))
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

### starting from PSEM COVARY: test the models: without biotic interactions; without effect of micro-habitat, without effect of time.
psem_co_noh = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  n.s%~~%ph.s,
  p.s%~~%ph.s,
  
  n.s%~~%p.s,
  
  ph.s%~~%ph2,
  ndvi.s%~~%ph.s,
  ndvi.s%~~%ph2,
  ndvi.s%~~%n.s,
  ndvi.s%~~%p.s,
  
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

summary(psem_co_noh)
# Fisher's C = 152.163 with P-value = 0 and on 54 degrees of freedom
# AIC      BIC
# 362.163   853.124
# Response method Marginal Conditional
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

psem_co_nobio = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +time.s+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
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
  ndvi.s%~~%fung.q1.l
)
summary(psem_co_nobio)
# Fisher's C = 1194.119 with P-value = 0 and on 64 degrees of freedom
#   AIC      BIC
# 1394.119   1861.701
# Response method Marginal Conditional
# euka.uni.q1.l   none     0.03        0.29
# euka.ani.q1.l   none     0.15        0.39
# coll.q1.l   none     0.19        0.52
# inse.q1.l   none     0.27        0.57
# olig.q1.l   none     0.10        0.52
# bact.q1.l   none     0.09        0.39
# fung.q1.l   none     0.29        0.59
# sper.q1.l   none     0.09        0.45
# ph.s   none     0.05        0.82
# n.s   none     0.26        0.73
# p.s   none     0.13        0.71
# ndvi.s   none     0.30        0.92

psem_co_notime = psem(
  lmer_euka.uni=  lmer(euka.uni.q1.l~ph.s+ph2+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_euka.multi=lmer(euka.ani.q1.l~ph.s+ph2+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_coll=      lmer(coll.q1.l    ~ph.s+ph2+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_inse=      lmer(inse.q1.l    ~ph.s+ph2+twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_olig=      lmer(olig.q1.l    ~ph.s    +twi.s+t.s+sper.q1.l+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_bact=      lmer(bact.q1.l    ~ph.s+ph2+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_fung=      lmer(fung.q1.l    ~ph.s+ph2+twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_sper=      lmer(sper.q1.l    ~ph.s    +twi.s+t.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  
  lmer_ph=lmer(ph.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F), 
  lmer_n=lmer(n.s~time.s*t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_p=lmer(p.s~time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
  lmer_ndvi=lmer(ndvi.s~sper.q1.l+time.s+t.s+twi.s+(1|year_g)+(1|Glacier), data=full_nona, REML=F),
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

summary(psem_co_notime)
# Fisher's C = 260.51 with P-value = 0 and on 24 degrees of freedom
#  AIC      BIC
#  454.510   908.065
# Response method Marginal Conditional
# euka.uni.q1.l   none     0.07        0.29
# euka.ani.q1.l   none     0.14        0.40
# coll.q1.l   none     0.12        0.51
# inse.q1.l   none     0.20        0.57
# olig.q1.l   none     0.10        0.51
# bact.q1.l   none     0.08        0.39
# fung.q1.l   none     0.20        0.59
# sper.q1.l   none     0.06        0.48
# ph.s   none     0.05        0.82
# n.s   none     0.26        0.73
# p.s   none     0.13        0.71
# ndvi.s   none     0.30        0.92


