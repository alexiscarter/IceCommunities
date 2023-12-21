## SEM beta diversity
## In this script there is the code for the structural equation modeling for the beta diversity

## Packages
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(car)
library(MASS)
library(MuMIn)
library(parallel)
library(raster)
library(lavaan)
library(survey)
library(lavaan.survey)
library(ggplot2)
library(gridExtra)

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
beta.tot=read.csv("data/beta.biotic.plot.csv",header=T,stringsAsFactors=F)   # this file includes the beta-diversity values between plots
div=read.table("data/div.clim.chem.csv", sep=",",header=T,stringsAsFactors=F)		# load the table with habitat features and the values of alpha-diversity

div$all.ani.q0=div$coll.q0+div$olig.q0+div$inse.q0+div$euka.ani.q0  # the total diversity of animals, to remove (later) plots without detected animals

####################
## data selection ## - select the markers of interest and filter the table
####################

nrow(beta.tot)					# 18,429 comparisons
length(unique(beta.tot$glacier1))		# 46 glaciers


## remove glaciers / samples with incomplete soil data

div=div[!is.na(div$n),] 

tokeep=unique(c(beta.tot$Plot1[which(beta.tot$Plot1 %in% div$plot==T)],
		beta.tot$Plot2[which(beta.tot$Plot2 %in% div$plot==T)]))

beta.tot=beta.tot[which(beta.tot$Plot1 %in% tokeep & beta.tot$Plot2 %in% tokeep),]

nrow(beta.tot)					# 10,232 comparisons
length(unique(beta.tot$glacier1))		# 32 glaciers


all(beta.tot$Plot1 %in% div$plot==T)		# ok
all(beta.tot$Plot2 %in% div$plot==T)


## remove samples with 0 species

zero_plus=apply(div[,which(names(div) %in% c("sper.q0","bact.q0","fung.q0","all.ani.q0"))],
		MARGIN=1,FUN=function(x){all(x>0)})

table(zero_plus)
	# zero_plus
	# FALSE  TRUE 
	#   230   563

tokeep=div$plot[which(zero_plus==T)]
beta.tot=beta.tot[which(beta.tot$Plot1 %in% tokeep & beta.tot$Plot2 %in% tokeep),]

nrow(beta.tot)			# 5,741 - or 3,958 including euka.uni.q0


# check for square matrices (within each glacier)

comp=tapply(beta.tot$glacier1,INDEX=beta.tot$glacier1,FUN=function(x){length(x)})
samp=tapply(c(beta.tot$Plot1,beta.tot$Plot2),INDEX=rep(beta.tot$glacier,2),FUN=function(x){length(unique(x))})

comp2=unlist(lapply(samp,FUN=function(x){max(cumsum(1:(x-1)))}))	# this equals (samp^2-samp)/2 or unlist(lapply(samp,FUN=function(x){ncol(combn(1:x,2))}))
all(comp==comp2)							# ok

##################
## beta abiotic ##
##################

# Note: apart for geographic distances, all the other have to rely on scaled values

beta.tot$time_diff=NA
beta.tot$beta.geo=NA
beta.tot$beta.microclim=NA
beta.tot$ndvi_diff=NA
beta.tot$ph_diff=NA
beta.tot$beta.npc=NA

time=scale(div$time)				# variables (except geographic coordinates) must be scaled before calculating multivariate distances elsewhere the ones more variables will overly affect them
ph=scale(div$ph)
p=scale(div$p)
n=scale(div$n)
c=scale(div$c)
meanT=scale(div$mean.temp)
twi=scale(div$twi)
ndvi=scale(div$ndvi)

for(i in 1:nrow(beta.tot)){

	sel=which(div$plot %in% c(beta.tot$Plot1[i],beta.tot$Plot2[i]))
	
	beta.tot$time_diff[i]=dist(time[sel])
	beta.tot$beta.microclim[i]=dist(cbind(meanT[sel],twi[sel]))
	beta.tot$ndvi_diff[i]=dist(ndvi[sel])
	beta.tot$ph_diff[i]=dist(ph[sel])
	beta.tot$beta.npc[i]=dist(cbind(n[sel],p[sel],c[sel]))		# CHANGE!! -> c added
	beta.tot$beta.geo[i]=raster::pointDistance(p1=c(div$lon[sel][1],div$lat[sel][1]),p2=c(div$lon[sel][2],div$lat[sel][2]),lonlat=T)
	
}


summary(beta.tot)			# NAs remaining in beta.coll, beta.euka.uni, beta.euka.animals, beta.inse, beta.olig (not to be used)



#########################
## data transformation ##
#########################
### LOG-TRANSFORM VARIABLES (for all of them, the log transformation is the transformation best improving normality)
beta.tot$beta.microclim.lg=log(beta.tot$beta.microclim+0.01)
beta.tot$beta.npc.lg=log(beta.tot$beta.npc)	
beta.tot$ph_diff.lg=log(beta.tot$ph_diff+0.01)
beta.tot$beta.geo.lg=log(beta.tot$beta.geo)
beta.tot$ndvi.diff.lg=log(beta.tot$ndvi_diff+0.01)
beta.tot$time_diff.lg=log(beta.tot$time_diff+0.5)



## CHECK: AN EXAMPLE OF MIXED MODELS CONSIDERING PLANT BETA DIVERSITY

lmer_sper=lmer(beta.sper ~ time_diff.lg+beta.geo.lg+ph_diff.lg+beta.npc.lg+beta.microclim.lg + (1|glacier1)+(1|Plot1) + (1|Plot2) , data=beta.tot, na.action = na.omit)

summary(lmer_sper)

r.squaredGLMM(lmer_sper)
# R2m       R2c
# [1,] 0.07249042 0.5340582

######
## STRUCTURAL EQUATION MODEL USING PSEM:
######

### WE USED A MODEL ASSUMING CO-VARIATION, AS IT WAS THE BEST SUPPORTED IN THE RICHNESS ANALYSIS

psem_covary=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_npc=lmer(beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  # time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.npc.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.npc.lg,
  
  # Biodiv covariation with nutrients
  beta.npc.lg %~~% beta.all.ani,
  beta.npc.lg %~~% beta.fung,
  beta.npc.lg %~~% beta.bact,
  beta.npc.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

# NOTE THAT I REMOVED THE RELATIONSHIPS: TIME-MICROCLIMATE; MICROCLIM-PH

fisherC(psem_covary)
# Fisher.C df P.Value
# 1    4.225  4   0.376

## the coefficients of the PSEM are here:
(s=summary(psem_covary))

###### NOW USE LAVAAN SURVEY TO CALCULATE BIC VALUES

library(lavaan.survey)

lavaan.covary<-'
beta.microclim.lg ~ beta.geo.lg
ph_diff.lg ~ time_diff.lg + beta.geo.lg+beta.microclim.lg
beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg
ndvi.diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg
  
beta.all.ani ~   beta.sper+time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg
beta.bact ~                time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg
beta.fung ~                time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg
beta.sper ~                time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg
  
  time_diff.lg ~~ beta.geo.lg
  ndvi.diff.lg ~~ beta.all.ani
  ndvi.diff.lg ~~ beta.fung
  ndvi.diff.lg ~~ beta.bact
  ndvi.diff.lg ~~ beta.npc.lg
  ndvi.diff.lg ~~ beta.sper 
  ph_diff.lg ~~ beta.npc.lg
  beta.npc.lg ~~ beta.all.ani
  beta.npc.lg ~~ beta.fung
  beta.npc.lg ~~ beta.bact
  beta.npc.lg ~~ beta.sper
  beta.bact ~~ beta.sper
  beta.fung ~~ beta.sper
  beta.bact ~~ beta.all.ani
  beta.fung ~~ beta.all.ani
  beta.fung ~~ beta.bact'


fit_lavaan_covary <- sem(lavaan.covary, data=beta.tot)
survey.design <- svydesign(ids=~date1+date2+glacier1, prob=~1, data=beta.tot)   ## the 2 sites involved in the comparison are used to take into account non-independence. This is different from pSEM because here we cannot use permutation tests. all fit measures are unchanged if these random factors are not included
survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)

fitMeasures(survey_covary, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "SRMR", "BIC"))
# chisq.scaled     df.scaled pvalue.scaled  rmsea.scaled          srmr           bic   
# 0.030         1.000         0.861         0.000         0.001     50158.158

###########################################################################################
### ASSESS THE RELATIVE IMPORTANCE OF DIFFERENT PATHS USING BIC
########################################################################################### 

### rationale: starting from the covariation  model (fitted in lavaan.survey), we build nodels iteratively removing all the paths
### BIC is used to assess the fit of each model
### it is thus possible testing the average importance of paths representing the effects of habitat, time and biotic interactions

formlist=list(
"beta.all.ani~beta.sper+time_diff.lg+beta.microclim.lg+beta.geo.lg+ph_diff.lg",
"beta.bact~time_diff.lg+beta.microclim.lg+beta.geo.lg+ph_diff.lg",
"beta.fung~time_diff.lg+beta.microclim.lg+beta.geo.lg+ph_diff.lg",
"beta.sper~time_diff.lg+beta.microclim.lg+beta.geo.lg+ph_diff.lg",
"beta.npc.lg~~beta.all.ani",
"beta.npc.lg~~beta.fung",
"beta.npc.lg~~beta.bact",
"beta.npc.lg~~beta.sper",
"beta.bact~~beta.sper",
"beta.fung~~beta.sper",
"beta.bact~~beta.all.ani",
"beta.fung~~beta.all.ani",
"beta.fung~~beta.bact",
"ph_diff.lg~time_diff.lg+beta.geo.lg+beta.microclim.lg",
"beta.npc.lg~time_diff.lg+beta.microclim.lg+beta.geo.lg",
"ndvi.diff.lg~time_diff.lg+beta.microclim.lg+beta.geo.lg+ph_diff.lg",
"beta.microclim.lg~beta.geo.lg",
"time_diff.lg~~beta.geo.lg",
"ndvi.diff.lg~~beta.all.ani",
"ndvi.diff.lg~~beta.fung",
"ndvi.diff.lg~~beta.bact",
"ndvi.diff.lg~~beta.npc.lg",
"ndvi.diff.lg~~beta.sper",
"ph_diff.lg~~beta.npc.lg",
"beta.sper~~0*beta.all.ani"
)
  
  
out=vector("list",0)
i=1

## remove iteratively all the paths representing potential effects on biodiversity:
for(i in 1:13){
  print(i)
  split=unlist(strsplit(formlist[[i]],split="\\~{1}"))		# this produces three elements when the formula is a co-variation
  
  if(length(split)==3){						# if this line represents a co-variation, completely remove i-th covariation (by using "0*")
    
    term1=split[1]
    term2=split[3]
    
    formula=formlist
    formula[[i]]=paste0(term1,"~~0*",term2)
    
    fit_lavaan_covary <- sem(paste0(formula,collapse="\n "), data=beta.tot)
    survey.design <- svydesign(ids=~glacier1, prob=~1, data=beta.tot)
    survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)
    bic=BIC(survey_covary)
    out[[length(out)+1]]=c(formlist[[i]],formlist[[i]],bic)
    
  }else{								#  if this line represents a directional effect, iteratively removes each of the independent variables
    
    
    dep=split[1]
    ind=split[2]
    
    var=unlist(strsplit(ind,split="\\+|\\*"))
    
    for(j in 1:length(var)){
      
      jind=paste(var[-j],collapse="+")
      
      formula=formlist
      formula[[i]]=paste(dep,jind,sep="~")
      
      fit_lavaan_covary <- sem(paste0(formula,collapse="\n "), data=beta.tot)
      survey.design <- svydesign(ids=~glacier1, prob=~1, data=beta.tot)
      survey_covary <- lavaan.survey(lavaan.fit=fit_lavaan_covary, survey.design=survey.design)
      bic=BIC(survey_covary)
      out[[length(out)+1]]=c(formlist[[i]],var[j],bic)
    }
  }
}       

res=matrix(unlist(out),length(out),3,byrow=T)
colnames(res)=c("model","removed_variable","bic")
summary(res)
res=data.frame(res)
res$bic=as.numeric(as.character(res$bic))
res$delta_bic=res$bic-50158.158     ### calculate the change in BIC resulting from the removal of each variable. 50158.158 is the BIC of the full mdoel

## for each path, define if it represents biotic effects (b), habitat (h) or time (t)
type=c("Biot", "Time", "Habitat", " Geog", "Habitat", "Time", "Habitat", " Geog", "Habitat", "Time", "Habitat", " Geog", "Habitat", "Time", "Habitat", " Geog", "Habitat", "Habitat", "Habitat", "Habitat", "Habitat", "Biot", "Biot", "Biot", "Biot", "Biot")

######################################################
#  evaluate if the importance of processes is similar between boreal, temperate and tropical
######################################################
### CREATE 3 GEOGRAPHICALLY-RESTRICTED SUBSETS:

bor=subset(div, div$lat>60|div$lat< -46)
bb=unique(bor$Glacier)
temperate=subset(div, div$lat<60&div$lat> -45)
temperate=subset(temperate, temperate$lat>33|temperate$lat< -20)
te=unique(temperate$Glacier)
tropical=subset(div, div$lat<33&div$lat> -20)
tr=unique(tropical$Glacier)

boreal=beta.tot[beta.tot$glacier1 %in% bb, ]
temp=beta.tot [beta.tot$glacier1 %in% te, ]
trop=beta.tot [beta.tot$glacier1 %in% tr, ]

####################################
### ANALYSIS OF GEOGRAPHICAL SUBSETS USING PSEM
####################################

### BOREAL

psem_covary_boreal=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=boreal, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=boreal, na.action = na.omit),
  lmer_npc=lmer(beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=boreal, na.action = na.omit),
  
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=boreal, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=boreal, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=boreal, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=boreal, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=boreal, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  # time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.npc.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.npc.lg,
  
  # Biodiv covariation with nutrients
  beta.npc.lg %~~% beta.all.ani,
  beta.npc.lg %~~% beta.fung,
  beta.npc.lg %~~% beta.bact,
  beta.npc.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)


### temperate

psem_covary_temperate=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=temp, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=temp, na.action = na.omit),
  lmer_npc=lmer(beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=temp, na.action = na.omit),
  
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=temp, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=temp, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=temp, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=temp, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=temp, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  # time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.npc.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.npc.lg,
  
  # Biodiv covariation with nutrients
  beta.npc.lg %~~% beta.all.ani,
  beta.npc.lg %~~% beta.fung,
  beta.npc.lg %~~% beta.bact,
  beta.npc.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

### tropical

psem_covary_tropical=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=trop, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=trop, na.action = na.omit),
  lmer_npc=lmer(beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=trop, na.action = na.omit),
  
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=trop, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=trop, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=trop, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=trop, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=trop, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  # time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.npc.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.npc.lg,
  
  # Biodiv covariation with nutrients
  beta.npc.lg %~~% beta.all.ani,
  beta.npc.lg %~~% beta.fung,
  beta.npc.lg %~~% beta.bact,
  beta.npc.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

### PUT TOGHETHER COEFFICIENTS TO MAKE COMPARISONS

tab_res_psem_all=summary(psem_covary)$coefficients
tab_res_psem_boreal=summary(psem_covary_boreal)$coefficients
tab_res_psem_temperate=summary(psem_covary_temperate)$coefficients
tab_res_psem_tropical=summary(psem_covary_tropical)$coefficients

typology_psem=c("other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "b", "t", "h", "other", "h", "t", "h", "other", "h", "t", "h", "other", "h", "t", "h", "other", "h", "other", "b", "b", "h", "h", "h", "other", "other", "h", "h", "h", "h", "b", "b", "b")

tab_res_psem=cbind(abs(tab_res_psem_all$Std.Estimate), abs(tab_res_psem_boreal$Std.Estimate), abs(tab_res_psem_temperate$Std.Estimate), abs(tab_res_psem_tropical$Std.Estimate), typology_psem)
colnames(tab_res_psem)<-c("All", "Boreal", "Temperate", "Tropical", "Type")
tab_res_psem=tab_res_psem[!tab_res_psem[,5]=="other",]
tab_res_psem=data.frame(tab_res_psem)
tab_res_psem$All=as.numeric(as.character(tab_res_psem$All))
tab_res_psem$Boreal=as.numeric(as.character(tab_res_psem$Boreal))
tab_res_psem$Temperate=as.numeric(as.character(tab_res_psem$Temperate))
tab_res_psem$Tropical=as.numeric(as.character(tab_res_psem$Tropical))

### PLOTS
################################ 
### PLOTS WITH THE STANDARDIZED EFFECTS
################################ 

## all the landscapes
gg_all= ggplot(tab_res_psem, aes(x = Type, y =  All, fill = Type, col=Type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  ylim(0,0.5)+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
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
  ylim(0,0.5)+
  ylab("")+
  xlab("")+
  theme_minimal()

## temperate only
gg_temp=ggplot(tab_res_psem, aes(x = Type, y =  Temperate, fill = Type, col=Type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.5)+
  ylab("")+
  xlab("")+
  theme_minimal()

## boreal only
gg_bor=ggplot(tab_res_psem, aes(x = Type, y =  Boreal, fill = Type, col=Type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c( "#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylim(0,0.5)+
  ylab("")+
  xlab("")+
  theme_minimal()

# PLOT THE GEOGRAPHICALLY-RESTRICTED DATASETS TOGHETHER:
x11(width=5,heigh=3)
grid.arrange(gg_all, gg_bor, gg_temp, gg_trop, nrow = 1)

### PLOT THE RELATIVE IMPORTANCE ASSESSED BY BIC DROP
res_bic=cbind(res, type)
res_bic=res_bic[!res_bic$type==" Geog",]  ## only keep paths representing the effects of biotic variables, habitat or time

x11(width=5,height = 3)
ggplot(res_bic, aes(x = type, y =  delta_bic, fill = type, col=type)) +
  geom_boxplot(alpha=0.6, fatten = NULL, show.legend = FALSE,   outlier.shape = NA) +#fatten avoids plotting the median
  stat_summary(fun.y = median, col="black", geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
               width = 0.75, size = 1, linetype = "solid", show.legend = FALSE)+
  scale_colour_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  scale_fill_manual(values=c("#44bb99",  "#EEDD88","#FFAABB"))+
  geom_point(    position = position_jitter(width = 0.1), size = 3,   shape = 19,  show.legend = FALSE)+
  ylab("")+
  xlab("")+
  theme_minimal()

### EXPORT TABLES OF DATA USED TO BILD FIG 3C AND FIG 3D
write.table(res_bic, "bic_iterative_exclusion_beta.txt", row.names = F)  # to export and explore the results
write.table(tab_res_psem, "psem_standardized_coefficients_beta.txt", row.names = F)  # to export and explore the results


#################################################################
### SCRIPT TO CALCULATE SIGNIFICANCE OF COEFFICIENTS USING PERMUTATIONS
##################################################################

# please note that it requires several hours to be run with 10,000 permutations. Only run it if you want to produce all the significance values

# retrieve observed coefficients (corrected by "/sqrt(1-R2)"), covariances and R2

resp=c("beta.microclim.lg","ph_diff.lg","beta.npc.lg","ndvi.diff.lg","beta.bact","beta.fung","beta.sper", "beta.all.ani")

obs=numeric(nrow(s$coefficients)+nrow(s$R2))
names(obs)=c(paste0(s$coefficients$Response,"|",s$coefficients$Predictor),			# coefficients and covariances
             paste0("R2_",s$R2$Response))

for(i in 1:length(resp)){
  
  sel.coef=grep(paste0("^",resp[i]),s$coefficients$Response)
  sel.cov=grep(paste0("^\\~\\~",resp[i]),s$coefficients$Response)
  sel.r=grep(resp[i],s$R2$Response)
  
  obs[sel.coef]=s$coefficients$Estimate[sel.coef]/sqrt(1-s$R2$Conditional[sel.r])		# coefficients - pseudo-t test
  obs[sel.cov]=s$coefficients$Estimate[sel.cov]						# covariace - to be treated as a typical Mantel (use as-is)
  obs[nrow(s$coefficients)+sel.r]=s$R2$Conditional[sel.r]					# Rsquared - to be treated as a typical Mantel (use as-is)
  
}

#############
## Permute ##
#############

nperm=9                                 # running permutations is very long. we suggest using 10,000 permutations
sink("permutations.txt",append=T)				# open permutations.txt in your wd (and scroll down to the end) to see progress

Ncores=detectCores()-3						# 3 cores used
cl=makeCluster(Ncores,outfile="permutations.txt")

clusterSetRNGStream(cl,1234567)					# set.seed

clusterExport(cl,						# allow parallel to see *data* of interest for the function
              c("nperm",
                "Ncores",
                #"pb",
                "resp",
                "s",
                "beta.tot"))

clusterEvalQ(cl,						# allow parallel to see *packages* of interest for the function
             c(library(piecewiseSEM),
               library(lme4),
               library(lmerTest),
               library(car),
               library(MASS),
               library(MuMIn))) 

start.time=Sys.time()

perm=parLapply(cl,
               1:nperm,
               fun=function(x){
                 
                 # print (to sink file) the % work executed by the first core
                 
                 if(x<(nperm/Ncores)){cat(paste0("Executed: ",round((x/(nperm/Ncores))*100,2),"%"),sep="\n")}
                 
                 if(x==1){
                   
                   # observed values are on the first row - conservative: cf. Legendre et al., 1994 and Hope, 1968
                   # we treated them outside the parLapply to simplify the function (i.e., avoid to repeat the psem() function twice)
                   
                 }else{
                   
                   exp=numeric(nrow(s$coefficients)+nrow(s$R2))
                   names(exp)=c(paste0(s$coefficients$Response,"|",s$coefficients$Predictor),
                                paste0("R2_",s$R2$Response))
                   
                   
                   for(i in 1:length(resp)){
                     
                     data=beta.tot
                     
                     for(j in 1:length(unique(data$glacier1))){						# within glacier permutation
                       
                       sub=data[which(data$glacier1==unique(data$glacier1)[j]),]
                       
                       lev=unique(c(sub$Plot1,sub$Plot2))
                       
                       r=as.integer(factor(sub$Plot1,levels=lev))
                       c=as.integer(factor(sub$Plot2,levels=lev))
                       D=sub[,which(names(sub)==resp[i])]
                       
                       m=matrix(NA,length(lev),length(lev))
                       m[cbind(r,c)]=m[cbind(c,r)]=D
                       
                       samp=sample(nrow(m))
                       m=m[samp,samp]
                       
                       D1=m[cbind(r,c)]
                       
                       # cbind(sub$Plot1,sub$Plot2,D,D1)
                       # just to check (i=1,j=1 - amola_1850) 
                       # ok, values are transferred between samples
                       
                       data[which(data$glacier1==unique(data$glacier1)[j]),which(names(data)==resp[i])]=D1
                       
                     }
                     
                     
                     psem_x = psem(
                       lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=data, na.action = na.omit),
                       
                       # env
                       lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=data, na.action = na.omit),
                       lmer_np=lmer(beta.npc.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=data, na.action = na.omit),
                       
                       lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
                       
                       # biodiv
                       lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
                       lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
                       lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
                       lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
                       
                       # Covariation env
                       time_diff.lg %~~% beta.geo.lg,
                       
                       # covariation between micro-organisms and plants
                       beta.bact %~~% beta.sper,
                       beta.fung %~~% beta.sper,
                       
                       # Biodiv covariation with NDVI
                       ndvi.diff.lg %~~% beta.all.ani,
                       ndvi.diff.lg %~~% beta.fung,
                       ndvi.diff.lg %~~% beta.bact,
                       ndvi.diff.lg %~~% beta.npc.lg,
                       
                       # Covariation with pH
                       ph_diff.lg %~~% beta.npc.lg,
                       
                       # Biodiv covariation with nutrients
                       beta.npc.lg %~~% beta.all.ani,
                       beta.npc.lg %~~% beta.fung,
                       beta.npc.lg %~~% beta.bact,
                       beta.npc.lg %~~% beta.sper,
                       
                       # Covariation between biotic
                       beta.bact %~~% beta.all.ani,
                       beta.fung %~~% beta.all.ani,
                       beta.fung %~~% beta.bact
                     )
                     
                     
                     s_x=summary(psem_x)
                     
                     sel.coef=grep(paste0("^",resp[i]),s_x$coefficients$Response)
                     sel.cov=grep(paste0("^\\~\\~",resp[i]),s_x$coefficients$Response)
                     sel.r=grep(resp[i],s_x$R2$Response)
                     
                     exp[sel.coef]=s_x$coefficients$Estimate[sel.coef]/sqrt(1-s_x$R2$Conditional[sel.r])
                     exp[sel.cov]=s_x$coefficients$Estimate[sel.cov]
                     exp[nrow(s$coefficients)+sel.r]=s_x$R2$Conditional[sel.r]
                     
                   }
                   
                   
                   return(exp)
                   
                 }
                 
               })

end.time=Sys.time()
end.time-start.time

stopCluster(cl)
sink()

perm[[1]]=obs

names=names(obs)
save(s,names,perm,file="permutations_out.RData")



## Produce usual summary - can be run in a second moment

print(load("permutations_out.RData"))

print(load("permutations_out_npc.RData"))
# [1] "s"     "names" "perm"

perm=matrix(unlist(perm),nperm,length(perm[[1]]),byrow=T)		# list to matrix

pval=apply(perm,MARGIN=2,FUN=function(x){length(which(abs(x)>=abs(x[1])))/nperm})

class=cut(pval,breaks=c(1,0.1,0.05,0.01,0.001,0),labels=c("***","**","*",".",""),include.lowest=T)


# check and assemble

resp=unlist(lapply(names[-grep("^R2",names)],FUN=function(x){unlist(strsplit(x,split="\\|"))[1]}))
pred=unlist(lapply(names[-grep("^R2",names)],FUN=function(x){unlist(strsplit(x,split="\\|"))[2]}))

all(resp==s$coefficients$Response & pred==s$coefficients$Predictor)		# ok - same order

r2=gsub("^R2\\_","",names[grep("^R2",names)])

all(r2==s$R2$Response)								# ok - same order


sum.coef=s$coefficients
sum.R2=s$R2

sum.coef$P_ours=paste0(round(pval[-grep("^R2",names)],3)," ",class[-grep("^R2",names)])
sum.R2$P_ours=paste0(round(pval[grep("^R2",names)],3)," ",class[grep("^R2",names)])

sum.coef
sum.R2
