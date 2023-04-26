## SEM beta diversity
## In this script there is the code for the structural equation modeling for the beta diversity

## Rationale ##
# We aim at comparing distance matrices, so the assumption of independence between samples is violated -> Mantel / Partial Mantel / MRM

# Given the complex structure of the model (see below), we'll have several responses and several R2 (marginal and conditional) - one for each response

# We're going to permute, one at a time, all the responses, by keeping fixed the remaining model structure, 
#	and inspecting the p-values under the null hypothesis of no effect of the predictors on the given response

# We're going to generate p-values for:
	# coefficients: pseudo-t test (b/sqrt(1-r2c))
	# covariances: Mantel
	# R2c:Mantel

# The typical permutation scheme for both Mantel and MRM is the row and column shuffling, 
	
	# this reallocates the observations between samples (name shuffling), while keeping fixed the distance structures
	# e.g., if sampleA has d=c(3,5,9), after shuffling will see those values associated to another sample, but the values move as a whole through the matrix

# In turn this involves that some combinations of distances (those never observed) are never generated during permutation
	# i.e., completely random sampling of the matrix values usually leads to low and low Pvalues

# Estimated distances using Sorensen for all the combination of samples - but "0motus vs 0motus" produced NA
	
	# This generated asymmetric and incomplete matrices (i.e., not-square and with lots of NAs) 
	# This involves we cannot generate row-and-columnns permutations
		# with a symmetrix matrix all samples have the same number of pairwise comparisons
		# with an asymmetric matrix (e.g., inter-species only comparisons) all samples from the same species have the same number of pairwise comparisons

	# Simply transforming asymmetric and incomplete matrices to vector and removing NAs equals performing a random sample, that usually boosts P-values (see above)

# Two approaches to solve the problem:

	# Remove from the pairwise comparisons samples with 0 motus and produce square matrices - but saying it with Alexis "absences may be biologically meaningful"
	# Remove NAs after permutation (so we need a complete / square matrix - with NAs when needed - to run the analysis)
	

# Still, permutations may answer to slightly different questions

	# In the first case, we are simply looking for the lack of associations between distance matrices (null hypothesis)
	# in the second case, we look for both:
	
		# the lack of association
		# the lack of effect on co-absence (that may be biologically due to the value of the variable(s), rather than to their distance)
	
# Additionally: *when combining several responses and removing NAs, may generate datasets with varying number of rows, unless we shuffle all the bio part together*


## Example

	# 5 sites with different species composition (BD) and at different distances along an equally-spaced transect (GD)
		
		# two sites (4 and 5) with zero species (e.g., recently deglaciated / close to glacier front)
		

	# GD			BD
	#     [1] [2] [3] [4]	    [1] [2] [3] [4]
	# [2] 1.0		[2] 0.2
	# [3] 2.0 1.0		[3] 0.4 0.2
	# [4] 3.0 2.0 1.0	[4] 1.0 1.0 1.0
	# [5] 4.0 3.0 2.0 1.0	[5] 1.0 1.0 1.0  NA
	

# First approach
	# GD 1.0  2.0  1.0
	# BD 0.2  0.4  0.2
	
# Second approach
	# GD 1.0  2.0  3.0  4.0  1.0  2.0  3.0  1.0  2.0  1.0		# sort of threshold effect - non-dependence on GB
	# BD 0.2  0.4  1.0  1.0  0.2  1.0  1.0  1.0  1.0  NA

## Load packages and data
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(car)
library(MASS)
library(MuMIn)
library(parallel)
library(raster)

beta.tot=read.csv("data/beta.biotic.plot.csv",header=T,stringsAsFactors=F)
all_data=read.delim("data/div.clim.chem.csv", sep=";")
div=all_data[!is.na(all_data$lg_n),] # remove rows with NA

div$all.ani.q0=div$coll.q0+div$olig.q0+div$inse.q0+div$euka.ani.q0

####################
## data selection ## - select the markers of interest and filter the table
####################

nrow(beta.tot)					# 18,429 comparisons
length(unique(beta.tot$glacier1))		# 46 glaciers

## remove glaciers / samples with incomplete soil data
tokeep=unique(c(beta.tot$Plot1[which(beta.tot$Plot1 %in% div$uniqPlot==T)],
		beta.tot$Plot2[which(beta.tot$Plot2 %in% div$uniqPlot==T)]))

beta.tot=beta.tot[which(beta.tot$Plot1 %in% tokeep & beta.tot$Plot2 %in% tokeep),]

nrow(beta.tot)					# 10,232 comparisons
length(unique(beta.tot$glacier1))		# 32 glaciers


all(beta.tot$Plot1 %in% div$uniqPlot==T)		# ok
all(beta.tot$Plot2 %in% div$uniqPlot==T)

## remove samples with 0 species
zero_plus=apply(div[,which(names(div) %in% c("sper.q0","bact.q0","fung.q0","all.ani.q0"))],		# eventually add ,"euka.uni.q0"
		MARGIN=1,FUN=function(x){all(x>0)})

table(zero_plus)
	# zero_plus
	# FALSE  TRUE 
	#   230   563

tokeep=div$uniqPlot[which(zero_plus==T)]
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

# Note: log-transform time, ndvi, n, p and twi before calculating distances / differences
	# in the simpler case, this involves we are working with log(ratios) (log(time1)-log(time2) = log(time1/time2))

beta.tot$time_diff=NA
beta.tot$beta.geo=NA
beta.tot$beta.microclim=NA
beta.tot$ndvi_diff=NA
beta.tot$ph_diff=NA
beta.tot$beta.np=NA
beta.tot$beta.soil=NA

time=scale(div$time)				# variables (except geographic coordinates) must be scaled before calculating multivariate distances elsewhere the ones more variables will overly affect them
ph=scale(div$ph)
p=scale(div$p)
n=scale(div$n)
meanT=scale(div$meanT_new)
twi=scale(div$twi)
ndvi=scale(div$ndvi)

for(i in 1:nrow(beta.tot)){

	sel=which(div$uniqPlot %in% c(beta.tot$Plot1[i],beta.tot$Plot2[i]))
	
	beta.tot$time_diff[i]=dist(time[sel])
	beta.tot$beta.microclim[i]=dist(cbind(meanT[sel],twi[sel]))
	beta.tot$ndvi_diff[i]=dist(ndvi[sel])
	beta.tot$ph_diff[i]=dist(ph[sel])
	beta.tot$beta.np[i]=dist(cbind(n[sel],p[sel]))
	beta.tot$beta.soil[i]=dist(cbind(n[sel],p[sel],ph[sel]))

	beta.tot$beta.geo[i]=raster::pointDistance(p1=c(div$lon[sel][1],div$lat[sel][1]),p2=c(div$lon[sel][2],div$lat[sel][2]),lonlat=T)
	
}

summary(beta.tot)			# NAs remaining in beta.coll, beta.euka.uni, beta.euka.animals, beta.inse, beta.olig (not to be used)

### LOG-TRANSFORM VARIABLES (for all of them, the log transformation is the one allowing to best approach normality)
beta.tot$beta.microclim.lg=log(beta.tot$beta.microclim+0.01)
beta.tot$beta.np.lg=log(beta.tot$beta.np)
beta.tot$ph_diff.lg=log(beta.tot$ph_diff+0.01)
beta.tot$beta.soil.lg=log(beta.tot$beta.soil)
beta.tot$beta.geo.lg=log(beta.tot$beta.geo)
beta.tot$ndvi.diff.lg=log(beta.tot$ndvi_diff+0.01)
beta.tot$time_diff.lg=log(beta.tot$time_diff+0.5)

## the PSEM MODEL ASSUMING CO-VARIATION BETWEEN BIOTIC VARIABLES AND SOIL
psem=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # Biodiv covariation with nutrients
  beta.np.lg %~~% beta.all.ani,
  beta.np.lg %~~% beta.fung,
  beta.np.lg %~~% beta.bact,
  beta.np.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

s=summary(psem)
BIC(psem)
# 458.736

# retrieve observed coefficients (corrected by "/sqrt(1-R2)"), covariances and R2

resp=c("beta.microclim.lg","ph_diff.lg","beta.np.lg","ndvi.diff.lg","beta.bact","beta.fung","beta.sper", "beta.all.ani")

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
## PERMUTATION TEST ##
#############

nperm=9999
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
					  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=data, na.action = na.omit),

					  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=data, na.action = na.omit),
					  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),

					  # biodiv
					  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
					  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
					  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),
					  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=data, na.action = na.omit),

					  # Covariation env
					  time_diff.lg %~~% beta.geo.lg,
					  time_diff.lg %~~% beta.microclim.lg,

					  # covariation between micro-organisms and plants
					  beta.bact %~~% beta.sper,
					  beta.fung %~~% beta.sper,

					  # Biodiv covariation with NDVI
					  ndvi.diff.lg %~~% beta.all.ani,
					  ndvi.diff.lg %~~% beta.fung,
					  ndvi.diff.lg %~~% beta.bact,
					  ndvi.diff.lg %~~% beta.np.lg,

					  # Covariation with pH
					  ph_diff.lg %~~% beta.np.lg,

					  # Biodiv covariation with nutrients
					  beta.np.lg %~~% beta.all.ani,
					  beta.np.lg %~~% beta.fung,
					  beta.np.lg %~~% beta.bact,
					  beta.np.lg %~~% beta.sper,

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

print(load("permutations_out_np.RData"))
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

write.csv(sum.coef,"coefs.csv",row.names=F)
write.csv(sum.R2,"R2.csv",row.names=F)


#######################################################
#######################################################
#######################################################
### MODELS WITHOUT KEY COMPONENTS OF THE SYSTEM

### THE MODEL WITHOUT BIOTIC INTERACTIONS

psem_nobio=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,

  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # Biodiv covariation with nutrients
  beta.np.lg %~~% beta.all.ani,
  beta.np.lg %~~% beta.fung,
  beta.np.lg %~~% beta.bact,
  beta.np.lg %~~% beta.sper
)
summary(psem_nobio)
# Fisher's C = 2187.683 with P-value = 0 and on 12 degrees of freedom

BIC(psem_nobio)
# 2637.763

###############
## THE MODEL WITHOUT THE EFFECTS OF HABITAT ON BIODIVERSITY
###############

psem_nohabitat=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  # ndvi.diff.lg %~~% beta.all.ani,
  # ndvi.diff.lg %~~% beta.fung,
  # ndvi.diff.lg %~~% beta.bact,
   ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # Biodiv covariation with nutrients
  # beta.np.lg %~~% beta.all.ani,
  # beta.np.lg %~~% beta.fung,
  # beta.np.lg %~~% beta.bact,
  # beta.np.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

summary(psem_nohabitat)
# Fisher's C = 657.433 with P-value = 0 and on 30 degrees of freedom


BIC(psem_nohabitat)
# 1046.925


###############
## THE MODEL WITHOUT THE direct EFFECTS OF TIME ON BIODIVERSITY
###############

psem_notime=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # Biodiv covariation with NDVI
  ndvi.diff.lg %~~% beta.all.ani,
  ndvi.diff.lg %~~% beta.fung,
  ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # Biodiv covariation with nutrients
  beta.np.lg %~~% beta.all.ani,
  beta.np.lg %~~% beta.fung,
  beta.np.lg %~~% beta.bact,
  beta.np.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)
summary(psem_notime)
# Fisher's C = 390.453 with P-value = 0 and on 8 degrees of freedom

BIC(psem)
# 458.736

BIC(psem_nobio)
# 2637.763

BIC(psem_nohabitat)
# 1046.925

BIC(psem_notime)
# 814.567


####################################
###  ALTERNATIVE MODELS:
####################################

# NUTRIENT-LED MODEL: SOIL FEATURES AND PRODUCTIVITY DETERMINE BIODIVERSITY (model A in Fig S4)



psem_soil_affects_bio=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  lmer_ndvi=lmer(ndvi.diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + ndvi.diff.lg+beta.np.lg+time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ ndvi.diff.lg+beta.np.lg+time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ ndvi.diff.lg+beta.np.lg+time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ ndvi.diff.lg+beta.np.lg+time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  # covariations with NDVI
  # ndvi.diff.lg %~~% beta.all.ani,
  # ndvi.diff.lg %~~% beta.fung,
  # ndvi.diff.lg %~~% beta.bact,
  ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # # Biodiv covariation with nutrients
  # beta.np.lg %~~% beta.all.ani,
  # beta.np.lg %~~% beta.fung,
  # beta.np.lg %~~% beta.bact,
  # beta.np.lg %~~% beta.sper,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)
summary(psem_soil_affects_bio)
BIC(psem_soil_affects_bio)
# 519.323

# BIODIVERSITY-LED MODEL: BIODIVERSITY AFFECTS SOIL FEATURES AND PRODUCTIVITY. (model B in Fig S4)
### NOTE THAT IN THIS MODEL WE REMOVED THE POSSIBLE EFFECT OF ANIMALS ON NP AND NDVI BECAUSE OF 1) UNLIKELY EFFECT; 2) THE OBSERVED EFFECT IS NEGATIVE (COUNTERINTUITIVE)
psem_biodiv_affects_soil=psem(
  lmer_microclim=lmer(beta.microclim.lg ~ beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # env
  lmer_ph=lmer(ph_diff.lg ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_np=lmer(beta.np.lg ~ beta.sper+beta.bact+beta.fung+time_diff.lg + beta.microclim.lg + beta.geo.lg + (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # lmer_soil=lmer(beta.soil.sc ~ time_diff.sc + beta.microclim.sc + beta.geo.sc +   (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_ndvi=lmer(ndvi.diff.lg ~ beta.sper + beta.bact+beta.fung+time_diff.lg + beta.microclim.lg + beta.geo.lg +ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # biodiv
  lmer_animals=lmer(beta.all.ani ~ beta.sper + time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_bact=   lmer(beta.bact ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_fung=   lmer(beta.fung ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  lmer_sper=   lmer(beta.sper ~ time_diff.lg + beta.microclim.lg + beta.geo.lg + ph_diff.lg+ (1|glacier1), data=beta.tot, na.action = na.omit),
  
  # Covariation env
  time_diff.lg %~~% beta.geo.lg,
  time_diff.lg %~~% beta.microclim.lg,
  
  # covariation between micro-organisms and plants
  beta.bact %~~% beta.sper,
  beta.fung %~~% beta.sper,
  
  #  covariation with NDVI
  ndvi.diff.lg %~~% beta.np.lg,
  
  # Covariation with pH
  ph_diff.lg %~~% beta.np.lg,
  
  # Covariation between biotic
  beta.bact %~~% beta.all.ani,
  beta.fung %~~% beta.all.ani,
  beta.fung %~~% beta.bact
)

BIC(psem)
# 458.736

BIC(psem_biodiv_affects_soil)
# 517.087

BIC(psem_soil_affects_bio)
# 519.323

