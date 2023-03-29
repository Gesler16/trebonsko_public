# Loading packages
library(phytools)
library(ape)
library(geiger)
library(dplyr)
library(nlme)
library(car)
library(MuMIn)
library(Hmisc)

# Maybe need to load:
library(lmerTest)
library(tseries)
library(lmtest)
library(DescTools)
library(cAIC4)
library(AICcmodavg)


### Temporal model
## Importing dataframes
birds <- read.table("bird_data.csv",
										header = TRUE,sep=",") #abundance per pond per year
water_df <- read.table("water_transp.csv", header = TRUE,sep=",") 
full_df <- read.csv("data/processed_data/full_df.csv",header=TRUE)

# Merging with traits
merged_full <- merge(full_df,traits, by="spp")
merged_full$mass <- log(merged_full$mass)
merged_full$month <- as.factor(merged_full$month)
merged_full$year <- as.factor(merged_full$year)

# Only ponds measured in both months
surf_df <- data.frame(birds$pond,birds$surface)
colnames(surf_df) <- c("pond","surf")
surf_df <- surf_df[!duplicated(surf_df), ]


## Fishpond water transparency SD
water_df2 <- merge(water_df,surf_df, by ="pond" )
season_sd <- data.frame(year=1,month=2,sd=3, pond_no=1, med_sa=1, mean_sa=1)

month1 <- c(5,7)
years <- year1

for (y in years){
	subs_sd <- subset(water_df2, year==y)
	for (m in month1){
		subs_sd2 <- subset(subs_sd, month==m)
		pond_no <- nrow(subs_sd2)
		sd <- sd(subs_sd2$transp)
		med_sa <- median(subs_sd2$surf)
		mean_sa<- mean(subs_sd2$surf)
		vec <- list(y,m,sd,pond_no, med_sa, mean_sa)
		season_sd[nrow(season_sd)+1,] <- vec
	} 
}


season_sd <- season_sd[-1,]
sd_all <- season_sd

## Creating full dataset for temporal model
# Only species in both months
slope_df <- read.csv("data/processed_data/common_spp_slope_df.csv",header=TRUE)
temp_df <- slope_df
temp_df <- temp_df[,1:3]
colnames(temp_df) <- c("year","month","lm_slope")
temp_df <- arrange(temp_df,month,year)
temp_df$year <- as.factor(temp_df$year)
temp_df$month <- as.factor(temp_df$month)

## Adding mean water transparency per year
mean_tr <- data.frame(month=1,year=1,transp=1)

year1 <- c(2005:2016)
month1 <- c(5,7)

for(y in year1){
	year2 <- y
	sub1 <- subset(water_df,year==year2)
	for(m in month1){
		month2 <- m
		sub2 <- subset(sub1, month==month2) 
		transp <- mean(sub2$transp)
		vec <- list(month2,year2,transp)
		mean_tr[nrow(mean_tr)+1,] <- vec
	}
}
mean_tr <- mean_tr[-1,]

mean_tr <- arrange(mean_tr,month,year)

mean_tr <- mean_tr[,3]

#### Adding water transparency SD per year ####
# Already generated in water_transparency
water_sd <-  arrange(sd_all,month,year)
water_sd <- water_sd[,3]

# Adding NAO index
nao <- read.table("nao.csv", header = TRUE,sep=",")
nao <- arrange(nao,month,year)
nao <- nao[,3]


# Adding temperature
tempr <- read.table("temp_data.csv", header = TRUE,sep=",")
tempr <- arrange(tempr,month,year)
tempr <- tempr[,3]

# Merging df
env_df <- cbind(temp_df,mean_tr,water_sd,nao,tempr)
write.csv(env_df,"output/env_df.csv", row.names = F)

#### Temporal model ####
lm_ym <- lm(lm_slope~month+year,env_df)
summary(lm_ym)
anova(lm_ym)
lm_df <- data.frame(lm_ym$coefficients)

write.csv(lm_df,"output/temp_df.csv",row.names=FALSE)

#### Environmental Models ####
# Standardising/normalising continuous variables
env_df$mean_tr <- (env_df$mean_tr-mean(env_df$mean_tr))/sd(env_df$mean_tr)
env_df$nao <- (env_df$nao-mean(env_df$nao))/sd(env_df$nao)
env_df$lm_slope <- (env_df$lm_slope-mean(env_df$lm_slope))/sd(env_df$lm_slope)
env_df$tempr <- (env_df$tempr-mean(env_df$tempr))/sd(env_df$tempr)
env_df$water_sd <- (env_df$water_sd-mean(env_df$water_sd))/sd(env_df$water_sd)

### Effect of environmental variables on DAR slope
temp_mod <- lm(lm_slope~month+mean_tr*water_sd+month:water_sd+month:mean_tr+year+nao,env_df, na.action = "na.fail")
summary(temp_mod)
anova(temp_mod)
temp_mod$coefficients

# Checking for multicolinearity, normality, heteroscedacity
vif(temp_mod) 
plot(temp_mod$residuals)
acf(temp_mod$residuals)
resd <- resid(temp_mod)
ggpubr::ggqqplot(resd)
hist(resd, xlab = "Residuals", main = "")


par(mfrow=c(2,2))
plot(temp_mod)

# Running model comparison
temp_all <- dredge(temp_mod)
temp_all
avg <- model.avg(temp_all,subset=cumsum(weight) <= .95)
avg$sw
confint(avg)
summary(avg)

### Partitioning variation between density and occupancy

year <- c(2005:2016,2005:2016)
month <- c(rep(5,times=12),rep(7,times=12))
dens_mods <- read.csv("data/processed_data/dens_mods.csv", header=T)
dens_mods <- dens_mods[,-c(4:6)]
dens_mods <- cbind(dens_mods,month,year)
dens_mods <- dens_mods[,-c(1:2)]

part_df <- merge(merged_full,dens_mods,no.dups=F,by=c("month","year"))

part_lm <- lm(lm_slope~log(den+1)*log(occ+1), part_df)
summary(part_lm)
par(mfrow=c(2,2))
plot(part_lm)
dev.off()
resd <- resid(part_lm)
ggpubr::ggqqplot(resd)

hist(resd, xlab = "Residuals", main = "")



### Disentangling the relationships of environmental and sampling variables  
## Influence on water transparency SD
season_sd$year <- as.factor(season_sd$year)
season_sd$month <- as.factor(season_sd$month)

cv_mod <- lm(sd~pond_no+med_sa+mean_sa+year+month,season_sd,
             na.action = "na.fail")
summary(cv_mod)
anova(cv_mod)

## Transparency against SA
water_df3 <- water_df2
water_df3$surf <- log(water_df3$surf)
water_df3$transp <- log(water_df3$transp)
water_df3$month <- as.factor(water_df3$month)
water_df3$year <- as.factor(water_df3$year)
watermod <- lm(transp~month*year+surf,water_df2)
summary(watermod)
anova(watermod)


## Bird density against pond water transparency
abu_pond <- merge(birds,water_df3, by=c("pond","year","month"))
abu_pond$total <- apply(abu_pond[5:54],1,sum)
abu_pond$tot_dens <- abu_pond$total/abu_pond$surface
abu_pond$pond <- as.factor(abu_pond$pond)
abu_pond$year <- as.factor(abu_pond$year)

lm4 <- lme(log(tot_dens+1)~log(transp),data=abu_pond[abu_pond$month==5,],random=~1|pond/year)
summary(lm4)
anova(lm4)

lm4 <- lme(log(tot_dens+1)~log(transp),data=abu_pond[abu_pond$month==7,],random=~1|pond/year)
summary(lm4)
anova(lm4)


# Bird density against pond surface area
lm4 <- lme(log(tot_dens+1)~log(surface),data=abu_pond[abu_pond$month==5,],random=~1|pond/year)
summary(lm4)
anova(lm4)

lm4 <- lme(log(tot_dens+1)~log(surface),data=abu_pond[abu_pond$month==7,],random=~1|pond/year)
summary(lm4)
anova(lm4)




#### PGLS ####
## Preparing dataset for PGLS
### Dataframe for models using slope instead of density or occupancy as response variable

com_spp <- c(unique(merged_full$spp))
months <- c(5,7)

slope_df <- data.frame(spp="first",month=1,slope=1)

for (i in com_spp){
	subset1 <- merged_full[merged_full$spp==i,]
	for (j in months){
		subset2 <- subset1[subset1$month==j,]
		slope1 <- summary(lm(occ~den,subset2))
		slope <- slope1$coefficients[2]
		vec1 <- c(i,j,slope)
		slope_df <- rbind(slope_df,vec1)
	} 
}
slope_df <- slope_df[-1,]

write.csv(slope_df,"data/processed_data/slope_df.csv",row.names=FALSE)


# Merging new slope df with traits
slope_full <- merge(slope_df,traits, by="spp")
slope_full$mass <- log(slope_full$mass)
slope_full$slope <- as.numeric(slope_full$slope)
slope_full$month <- as.factor(slope_full$month)
slope_full$phylo <- slope_full$spp

slope_full <- slope_full[-c(39:40),] # Removing outlier spp (Larus cachinnans; see main text)


# Standardising continuous variables
slope_full$slope <- (slope_full$slope-mean(slope_full$slope))/sd(slope_full$slope)
slope_full$mass <- (slope_full$mass-mean(slope_full$mass))/sd(slope_full$mass)
slope_full$latitude <- (slope_full$latitude-mean(slope_full$latitude))/sd(slope_full$latitude)
slope_full$hssi <- (slope_full$hssi-mean(slope_full$hssi))/sd(slope_full$hssi)
slope_full$dssi <- (slope_full$dssi-mean(slope_full$dssi))/sd(slope_full$dssi)

## Phylogenetic correction 
# Preparing final tree
trees <- read.nexus("output.nex")
cons_tree <- ls.consensus(trees)

# Rename in consensus tree to acronyms and match them to dataset

acro_spp <- c("cyg_olo","ans_ans","buc_cla","tad_tad","ayt_fer","ayt_ful","ayt_mar","ayt_nyr","net_ruf",
							"ana_pla","ana_cre","ana_str","ana_cly","ana_que","ard_cin","cas_alb","nyc_nyc","cic_nig",
							"cic_cic","pha_car","gal_gal","tri_och","tri_ery","tri_tot","tri_gla","act_hyp","phi_pug",
							"ste_hir","chl_hyb","chl_nig","lar_mel","lar_can","lar_arg","lar_cac","lar_rid","lar_min",
							"van_van","cha_dub","tac_ruf","pod_nig","pod_cri","ful_atr","gal_chl","por_por","ral_aqu",
							"pan_bia","loc_nae","loc_flu","acr_sch","acr_pal","acr_sci","acr_aru","rem_pen","lus_sve",
							"mot_alb","mot_fla","mot_cin","emb_sch","alc_ath","hal_alb","cir_aer")

for(i in 1:length(acro_spp)){
	cons_tree$tip.label[i] <- acro_spp[i]
}

dataset_sp <- unique(merged_full$spp)

# Pruning tree
missing <- name.check(cons_tree, slope_full,slope_full$spp)
to_drop <- missing$tree_not_data 

missing <- name.check(cons_tree, merged_full,merged_full$spp)
to_drop <- missing$tree_not_data 


for (i in 1:length(to_drop)){
	cons_tree <- drop.tip(cons_tree,to_drop[i])
}


# Make brownian correlation matrix
bm.mat <- corBrownian(1,phy=cons_tree,form=~phylo)


# adding phylo column
slope_full$phylo <- slope_full$spp
merged_full$phylo <- merged_full$spp



## May PGLS Model
slope_may <- slope_full[slope_full$month==5,]
fullmod <-gls(slope~diet+hssi+dssi+eu_trend+latitude+mass+breeding+migration,correlation=bm.mat,
							data=slope_may,na.action="na.fail",control = list(singular.ok = TRUE))

# Checking for multicolinearity, normality, heteroscedacity
plot(fullmod)
vif(fullmod, type="predictor")
plot(fullmod$residuals)
acf(fullmod$residuals)
resd <- resid(fullmod)
ggpubr::ggqqplot(resd)

summary(fullmod)

## Model selection
# Full model dredging and averaging, getting SWs
mod_sel <- dredge(fullmod, trace=2)
avg <- model.avg(mod_sel)
avg <- model.avg(mod_sel,subset=cumsum(weight) <= .95)
avg$sw
confint(avg)
summary(avg)

## July PGLS Model
slope_july <- slope_full[slope_full$month==7,]
fullmod <-gls(slope~diet+hssi+dssi+eu_trend+latitude+mass+breeding+migration,correlation=bm.mat,
							data=slope_july,na.action="na.fail",control = list(singular.ok = TRUE))

# Checking for multicolinearity, normality, heteroscedacity
plot(fullmod)
vif(fullmod)
plot(fullmod$residuals)
acf(fullmod$residuals)
resd <- resid(fullmod)
ggpubr::ggqqplot(resd)

## Model selection
# Full model dredging and averaging, getting SWs
mod_sel <- dredge(fullmod, trace=2)
avg <- model.avg(mod_sel)
avg <- model.avg(mod_sel,subset=cumsum(weight) <= .95)
avg$sw
confint(avg)
summary(avg)




## All variables per spp per sampling period
# Code
birds3 <- merge(birds,water_df,by = c("pond","month","year"))
want2 <- c("pond","month","year","surface",com_spp,"transp")
want2 <- select(birds3, any_of(want2))
spp_names <- com_spp
numbrs <- c(5:35) #spp columns only
month1 <- c(5,7)
years <- c(2005:2016)

pond_spp <- data.frame(month=1,year=1,pond_no=1,mean_sa=1,spp="a",w_sd=1, w_mean=1)

for (y in years){
	ponds1 <- want2
	ponds1 <- ponds1[ponds1$year==y,]
	for (m in month1){
		ponds2 <- ponds1[ponds1$month==m,]
		for (i in numbrs){
			ponds3 <- ponds2[,c(1:4,i,36)]
			colnames(ponds3) <- c("pond","year","month","surface","spp","transp")
			ponds3 <- ponds3[ponds3$spp>0,]
			mean_sa <- mean(ponds3$surface)
			pond_no <- nrow(ponds3)/134
			name1 <- spp_names[i-4]
			sd <- sd(ponds3$transp)
			mean <- mean(ponds3$transp)
			vec <- list(m,y,pond_no, mean_sa, name1,sd,mean)
			pond_spp[nrow(pond_spp)+1,] <- vec
		}
	} 
}

pond_spp <- pond_spp[-1,]
pond_spp <- pond_spp[pond_spp$pond_no>0,]

## Add density for each spp to include in the analysis
# Common spp in each year
years <- c(2005:2016)
months <- c(5,7)

common_spp <- data.frame(year=1,spp="one")
for(y in years){
	subset1 <- full_df[full_df$year==y,]
	vec1 <- subset1$spp[subset1$month==5]
	vec2 <- subset1$spp[subset1$month==7]
	com_spp2 <- unique(vec1[vec1 %in% vec2])
	year_no <- length(com_spp2)
	year1 <- c(rep(y,time=year_no))
	df1 <- data.frame(year=year1,spp=com_spp2)
	common_spp <- merge(common_spp,df1, all=T)
}
common_spp <- common_spp[-1,]		
merge_dens <-  full_df[,c(1,2,4,5)]
merge_dens <- merge(merge_dens,common_spp, by=c("spp","year"))
env_spp <- merge(pond_spp,merge_dens,by = c("spp","month","year"))


## Models
env_spp$month <- as.factor(env_spp$month)
env_spp$year <- as.factor(env_spp$year)

# Occupancy 
env_spp$pond_no <- log(env_spp$pond_no)
pond_mod <- lme(pond_no~month,data=env_spp,random=~1|spp/year)
summary(pond_mod)
anova(pond_mod)

# Density
env_spp$den <- log(env_spp$den+1)
pond_mod <- lme(den~month,data=env_spp,random=~1|spp/year)
summary(pond_mod)
anova(pond_mod)

# Mean Surface Area
env_spp$mean_sa <- log(env_spp$mean_sa)
pond_mod <- lme(mean_sa~month,data=env_spp,random=~1|spp/year)
summary(pond_mod)
anova(pond_mod)

# Mean Water Transparency
env_spp$w_mean <- log(env_spp$w_mean)
pond_mod <- lme(w_mean~month,data=env_spp,random=~1|spp/year)
summary(pond_mod)
anova(pond_mod)

# Water Transparency SD
pond_mod <- lme(w_sd~month,data=env_spp,random=~1|spp/year,na.action = "na.omit")
summary(pond_mod)
anova(pond_mod)

plot(pond_mod)
resd <- resid(pond_mod)
ggpubr::ggqqplot(resd)
hist(resd, xlab = "Residuals", main = "")


####Does the abundance of diff spp vary between years to a different extent?####
may_abund <- abund_spp[abund_spp$month==5,]
may_abund <- may_abund %>% 
	group_by(spp) %>% 
		summarise(sd= sd(abu), mean=mean(abu))
may_abund$cv <- may_abund$sd/may_abund$mean

july_abund <- abund_spp[abund_spp$month==7,]
july_abund <- july_abund %>% 
	group_by(spp) %>% 
		summarise(sd= sd(abu), mean=mean(abu))
july_abund$cv <- july_abund$sd/july_abund$mean


#### Mean total abundance and occupancy per year ####
## Abundance
# All spp
sum(abund_spp$abu[abund_spp$month==5])/12
sum(abund_spp$abu[abund_spp$month==7])/12

# Only analysed
abu_sum_comm <- merge(abund_spp,common_spp,by=c("spp", "year"))
sum(abu_sum_comm$abu[abu_sum_comm$month==5])/12
sum(abu_sum_comm$abu[abu_sum_comm$month==7])/12

## Occupancy
# All spp
nrow(abund_spp[abund_spp$month==5,])/12 # 31.5
sum(abund_spp$occ[abund_spp$month==5])/31.5/12
nrow(abund_spp[abund_spp$month==7,])/12 # 28.75
sum(abund_spp$occ[abund_spp$month==7])/28.75/12

# Only analysed
abu_sum_comm <- merge(abund_spp,common_spp,by=c("spp", "year"))
nrow(abu_sum_comm[abu_sum_comm$month==5,])/12 # 26.75
sum(abu_sum_comm$occ[abu_sum_comm$month==5])/26.75/12
nrow(abu_sum_comm[abu_sum_comm$month==7,])/12 # 26.75
sum(abu_sum_comm$occ[abu_sum_comm$month==7])/26.75/12





		

#### Variation in mean total abundance and occupancy annually ####	

pond_mod <- lm(abu~year,abund_sums[abund_sums$month==5,])
pond_mod <- lm(abu~year,abund_sums[abund_sums$month==7,])
pond_mod <- lme(abu~year,abund_sums[abund_sums$month==5,],random=~1|spp_no)
pond_mod <- lme(abu~year,abund_sums[abund_sums$month==7,],random=~1|spp_no)
pond_mod <- lme(abu~year,data=abund_sums,random=~1|month)


pond_mod <- lm(occ~year,occ_sums[occ_sums$month==5,])
pond_mod <- lm(occ~year,occ_sums[occ_sums$month==7,])
pond_mod <- lme(occ~year,occ_sums[occ_sums$month==5,],random=~1|spp_no)
pond_mod <- lme(occ~year,occ_sums[occ_sums$month==7,],random=~1|spp_no)
pond_mod <- lme(occ~year,data=occ_sums,random=~1|month)




summary(pond_mod)
anova(pond_mod)


#### Does the population of diving and dabbling ducks vary per year? #### 

ducks <- c("net_ruf","ayt_fer","ayt_ful","ana_cly","ana_cre",
					 "ana_pla","ana_que","ana_str")

may_ducks <- abund_spp[abund_spp$month==5,]
may_ducks <- filter(may_ducks,spp %in% ducks)
pond_mod <- lme(log10(occ)~year,may_ducks,random=~1|spp)
pond_mod <- lme(log10(abu)~year,may_ducks,random=~1|spp)

july_ducks <- abund_spp[abund_spp$month==7,]
july_ducks <- filter(july_ducks,spp %in% ducks)
pond_mod <- lme(log10(occ)~year,july_ducks,random=~1|spp)
pond_mod <- lme(log10(abu)~year,july_ducks,random=~1|spp)

plot(pond_mod)
resd <- resid(pond_mod)
ggpubr::ggqqplot(resd)
hist(resd, xlab = "Residuals", main = "")

summary(pond_mod)
anova(pond_mod)


