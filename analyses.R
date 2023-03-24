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

# Only ponds measured in both months
surf_df <- data.frame(birds$pond,birds$surface)
colnames(surf_df) <- c("pond","surf")
surf_df <- surf_df[!duplicated(surf_df), ]
water_df <- merge(water_df,surf_df, by ="pond" )

## Fishpond water transparency SD
season_sd <- data.frame(year=1,month=2,sd=3, pond_no=1, med_sa=1, mean_sa=1)

month1 <- c(5,7)
years <- year1

for (y in years){
	subs_sd <- subset(water_df, year==y)
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
temp_df <- temp_df[,2:4]
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

write.csv(lm_df,"output/temp_df.csv")

#### Environmental Models ####
# Standardising/normalising continuous variables
env_df$mean_tr <- (env_df$mean_tr-mean(env_df$mean_tr))/sd(env_df$mean_tr)
env_df$nao <- (env_df$nao-mean(env_df$nao))/sd(env_df$nao)
env_df$lm_slope <- (env_df$lm_slope-mean(env_df$lm_slope))/sd(env_df$lm_slope)
env_df$tempr <- (env_df$tempr-mean(env_df$tempr))/sd(env_df$tempr)
env_df$water_sd <- (env_df$water_sd-mean(env_df$water_sd))/sd(env_df$water_sd)

### Environmental model
temp_mod <- lm(lm_slope~month+mean_tr*water_sd+month:water_sd+month:mean_tr+year+nao,env_df, na.action = "na.fail")
summary(temp_mod)
anova(temp_mod)
temp_mod$coefficients


# Checking for multicolinearity, normality, heteroscedacity
vif(temp_mod, type="marginal") 
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


with(env_df, plot(lm_slope,year))

### Partitioning variation between density and occupancy

year <- c(2005:2016,2005:2016)
month <- c(rep(5,times=12),rep(7,times=12))
dens_mods <- read.csv("data/processed_data/dens_mods.csv", header=T)
dens_mods <- dens_mods[,-c(4:6)]
dens_mods <- cbind(dens_mods,month,year)
dens_mods <- dens_mods[,-c(1:2)]

part_df <- merge(merged_full,dens_mods,no.dups=F,by=c("month","year"))

# Simple linear model to test the impact of each variable (no need for random effects since slope is unique
# for each month-year combo and also there is only one spp in each)

# part_lm <- lm(lm_slope~occ, part_df)
# part_lm <- lm(lm_slope~den, part_df) 
# part_lm <- lm(lm_slope~den:occ, part_df) 

summary(lm(den~occ,part_df))
part_lm <- lm(lm_slope~den*occ, part_df)
#part_lm <- lm(lm_slope~den*occ, part_df[part_df$month==7,])
summary(part_lm)
par(mfrow=c(2,2))
plot(part_lm)
dev.off()
resd <- resid(part_lm)
ggpubr::ggqqplot(resd)

hist(resd, xlab = "Residuals", main = "")

#### PGLS ####
## Preparing dataset for PGLS
### Dataframe for models using slope instead of density or occupancy as response variable

# Merging with traits
merged_full <- merge(full_df,traits, by="spp")
merged_full$mass <- log(merged_full$mass)
merged_full$month <- as.factor(merged_full$month)
merged_full$year <- as.factor(merged_full$year)


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


write.csv(slope_df,"data/processed_data/slope_df.csv")


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


slope_may <- slope_full[slope_full$month==5,]
slope_july <- slope_full[slope_full$month==7,]

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
fullmod <-gls(slope~diet+hssi+dssi+eu_trend+latitude+migration+mass+waterbird,correlation=bm.mat,
							data=slope_may,na.action="na.fail",control = list(singular.ok = TRUE))

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

## July PGLS Model
fullmod <-gls(slope~diet+hssi+dssi+eu_trend+latitude+migration+mass+waterbird,correlation=bm.mat,
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




### Environmental Model
