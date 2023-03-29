# Creating folders
dir.create('data')
dir.create('data/processed_data')
dir.create('output')
dir.create('output/figures')

# Loading plugins
library(data.table)
library(dplyr)

### Importing and preparing bird census dataset
birds <- read.table("bird_data.csv",header = TRUE,sep=",")

#Transpose and add column and row names
species_rows <- transpose(birds[,5:52])
rownames(species_rows) <- colnames(birds[,5:52])
rownames1 <- as.character(birds$year)
rownames2 <- as.character(birds$month)
rownames3 <- as.character(birds$pond)
rownames4 <- paste(rownames1,rownames2,rownames3,sep = "_")
colnames(species_rows) <- rownames4

#converting categorical variables to factor
birds2 <- birds
birds2$year <- as.factor(birds2$year)
birds2$month <- as.factor(birds2$month)
birds2$pond <- as.factor(birds2$pond)


### Implementing thresholds
# Dataset abundance threshold function for a threshold of X individuals per spp per census on average
b_thresh <- function(df, abu.thr){
	b_thresh1 <- df[,1:4]
	b_thresh2 <- df[,5:ncol(df)]
	b_thresh2 <- b_thresh2[colSums(b_thresh2)>abu.thr*24]
	thresh_df <- cbind(b_thresh1,b_thresh2)
}
# Extracting for the selected threshold of 2 individuals per census
birds3 <- b_thresh(birds2, 2)

# Getting the subset species list
species <- birds3
species <- species[5:ncol(species)]
species <- as.character(colnames(species))

# Adding the spp that will be retained as its sampling is highly accurate (Anas clypaeta)

species <- append(species, "ana_cly")

# Combined selection
sel1 <- c("pond","surface","year","month", species)
birds3 <- select(birds2, any_of(sel1))

### Regression models for each census
# Function for extracting density and occupancy per species, with option of additional thresholds 

fun.de.oc.df <- function(df,month,year,abu.thr=0, occ.thr=0){
	df1 <- df[df$month==month,]
	df1 <- df1[df1$year==year,]
	birds_den1 <- df1[,1:4]
	birds_thr <- df1[,5:ncol(df1)]
	birds_thr <- birds_thr[colSums(birds_thr)>abu.thr]
	birds_den2 <- birds_thr/birds_den1$surface
	density <- colSums(birds_den2)/colSums(birds_den2!=0)
	occupancy <- colSums(birds_den2!=0)
	de_oc <- data.frame(colnames(birds_den2))
	de_oc <- cbind(de_oc,density,occupancy)
	colnames(de_oc ) <- c("species","den","occ")
	de_oc <- data.frame(de_oc[de_oc$occ>occ.thr,])
}

# Making dataframes per month/year

years <- c(2005:2016)
months <- c(5,7)

for(y in years){
	for(m in months){
		df1 <- fun.de.oc.df(birds3,m,y)
		df1 <- cbind(df1,m,y)
		colnames(df1) <- c("spp","den","occ","month","year")
		assign(paste0("DO_",m,"_",y),df1)
	}
}
# Model list and names
do_mod_names <- c("DO_5_2005","DO_5_2006","DO_5_2007","DO_5_2008","DO_5_2009","DO_5_2010",
								 "DO_5_2011","DO_5_2012","DO_5_2013","DO_5_2014","DO_5_2015","DO_5_2016",
								 "DO_7_2005","DO_7_2006","DO_7_2007","DO_7_2008","DO_7_2009","DO_7_2010",
								 "DO_7_2011","DO_7_2012","DO_7_2013","DO_7_2014","DO_7_2015","DO_7_2016")

do_mod_list <- list(DO_5_2005,DO_5_2006,DO_5_2007,DO_5_2008,DO_5_2009,DO_5_2010,DO_5_2011,DO_5_2012,
									DO_5_2013,DO_5_2014,DO_5_2015,DO_5_2016,DO_7_2005,DO_7_2006,DO_7_2007,DO_7_2008,
									DO_7_2009,DO_7_2010,DO_7_2011,DO_7_2012,DO_7_2013,DO_7_2014,DO_7_2015,DO_7_2016)



# Relationship regression statistics and testing for normality
norm_df <- data.frame( n_spp=1,slope=1,adj_r_rq = 1, p_value = 1,shap.pv = 1) 
for(i in do_mod_list){									#change list
	df <- i
	df$den <- log(df$den+1)		#change variable
	df$occ <- log(df$occ+1)
	df$den <- (df$den-mean(df$den))/sd(df$den)		#change variable
	df$occ <- (df$occ-mean(df$occ))/sd(df$occ)
	lm1 <- lm(occ~den,df)		#change variable
	slope <- lm1$coefficients[2]
	stat.coef  <- summary(lm1)$coefficients
	p_val <- stat.coef[2,4]
	r_sqr <- summary(lm1)$adj.r.squared
	shap <- shapiro.test(lm1$residuals)
	pv <- shap$p.value
	n_spp <- nrow(i)
	list1 <- list(n_spp,slope,r_sqr,p_val,pv)
	norm_df[nrow(norm_df)+1,] <- list1
}

norm_df <- data.frame(norm_df[-1,])
norm_df <- cbind(do_mod_names,norm_df) #change names
colnames(norm_df) <- c("Model","n_spp","lm_slope","adj_r_sq","p_value","shapiro_pv")


# Exporting
write.csv(norm_df, file="data/processed_data/dens_mods.csv",row.names=FALSE)

for (i in 1:length(do_mod_names)){
	write.csv(get(do_mod_names[i]),paste0("data/processed_data/",do_mod_names[i],".csv"),
						 row.names=FALSE)
}



### Preparing full dataset
## Importing data

# Relationship dataframes
do_mod_names <- c("DO_5_2005","DO_5_2006","DO_5_2007","DO_5_2008","DO_5_2009","DO_5_2010",
									"DO_5_2011","DO_5_2012","DO_5_2013","DO_5_2014","DO_5_2015","DO_5_2016",
									"DO_7_2005","DO_7_2006","DO_7_2007","DO_7_2008","DO_7_2009","DO_7_2010",
									"DO_7_2011","DO_7_2012","DO_7_2013","DO_7_2014","DO_7_2015","DO_7_2016")

years <- c(2005:2016,2005:2016)
months <- c(rep(5,times=12),rep(7,times=12))

for(i in 1:length(do_mod_names)) {
	assign(paste0("DO_",months[i],"_",years[i]),read.csv(paste0("data/processed_data/",
																															do_mod_names[i],".csv"), header=TRUE))
}

do_mod_list <- list(DO_5_2005,DO_5_2006,DO_5_2007,DO_5_2008,DO_5_2009,DO_5_2010,DO_5_2011,DO_5_2012,
										DO_5_2013,DO_5_2014,DO_5_2015,DO_5_2016,DO_7_2005,DO_7_2006,DO_7_2007,DO_7_2008,
										DO_7_2009,DO_7_2010,DO_7_2011,DO_7_2012,DO_7_2013,DO_7_2014,DO_7_2015,DO_7_2016)





### Merging datasets
# Adding month and year as variables
full_df <- bind_rows(do_mod_list, .id = "column_label")
full_df <- full_df[,2:6]
full_df$den <- log(full_df$den+1)
full_df$occ <- log(full_df$occ+1)

write.csv(full_df,"data/processed_data/full_df.csv", row.names=F)

# Traits dataframe
traits <- read.csv("species_traits.csv",header=TRUE)
traits$waterbird <- as.factor(traits$waterbird)
traits$eu_trend <- as.factor(traits$eu_trend)
traits$migration <- as.factor(traits$migration)


# Merging with traits
full_df <- read.csv("data/processed_data/full_df.csv",header=TRUE)
merged_full <- merge(full_df,traits, by="spp")
merged_full$mass <- log(merged_full$mass)
merged_full$month <- as.factor(merged_full$month)
merged_full$year <- as.factor(merged_full$year)

# Standardising continuous variables to prepare df for intra-annual slope graph
merged_full$occ <- (merged_full$occ-mean(merged_full$occ))/sd(merged_full$occ)
merged_full$den <- (merged_full$den-mean(merged_full$den))/sd(merged_full$den)
merged_full$mass <- (merged_full$mass-mean(merged_full$mass))/sd(merged_full$mass)
merged_full$latitude <- (merged_full$latitude-mean(merged_full$latitude))/sd(merged_full$latitude)
merged_full$hssi <- (merged_full$hssi-mean(merged_full$hssi))/sd(merged_full$hssi)
merged_full$dssi <- (merged_full$dssi-mean(merged_full$dssi))/sd(merged_full$dssi)

# Species found in both may and july (common spp) 
com_spp_slope <- data.frame(year=1,month=1,slope=1, rsq = 1)
years <- c(2005:2016)
months <- c(5,7)

for(y in years){
	subset1 <- merged_full[merged_full$year==y,]
	vec1 <- subset1$spp[subset1$month==5]
	vec2 <- subset1$spp[subset1$month==7]
	com_spp <- unique(vec1[vec1 %in% vec2])
	subset1 <- subset1[subset1$spp%in% com_spp,]
	for(m in months){
		subset2 <- subset1[subset1$month==m,]
		mod <- lm(occ~den,subset2)
		slope <- as.numeric(mod$coefficients[2])
		rsq <- summary(mod)$adj.r.squared
		vec <- c(y,m,slope, rsq)
		com_spp_slope <-rbind(com_spp_slope,vec)
	}
}
com_spp_slope <- com_spp_slope[-1,]		

write.csv(com_spp_slope,"data/processed_data/common_spp_slope_df.csv",, row.names=F)
