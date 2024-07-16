# Loading packages
library(ggplot2)
library(ggsignif)

# Figure 2
slopes <- read.csv("data/processed_data/common_spp_slope_df.csv", header=T)
years_axis <- c(2005:2016)

ggplot(slopes[slopes$month==5,],aes(year,slope))+
	geom_point(size=1.5, color = "#9348D9")+
	geom_smooth(method="lm",size=1, aes(color = "May"),weight = se,se=FALSE)+
	geom_point(data=slopes[slopes$month==7,],size=1.5, color = "#03723B",shape=24)+
	geom_smooth(data=slopes[slopes$month==7,],method="lm",weight = se,size=1, aes(color = "July"),se=FALSE)+
	labs(y = "Slope", x="Year")+
	scale_color_manual(name="Month",breaks=c("May", "July"),values=c('May'='#9348D9', 'July'='#03723B'))+
	theme_classic()+
	scale_x_continuous(labels=as.character(years_axis),breaks=years_axis)+
	theme(plot.title=element_text(hjust=0.5, size=12), axis.title = element_text(size=11), 
				axis.text=element_text(size=10), plot.background = element_blank(), panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())

ggsave(filename = "output/figures/common_spp_slopes.jpeg",width = 1920, height = 1080, units = "px",
			 dpi=320,bg = "white")

# R^2 values of DAR slopes in figure
summary(lm(slope~year, com_spp_slope[com_spp_slope$month==5,])) #May
summary(lm(slope~year, com_spp_slope[com_spp_slope$month==7,])) #July

# Figure 3
# Ponds vector
ponds <- unique(birds$pond)

# Water transparency annual mean for each pond per month 
ponds_wt <- data.frame(pond=1,month=2,mean=2)
month1 <- c(5,7)

for(p in ponds){
	p1 <- subset(water_df,pond==p)
	for(m in month1){
		month <- m
		p2 <- subset(p1, month==m)
		mean <- mean(p2$transp)
		vec <- list(p,month,mean)
		ponds_wt[nrow(ponds_wt)+1,] <- vec
	}
}
ponds_wt <- ponds_wt[-1,]

# Figure code
ponds_wt$month <- as.factor(ponds_wt$month)
monthlabs <- c("May", "July")
ggplot(ponds_wt, aes(x = month, y = mean, fill = month)) + 
	stat_boxplot(geom ='errorbar',width=0.5)+geom_boxplot()+
	geom_signif(comparisons = list(c("5","7")),map_signif_level = T)+
	theme_classic()+ scale_fill_manual(values=c("#9348D9","#03723B"))+
	theme_classic()+labs(y = "Fishpond Trasparency (cm)", x="Month")+ scale_x_discrete(labels= monthlabs)+
	theme(legend.position="none",plot.title=element_text(hjust=0.5, size=14), 
				axis.title = element_text(size=13), axis.text=element_text(size=12), plot.background = element_blank(), panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())

ggsave(filename = "output/figures/water_month.jpeg", width = 1920, height = 1080, units = "px",
			 dpi=320,bg = "white")

# Statistical values for boxplot
t.test(ponds_wt[ponds_wt$month==5,3],ponds_wt[ponds_wt$month==7,3],paired=T)

## Figure 5
#### Species moving around the interspecific DAR in different years/months ####
## 2005 and 2016 chosen as the first and last years of the dataset
# May
pdata <- merge(DO_5_2005[,-4],DO_5_2016[,-4],by=c("spp","year","den","occ"), all=T)
pdata$den <- log(pdata$den+1)
pdata$occ <- log(pdata$occ+1)
pdata$year <- as.factor(pdata$year)

ggplot(pdata,aes(x=occ,y=den))+
	geom_point(aes(color=year),size=2)+
	geom_path(aes(group=spp), arrow=arrow(length=unit(0.2,"cm")),colour="gray58", size=0.5) +
	geom_smooth(data=pdata[pdata$year==2005,],method="lm",color="steelblue3",size=0.8,se=FALSE)+
	geom_smooth(data=pdata[pdata$year==2016,],method="lm",color="brown3",size=0.8,se=FALSE)+
	labs(y ="Density", x="Occupancy")+theme_classic()+
	scale_colour_manual(name="Year", values=c("steelblue3","brown3"))+
	scale_x_continuous(limits=c(0.5,5),breaks=seq(0.5,5,0.5))+
	scale_y_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.3))+geom_text(x=0.6,y=1.4, label="A)", size=5)

ggsave(filename = "output/figures/may_comp.jpeg",width = 1920, height = 1080, units = "px",
			 dpi=320,bg = "white")

# R^2 values of DAR slopes in figure
summary(lm(den~occ, pdata[pdata$year==2005,])) #2005
summary(lm(den~occ, pdata[pdata$year==2016,])) #2016

## July
pdata <- merge(DO_7_2005[,-4],DO_7_2016[,-4],by=c("spp","year","den","occ"), all=T)
pdata$den <- log(pdata$den+1)
pdata$occ <- log(pdata$occ+1)
pdata$year <- as.factor(pdata$year)

ggplot(pdata,aes(x=occ,y=den))+
	geom_point(aes(color=year),size=2)+
	geom_path(aes(group=spp), arrow=arrow(length=unit(0.2,"cm")),colour="gray58", size=0.5) +
	geom_smooth(data=pdata[pdata$year==2005,],method="lm",color="steelblue3",size=0.8,se=FALSE)+
	geom_smooth(data=pdata[pdata$year==2016,],method="lm",color="brown3",size=0.8,se=FALSE)+
	labs(y ="Density", x="Occupancy")+theme_classic()+
	scale_colour_manual(name="Year", values=c("steelblue3","brown3"))+
	scale_x_continuous(limits=c(0.5,5),breaks=seq(0.5,5,0.5))+
	scale_y_continuous(limits=c(0,1.5),breaks=seq(0,1.5,0.3))+geom_text(x=0.6,y=1.4, label="B)", size=5)


ggsave(filename = "output/figures/july_comp.jpeg",width = 1920, height = 1080, units = "px",
			 dpi=320,bg = "white")

# R^2 values of DAR slopes in figure
summary(lm(den~occ, pdata[pdata$year==2005,])) #2005
summary(lm(den~occ, pdata[pdata$year==2016,])) #2016


