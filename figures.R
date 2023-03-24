library(ggplot2)


# Figure 1 - Final
slopes <- read.csv("data/processed_data/common_spp_slope_df.csv", header=T)
years_axis <- c(2005:2016)

ggplot(slopes[slopes$month==5,],aes(year,slope))+
	geom_point(size=1.5, color = "#9348D9")+
	geom_smooth(method="lm",size=1, aes(color = "May"),se=FALSE)+
	geom_point(data=slopes[slopes$month==7,],size=1.5, color = "#03723B",shape=24)+
	geom_smooth(data=slopes[slopes$month==7,],method="lm",size=1, aes(color = "July"),se=FALSE)+
	labs(y = "Slope", x="Year")+
	scale_color_manual(name="Month",breaks=c("May", "July"),values=c('May'='#9348D9', 'July'='#03723B'))+
	theme_classic()+
	scale_x_continuous(labels=as.character(years_axis),breaks=years_axis)+
	theme(plot.title=element_text(hjust=0.5, size=12), axis.title = element_text(size=11), 
				axis.text=element_text(size=10), plot.background = element_blank(), panel.grid.major = element_blank(),
				panel.grid.minor = element_blank())

ggsave(filename = "output/figures/common_spp_slopes.jpeg",width = 1920, height = 1080, units = "px",
			 dpi=320,bg = "white")
