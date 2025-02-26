setwd("/home/j_xiao/data/NCdemo/")
library(raster)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(cowplot)

normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}


gpp_sens_stk <- stack(paste0("./results/FLUXCOM_GPP/FLUXCOM_GPP_vpd_sensmap_alpha0_ERA5_monthly","_window",seq(1,228,by=12)[1:10],".tif"))
ndvi_sens_stk <- stack(paste0("./results/GIMMS_NDVI/GIMMS_NDVI_vpd_sensmap_alpha0_ERA5_monthly","_window",seq(1,228,by=12)[1:10],".tif"))
nirv_sens_stk <- stack(paste0("./results/MODIS_NIRv/NIRv_vpd_sensmap_alpha0_ERA5_monthly","_window",seq(1,228,by=12)[1:10],".tif"))

mask <- shapefile("/home/j_xiao/data/Terrestrial Ecoregions of the World_WWF/mask_shpfile.shp")



#total mean

r <- gpp_sens_stk
r <- ndvi_sens_stk
r <- nirv_sens_stk

r <- abs(r)
r <- mask(r,mask,inverse=T)

mean_gpp <- cellStats(r,"mean",na.rm=T)
sd_gpp <- cellStats(r,"sd",na.rm=T)
mean_ndvi <- cellStats(r,"mean",na.rm=T)
sd_ndvi <- cellStats(r,"sd",na.rm=T)
mean_nirv <- cellStats(r,"mean",na.rm=T)
sd_nirv <- cellStats(r,"sd",na.rm=T)

lower_limit <- mean_gpp - 1.5 * sd_gpp
upper_limit <- mean_gpp + 1.5 * sd_gpp
r_cleaned <- r
r_cleaned[r_cleaned < lower_limit | r_cleaned > upper_limit] <- NA
mean_gpp_cleaned <- cellStats(r_cleaned, "mean", na.rm = TRUE)

df1 <- data.frame(
  window = 1:10,
  mean = mean_gpp_cleaned,
  scaled_mean= normalize(mean_gpp_cleaned),
  vi ="GPP"
)



lower_limit <- mean_ndvi - 1.5 * sd_ndvi
upper_limit <- mean_ndvi + 1.5 * sd_ndvi
r_cleaned <- r
r_cleaned[r_cleaned < lower_limit | r_cleaned > upper_limit] <- NA
mean_ndvi_cleaned <- cellStats(r_cleaned, "mean", na.rm = TRUE)

df2 <- data.frame(
  window = 1:10,
  mean = mean_ndvi_cleaned,
  scaled_mean= normalize(mean_ndvi_cleaned),
  vi ="NDVI"
)

lower_limit <- mean_nirv - 1.5 * sd_nirv
upper_limit <- mean_nirv + 1.5 * sd_nirv
r_cleaned <- r
r_cleaned[r_cleaned < lower_limit | r_cleaned > upper_limit] <- NA
mean_nirv_cleaned <- cellStats(r_cleaned, "mean", na.rm = TRUE)

df3 <- data.frame(
  window = 1:10,
  mean = mean_nirv_cleaned,
  scaled_mean= normalize(mean_nirv_cleaned),
  vi ="NIRv"
)

df4 <- data.frame(
  window = 1:10,
  mean = NA,
  scaled_mean= (normalize(mean_nirv_cleaned)+normalize(mean_gpp_cleaned)+normalize(mean_ndvi_cleaned))/3,
  vi ="Average"
)
df <- rbind(df1,df2,df3,df4)
write.csv(df,"./results/globalmean_ridge_vpdsens_trend.csv")

p <- ggplot(df, aes(x = window, y = scaled_mean,group=vi)) +
  geom_line(aes(color = vi), size = 0.8) + 
  geom_line(data = subset(df, vi == "Average"), aes(color = vi), size = 1.2) +
  geom_smooth(data = subset(df, vi == "Average"), 
              method = "lm", size = 0.9, 
              linetype = "dashed", alpha = 0.1,fill="gray",color="black") +
  scale_x_continuous(breaks = c(1, 4, 7, 10),  
                     labels = c("2000-2009", "2003-2012", "2006-2015", "2009-2018")) +  
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75,1))+
  scale_color_manual(values=c("black","#CD1818","#F97300","#435585"))+
  labs(y = "Normalized VPD sensitivity",x="Window") +
  theme_few()+
  theme(legend.position =c(0.8,0.8))



dl <- df4[,c(1,3)]
lm <- lm(scaled_mean~window,df4)
summary(lm)$r.squared
mk <- trend::mk.test(df4$scaled_mean)
ggsave("./figs/globalmean_ridge_vpdsens_trend.pdf",p, width = 8,height = 2, dpi=300)



