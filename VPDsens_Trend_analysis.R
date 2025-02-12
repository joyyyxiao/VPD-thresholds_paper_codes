library(raster)
library(trend)
library(dplyr)


mask <- shapefile("/home/j_xiao/data/Terrestrial Ecoregions of the World_WWF/mask_shpfile.shp")
sens_window <- stack(paste0("./results/FLUXCOM_GPP/FLUXCOM_GPP_vpd_sensmap_alpha0_ERA5_monthly_window",seq(1,276,by=12)[1:10],".tif"))
mat_sensw <- as_array(sens_window)
colnames(mat_sensw) <- 2000:2009




#Sen+MK
fun_sen <- function(x){
  if(length(na.omit(x)) <10) return(c(NA, NA, NA))   
  MK_estimate <- trend::sens.slope(ts(na.omit(x), start = 2000, end = 2009, frequency = 1), conf.level = 0.95) 
  slope <- MK_estimate$estimate
  MK_test <- MK_estimate$p.value
  Zs <- MK_estimate$statistic
  results <- c(Zs, slope, MK_test)
  names(results) <-  c("Z", "slope", "p-value")
  return(results)
}

##Loop

library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(20)
registerDoSNOW(cl)

# Setup the progress bar
print(paste("Starting parellel processing at:", Sys.time()))

pb <- txtProgressBar(max =dim(mat_sensw)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

ptime  <- system.time(
  senstrend.df <- foreach(
    line=rownames(mat_sensw), .combine = rbind, .options.snow = opts,
    .packages = c('trend')) %dopar% {
      
      fun_sen(mat_sensw[line,])
      
    })
print(paste("\n Parellel processing complete at:", Sys.time()))

stopCluster(cl)

rownames(senstrend.df) <- lonlat
senstrend.df <- as.data.frame(senstrend.df)
r_slope <- as_raster(senstrend.df$slope)
r_slope <- mask(r_slope,mask, inverse = TRUE)
writeRaster(r_slope,'./results/FLUXCOM_GPP/FLUXCOM_GPP_monthly_window10_trends_slope.tif')

# senstrend.df$Level <- NA

##increase.lev1=1
##increase.lev2=2
##increase.lev3=3
##decrease.lev1=-1
##decrease.lev2=-2
##decrease.lev3=-3


# #
# senstrend.df[rownames(filter(senstrend.df,slope>0,abs(Z)>=1.65&abs(Z)<1.96)),4] <- 1
# senstrend.df[rownames(filter(senstrend.df,slope>0,abs(Z)>=1.96&abs(Z)<2.58)),4] <- 2
# senstrend.df[rownames(filter(senstrend.df,slope>0,abs(Z)>=2.58)),4] <- 3
# senstrend.df[rownames(filter(senstrend.df,slope==0)),4] <- 0
# senstrend.df[rownames(filter(senstrend.df,slope<0,abs(Z)>=1.65&abs(Z)<1.96)),4] <- -1
# senstrend.df[rownames(filter(senstrend.df,slope<0,abs(Z)>=1.96&abs(Z)<2.58)),4] <- -2
# senstrend.df[rownames(filter(senstrend.df,slope<0,abs(Z)>=2.58)),4] <- -3
# 
# 
# 
# 
# level <- data.frame(Level=c("in_1","in_2","in_3","de_1","de_2","de_3"),
#                     count= c(nrow(filter(senstrend.df,Level==1)),nrow(filter(senstrend.df,Level==2)),nrow(filter(senstrend.df,Level==3)),
#                              nrow(filter(senstrend.df,Level== -1)),nrow(filter(senstrend.df,Level== -2)),nrow(filter(senstrend.df,Level== -3))))
# library(ggplot2)
# ggplot(level,aes(x=Level,y=count,fill=Level))+
#   geom_bar(stat = "identity",width=0.8)+
#   scale_fill_manual(values=c("#fec980","#f17c4a","#d7191c","#c7e9ad","#80bfac","#2b83ba"))
  


##slope&trend
r23 <- raster("./results/MODIS_NIRv/NIRv_vpd_monthly_sensmap_alpha0_ERA5.tif")
r23 <- mask(r23,mask,inverse=T)
writeRaster(r23,"./results/MODIS_NIRv/NIRv_vpd_monthly_sensmap_alpha0_ERA5.tif",overwrite=T)
df <- data.frame(sens23years = as_array(vpd_sens_map)[,1],
                 trendslope = as_array(r_slope)[,1],
                 ltsens.vs.trend = NA)

df[rownames(filter(df,trendslope>0,sens23years<0)),3] <- -1 #"-+"
df[rownames(filter(df,trendslope<0,sens23years<0)),3] <- -2 #"--"
df[rownames(filter(df,trendslope>0,sens23years>0)),3] <- 2 #"++"
df[rownames(filter(df,trendslope<0,sens23years>0)),3] <- 1 #"+-"

r <- as_raster(df[,3])
plot(r)
writeRaster(r,'./results/FLUXCOM_GPP/FLUXCOM_GPP_monthly_window10_ltsens.vs.trend.tif')

class <- data.frame(ltsens.vs.trend=c("--","-+","++","+-"),
                    count= c(nrow(filter(df,ltsens.vs.trend== -2)),nrow(filter(df,ltsens.vs.trend== -1)),nrow(filter(df,ltsens.vs.trend== 2)),
                             nrow(filter(df,ltsens.vs.trend== 1))))
class$fraction <- class$count / sum(class$count)

# Compute the cumulative percentages (top of each rectangle)
class$ymax <- cumsum(class$fraction)

# Compute the bottom of each rectangle
class$ymin <- c(0, head(class$ymax, n=-1))

# Compute label position
class$labelPosition <- (class$ymax + class$ymin) / 2

# Compute a good label
class$label <- paste0(class$ltsens.vs.trend, "\n Fraction: ", scales::percent(class$fraction))



library(ggplot2)
# Make the donut plot
ggplot(class, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ltsens.vs.trend)) +
 geom_rect() +
 # geom_label( x=3.5, aes(y=labelPosition, label=label), size=5) +
  scale_fill_manual(values=c("#2F6B9D","#6A9ED4","#8A7356","#6B3C00")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

ggplot(class, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ltsens.vs.trend)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label, color=ltsens.vs.trend), size=4) + # x here controls label position (inner / outer)
  scale_fill_manual(values=c("#2F6B9D","#6A9ED4","#8A7356","#6B3C00")) +
  scale_color_manual(values=c("#2F6B9D","#6A9ED4","#8A7356","#6B3C00")) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")

#bar plot
ggplot(class,aes(x=ltsens.vs.trend,y=count,fill=ltsens.vs.trend))+
  geom_bar(stat = "identity",width=0.8)+
  scale_fill_manual(values=c("#d7191c","#fec980","#80bfac","#2b83ba"))


write.csv(senstrend.df,'./results/NIRv_annual_window10_ltsens.vs.trend.csv')