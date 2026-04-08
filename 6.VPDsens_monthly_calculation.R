library(glmnet)
library(raster)
library(dplyr)
load("suppl_data/data_monthly.RData")
meco2 <- read.csv("results_mr1/MODIS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- read.csv("results_mr1/GLASS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- read.csv("results_mr1/FLUXCOM_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- as.matrix(meco2)
meco2 <- meco2[,1:276]
mpar <- mpar [,1:228]
msm <- msm [,1:228]
mtmp <- mtmp [,1:228]
mvpd <- mvpd [,1:228]
mpre <- mpre [,1:228]
##############
sensi_mon_alpha0_win10 <- function(eco2,tmp,pre,sm,vpd,par,line){

  if(all(is.na(eco2))){
    
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)  
    rownames(temp) <- line
    return(temp)
  }
  
  
      mat <- matrix(data = c(eco2,tmp,pre,sm,vpd,par),nrow=length(eco2),ncol=6,
                  dimnames = list(1:length(eco2),c("y","tmp","pre","sm","vpd","par")))
      mat <- na.omit(mat)
      
  if(nrow(mat)<100){
        temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)
        rownames(temp) <- line
        return(temp)
  }
      x <- mat[,-1]
      y <- mat[,1]
      
  if(sd(y)==0){
        temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)
        rownames(temp) <- line
        return(temp)
  }
  if (any(apply(x, 2, sd) == 0)) {
        temp <- data.frame(lambda.min=NA, mse=NA, intercept=NA, 
                           coef.tmp=NA, coef.pre=NA, coef.sm=NA, 
                           coef.vpd=NA, coef.par=NA)
        rownames(temp) <- line
        return(temp)
      }
      alpha0.fit <- try(cv.glmnet(x, y, type.measure = "mse", alpha=0, family="gaussian"), silent = TRUE)
      if (inherits(alpha0.fit, "try-error")) {
        temp <- data.frame(lambda.min=NA, mse=NA, intercept=NA, 
                           coef.tmp=NA, coef.pre=NA, coef.sm=NA, 
                           coef.vpd=NA, coef.par=NA)
        rownames(temp) <- line
        return(temp)
      }
    
      lambda.min <- alpha0.fit$lambda.min
      mse <- alpha0.fit$cvm[[alpha0.fit$index[1]]]
      coef <- coef(alpha0.fit, s = "lambda.min")
      temp <- data.frame(lambda.min=lambda.min,mse=mse,intercept=coef["(Intercept)",],coef.tmp=coef["tmp",],coef.pre=coef["pre",],coef.sm=coef["sm",],coef.vpd=coef["vpd",],coef.par=coef["par",])
      rownames(temp) <- line
      return(temp)
}


################
library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(4)
registerDoSNOW(cl)

# Setup the progress bar


pb <- txtProgressBar(max =dim(meco2)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

print(paste("Starting parellel processing at:", Sys.time()))
ptime  <- system.time(
  sensi.df <- foreach(
    line=rownames(meco2), .combine = rbind, .options.snow = opts,
    .packages = c('glmnet')) %dopar% {
      
      sensi_mon_alpha0_win10(eco2 = zscore(meco2[line,]),
                             tmp = zscore(mtmp[line,]),
                             pre = zscore(mpre[line,]),
                             sm = zscore(msm[line,]),
                             vpd = zscore(mvpd[line,]),
                             par = zscore(mpar[line,]),line)

    })
print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)


###make raster###
mask <- raster("suppl_data/SYNMAP_biomass_4type_2000.tif")
vpd_sens_map <- as_raster(sensi.df$coef.vpd)
plot(vpd_sens_map,col=rainbow(20))
vpd_sens_map <- mask(vpd_sens_map,mask)
writeRaster(vpd_sens_map,"results_mr1/MODIS_GPP_sens_alpha0_00-22.tif")
writeRaster(vpd_sens_map,"results_mr1/GLASS_GPP_sens_alpha0_00-22.tif")
writeRaster(vpd_sens_map,"results_mr1/FLUXCOM_GPP_sens_alpha0_00-22.tif")


###############Trend analysis#############
##Moving window
#window == 10 years ==120 months

library(foreach)
library(doSNOW)
library(parallel)
library(glmnet)
library(trend)
library(raster)

cl <- makeCluster(4)
registerDoSNOW(cl)


n_lines <- nrow(meco2)
pb <- txtProgressBar(max = n_lines, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

window_indices <- seq(1, 228, by = 12)[1:10]
win_len <- 120 # 10 years

print(paste("Starting parellel processing at:", Sys.time()))

final_results <- foreach(
  line_idx = 1:n_lines, 
  .combine = rbind, 
  .options.snow = opts,
  .packages = c('glmnet', 'trend')
) %dopar% {
  
  line_name <- rownames(meco2)[line_idx]
  
  sens_series <- numeric(length(window_indices))
  
  for(k in seq_along(window_indices)) {
    i <- window_indices[k]
    
    eco2_sub <- meco2[line_name, i:(i + win_len - 1)]
    res <- sensi_mon_alpha0_win10(
      eco2 = zscore(eco2_sub),
      tmp  = zscore(mtmp[line_name, i:(i + win_len - 1)]),
      pre  = zscore(mpre[line_name, i:(i + win_len - 1)]),
      sm   = zscore(msm[line_name, i:(i + win_len - 1)]),
      vpd  = zscore(mvpd[line_name, i:(i + win_len - 1)]),
      par  = zscore(mpar[line_name, i:(i + win_len - 1)]),
      line = line_name
    )

    sens_series[k] <- res$coef.vpd
  }
  
  # --- Sen's Slope ---
  if(length(na.omit(sens_series)) < 10) {
    out <- c(Z = NA, slope = NA, p_value = NA)
  } else {
    tryCatch({
      MK_estimate <- trend::sens.slope(ts(na.omit(sens_series), start = 2000, frequency = 1))
      out <- c(
        Z = MK_estimate$statistic,
        slope = MK_estimate$estimate,
        p_value = MK_estimate$p.value
      )
    }, error = function(e) {
      out <- c(Z = NA, slope = NA, p_value = NA)
    })
  }
  
  return(out)
}
registerDoSEQ()
stopCluster(cl)

write.csv(final_results,"results_mr1/MODIS_GPP_sens_trend.csv")
write.csv(final_results,"results_mr1/FLUXCOM_GPP_sens_trend.csv")
write.csv(final_results,"results_mr1/GLASS_GPP_sens_trend.csv")


##slope&trend
r_sens <- raster("results_mr1/MODIS_GPP_sens_alpha0_00-22.tif")
m_trend <- read.csv("results_mr1/MODIS_GPP_sens_trend.csv",row.names = 1)
r_sens <- raster("results_mr1/FLUXCOM_GPP_sens_alpha0_00-22.tif")
m_trend <- read.csv("results_mr1/FlUXCOM_GPP_sens_trend.csv",row.names = 1)
r_sens <- raster("results_mr1/GLASS_GPP_sens_alpha0_00-22.tif")
m_trend <- read.csv("results_mr1/GLASS_GPP_sens_trend.csv",row.names = 1)
df <- data.frame(sens = as.numeric(as.matrix(r_sens)),
                 trendslope = as.numeric(m_trend$slope),
                 ltsens.vs.trend = NA)

df$ltsens.vs.trend[which(df$trendslope > 0 & df$sens < 0)] <- -1 #"Decreased +"
df$ltsens.vs.trend[which(df$trendslope < 0 & df$sens < 0)] <- -2 #"Increased -"
df$ltsens.vs.trend[which(df$trendslope > 0 & df$sens > 0)] <-  2 #"Increased +"
df$ltsens.vs.trend[which(df$trendslope < 0 & df$sens > 0)] <-  1 #"Decreased -"

r <- as_raster(df[,3])
plot(r)

writeRaster(r,'results_mr1/MODIS_GPP_window10_ltsens.vs.trend.tif')
writeRaster(r,'results_mr1/FLUXCOM_GPP_window10_ltsens.vs.trend.tif')
writeRaster(r,'results_mr1/GLASS_GPP_window10_ltsens.vs.trend.tif')

#########
library(ggplot2)
library(ggthemes)

g_1 <- raster('results_mr1/GLASS_GPP_window10_ltsens.vs.trend.tif')
f_1 <- raster('results_mr1/FLUXCOM_GPP_window10_ltsens.vs.trend.tif')
m_1 <- raster('results_mr1/MODIS_GPP_window10_ltsens.vs.trend.tif')
df <- as.data.frame(g_1, cells = FALSE, na.rm = TRUE)
colnames(df) <- "value"
stat_df <- df %>%
  group_by(value) %>%
  summarise(n = n()) %>%
  arrange(value)

stat_df

#Fig.5b, d, f
p <- ggplot(stat_df, aes(x = factor(value), y = n, fill =  factor(value))) +
  geom_col(width = 0.8) +
  scale_fill_manual(values=c("#2F6B9D","#6A9ED4","#8A7356","#6B3C00")) +
  labs(x = NULL, y = "Pixel count") +
#  scale_x_discrete(labels = c(
#    "-2" = "Increased -",
#    "-1" = "Decreased +",
#    "1"  = "Decreased -",
#    "2"  = "Increased +"
#  )) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank())
  
p
ggsave("plots/glass_sens_trend_pixelcount.pdf",p,width = 3.5,height = 2.2, dpi=300)
ggsave("plots/fluxcom_sens_trend_pixelcount.pdf",p,width = 3.5,height = 2.2, dpi=300)
ggsave("plots/modis_sens_trend_pixelcount.pdf",p,width = 3.5,height = 2.2, dpi=300)

