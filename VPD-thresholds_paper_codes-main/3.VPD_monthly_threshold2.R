library(chngpt)
library(gbm)
library(raster)
library(ggplot2)
library(ggthemes)


#load data
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
#######
lonlat

makeDF <- function(pd = pd, threshold = th, R2_train=R2_train, R2_test=R2_test){
  df <- data.frame(pd.vpd1=pd[1,"vpd"],
                   pd.vpd2=pd[2,"vpd"],
                   pd.vpd3=pd[3,"vpd"], 
                   pd.vpd4=pd[4,"vpd"], 
                   pd.vpd5=pd[5,"vpd"], 
                   pd.vpd6=pd[6,"vpd"], 
                   pd.vpd7=pd[7,"vpd"], 
                   pd.vpd8=pd[8,"vpd"], 
                   pd.vpd9=pd[9,"vpd"], 
                   pd.vpd10=pd[10,"vpd"], 
                   pd.vpd11=pd[11,"vpd"], 
                   pd.vpd12=pd[12,"vpd"], 
                   pd.vpd13=pd[13,"vpd"], 
                   pd.vpd14=pd[14,"vpd"], 
                   pd.vpd15=pd[15,"vpd"], 
                   pd.vpd16=pd[16,"vpd"], 
                   pd.vpd17=pd[17,"vpd"], 
                   pd.vpd18=pd[18,"vpd"], 
                   pd.vpd19=pd[19,"vpd"], 
                   pd.vpd20=pd[20,"vpd"], 
                   pd.mean_response1=pd[1,"mean_response"], 
                   pd.mean_response2=pd[2,"mean_response"], 
                   pd.mean_response3=pd[3,"mean_response"], 
                   pd.mean_response4=pd[4,"mean_response"], 
                   pd.mean_response5=pd[5,"mean_response"], 
                   pd.mean_response6=pd[6,"mean_response"], 
                   pd.mean_response7=pd[7,"mean_response"], 
                   pd.mean_response8=pd[8,"mean_response"], 
                   pd.mean_response9=pd[9,"mean_response"], 
                   pd.mean_response10=pd[10,"mean_response"], 
                   pd.mean_response11=pd[11,"mean_response"], 
                   pd.mean_response12=pd[12,"mean_response"], 
                   pd.mean_response13=pd[13,"mean_response"], 
                   pd.mean_response14=pd[14,"mean_response"], 
                   pd.mean_response15=pd[15,"mean_response"], 
                   pd.mean_response16=pd[16,"mean_response"], 
                   pd.mean_response17=pd[17,"mean_response"], 
                   pd.mean_response18=pd[18,"mean_response"], 
                   pd.mean_response19=pd[19,"mean_response"], 
                   pd.mean_response20=pd[20,"mean_response"],
                   threshold = threshold, R2_train=R2_train, R2_test=R2_test
                   )
  return(df)
}
#
emptydf <- data.frame(pd.vpd1=NA,
                      pd.vpd2=NA,
                      pd.vpd3=NA, 
                      pd.vpd4=NA, 
                      pd.vpd5=NA, 
                      pd.vpd6=NA, 
                      pd.vpd7=NA, 
                      pd.vpd8=NA, 
                      pd.vpd9=NA, 
                      pd.vpd10=NA, 
                      pd.vpd11=NA, 
                      pd.vpd12=NA, 
                      pd.vpd13=NA, 
                      pd.vpd14=NA, 
                      pd.vpd15=NA, 
                      pd.vpd16=NA, 
                      pd.vpd17=NA, 
                      pd.vpd18=NA, 
                      pd.vpd19=NA, 
                      pd.vpd20=NA, 
                      pd.mean_response1=NA, 
                      pd.mean_response2=NA, 
                      pd.mean_response3=NA, 
                      pd.mean_response4=NA, 
                      pd.mean_response5=NA, 
                      pd.mean_response6=NA, 
                      pd.mean_response7=NA, 
                      pd.mean_response8=NA, 
                      pd.mean_response9=NA, 
                      pd.mean_response10=NA, 
                      pd.mean_response11=NA, 
                      pd.mean_response12=NA, 
                      pd.mean_response13=NA, 
                      pd.mean_response14=NA, 
                      pd.mean_response15=NA, 
                      pd.mean_response16=NA, 
                      pd.mean_response17=NA, 
                      pd.mean_response18=NA, 
                      pd.mean_response19=NA, 
                      pd.mean_response20=NA,
                      threshold =NA, R2_train=NA, R2_test=NA
)

thres <- function(Loc){
  
  df <- data.frame(response=zscore(meco2[Loc,]),vpd=mvpd[Loc,],sm=msm[Loc,],pre=mpre[Loc,],tmp=mtmp[Loc,],par=mpar[Loc,])
  df <- na.omit(df)
  #check data
  if(nrow(df) < 200 |any(apply(df,2,sd)==0))
    return(emptydf)

  
  #
  set.seed(1234 + which(line==Loc))
  
  # train/test split
  idx <- sample(seq_len(nrow(df)), size = 0.8*nrow(df))
  train <- df[idx,]
  test  <- df[-idx,]
  
  # GBM model
  model <- gbm(
    formula = response ~ vpd + sm + pre + tmp + par,
    data = train,
    distribution = "gaussian",
    n.trees = 200,
    interaction.depth = 5,
    shrinkage = 0.05,
    bag.fraction = 0.5
  )
  # partial dependence
  vpd_grid <- seq(min(test$vpd), max(test$vpd), length = 20)
  
  pd <- data.frame(
    vpd = vpd_grid,
    mean_response = sapply(vpd_grid, function(x){
      tmp <- test
      tmp$vpd <- x
      mean(predict(model, tmp, n.trees = 200))
    })
  )
  
  fitst <- chngptm(formula.1=mean_response~1, formula.2=~vpd,pd,type="stegmented", family="gaussian")
  th <- fitst$coefficients[[5]]
  
  # predictions
  pred_train <- predict(model, train, n.trees = 200)
  pred_test  <- predict(model, test, n.trees = 200)
  
  # R2
  R2_train <- cor(train$response, pred_train)^2
  R2_test  <- cor(test$response, pred_test)^2
  
  mydata <- makeDF(pd = pd, threshold = th, R2_train=R2_train, R2_test=R2_test)
  return(mydata)
}

#loop
rr_gbm <- emptydf
for (i in 1:259200){
  Loc <- lonlat[i]
  rr <- thres(Loc)
  rownames(rr) <- Loc
  rr_gbm <- rbind(rr_gbm,rr)
  print(Loc)
}

#save data
write.csv(rr_gbm,"results_mr1/FLUXCOM_GPP_threshold2_h2o_gbm.csv")
write.csv(rr_gbm,"results_mr1/GLASS_GPP_threshold2_h2o_gbm.csv")
write.csv(rr_gbm,"results_mr1/MODIS_GPP_threshold2_h2o_gbm.csv")

####filter thresholds
fluxcom_th <- read.csv("results_mr1/FLUXCOM_GPP_threshold2_h2o_gbm.csv",row.names=1)[-1,]
glass_th <- read.csv("results_mr1/GLASS_GPP_threshold2_h2o_gbm.csv",row.names=1)[-1,]
modis_th <- read.csv("results_mr1/MODIS_GPP_threshold2_h2o_gbm.csv",row.names=1)[-1,]

vpd_cols <- paste0("pd.vpd", 1:20)
resp_cols <- paste0("pd.mean_response", 1:20)

get_slope <- function(x, y){
  
  df <- data.frame(x, y)
  df <- na.omit(df)
  
  if(nrow(df) < 5) return(c(NA, NA))
  
  m <- lm(y ~ x, data = df)
  
  slope <- coef(m)[2]
  pval <- summary(m)$coefficients[2,4]
  
  return(c(slope, pval))
}


res <- t(apply(modis_th, 1, function(row){
  
  x <- as.numeric(row[vpd_cols])
  y <- as.numeric(row[resp_cols])
  
  get_slope(x, y)
  
}))

res <- as.data.frame(res)
names(res) <- c("slope","pvalue")

sig_neg <- !is.na(res$slope) & !is.na(res$pvalue) & res$slope < 0 & res$pvalue < 0.05
other   <- !is.na(res$slope) & !is.na(res$pvalue) & !sig_neg

vals <- fluxcom_th$threshold
vals[other] <- 99
vals[is.na(res$slope) | is.na(res$pvalue)] <- NA
fluxcom_masked <- as_raster(vals)
glass_masked <- as_raster(vals)
modis_masked <- as_raster(vals)

#biome mask
biome <- raster("suppl_data/SYNMAP_biomass_4type_2000.tif")
f_th <- mask(fluxcom_masked,biome)
g_th <- mask(glass_masked,biome)
m_th <- mask(modis_masked,biome)

writeRaster(f_th,"results_mr1/FLUXCOM_threshold2_masked.tif",overwrite=TRUE)
writeRaster(g_th,"results_mr1/GLASS_threshold2_masked.tif",overwrite=TRUE)
writeRaster(m_th,"results_mr1/MODIS_threshold2_masked.tif",overwrite=TRUE)


###########relationship to multiyear mean VPD
mean_vpd <- read.csv("suppl_data/data_multiyearmean.csv")$VPD

df <- data.frame(threshold = as.numeric(as.matrix(m_th)),
                 meanvpd = as.numeric(as.matrix(mean_vpd)))

df <- na.omit(df)

fit <- lm(threshold ~ meanvpd, data = df)

slope <- coef(fit)[2]
R2 <- summary(fit)$r.squared
pval  <- summary(fit)$coefficients[2,4]
p_label <- ifelse(pval < 0.001,
                  "p < 0.001",
                  paste0("p = ", round(pval,3)))
lab <- paste0(
  "Slope = ", round(slope,3),
  "\nR² = ", round(R2,3),
  "\n", p_label
)

#Fig. 3b, d, f 
p3 <- ggplot(data=df,
            aes(x=meanvpd,
                y=threshold))+
  stat_density2d(aes(fill=after_stat(density)),n=500,geom = "raster",contour=FALSE)+
  scale_fill_gradientn(colours = c("white","grey","#36126C","#81317D","#CD546A","#F1A376","#FCFDC6"),
                       name= "Density")+
  geom_abline(slope=1,linetype="dashed",color="red")+
  geom_smooth(method="lm",linewidth=.4,color="black")+
  theme_few()+
#  scale_x_continuous(limits=c(0,2))+
  scale_y_continuous(limits=c(0,4.5))+
  xlab("Multiyear mean VPD (kPa)")+
  ylab("VPD threshold (kPa)")+
  annotate("text",
           x = Inf,
           y = -Inf,
           label = lab,
           hjust = 1.1,
           vjust = -0.5,
           size = 4)


ggsave("fluxcom_1y1ratio_Plot.pdf",p1, width = 5,height = 4, dpi=300)
ggsave("glass_1y1ratio_Plot.pdf",p2, width = 5,height = 4, dpi=300)
ggsave("modis_1y1ratio_Plot.pdf",p3, width = 5,height = 4, dpi=300)