###########Data preparing#########
library(raster)
library(dplyr)

setwd("/home/j_xiao/data/NCdemo")
co2ef <- matrix(read.csv("./results/FLUXCOM_GPP/FLUXCOM_GPP_total_CO2_effect.csv")[,3],259200,1)
co2ef <- matrix(read.csv("./results/MODIS_NIRv/NIRv_total_CO2_effect.csv")[,3],259200,1)
co2ef <- matrix(read.csv("./results/GIMMS_NDVI/GIMMS_NDVI_total_CO2_effect.csv")[,3],259200,1)
rownames(co2ef) <- read.csv("line.csv")[,-1]
co2ef[which(co2ef<=0)] <- NA

# 
# ##TMP## 
# rtmp <- stack("/home/j_xiao/data/CRUclimate/cru_ts4.06.1981.1990.tmp.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.1991.2000.tmp.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2001.2010.tmp.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2011.2020.tmp.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.07.2021.2022.tmp.dat.nc")
# rtmp <- mean(rtmp)
# tmp <- as.numeric(as.matrix(rtmp))
# write.csv(tmp,"/home/j_xiao/data/NCdemo/data/FLUXCOM_GPP/MeanTMP1981-2022.csv",na="")
# write.csv(tmp,"/home/j_xiao/data/NCdemo/data/GIMMS_NDVI/MeanTMP1981-2022.csv",na="")
# rm(rtmp)
# 
# 
# ##PRE##
# 
# rpre <- stack("/home/j_xiao/data/CRUclimate/cru_ts4.06.1981.1990.pre.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.1991.2000.pre.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2001.2010.pre.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2011.2020.pre.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.07.2021.2022.pre.dat.nc")
# rmap <- sum(rpre)/42
# map <- as.numeric(as.matrix(rmap))
# write.csv(map,"/home/j_xiao/data/NCdemo/data/FLUXCOM_GPP/MeanMAP1981-2022.csv",na="")
# write.csv(map,"/home/j_xiao/data/NCdemo/data/GIMMS_NDVI/MeanMAP1981-2022.csv",na="")
# rm(rmap,rpre)
# 
# ##PET##
# rpet <- stack("/home/j_xiao/data/CRUclimate/cru_ts4.06.1981.1990.pet.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.1991.2000.pet.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2001.2010.pet.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.06.2011.2020.pet.dat.nc",
#               "/home/j_xiao/data/CRUclimate/cru_ts4.07.2021.2022.pet.dat.nc")
# rpet <-mean(rpet) 
# pet <- as.numeric(as.matrix(rpet))
# write.csv(pet,"/home/j_xiao/data/NCdemo/data/FLUXCOM_GPP/MeanPET1981-2022.csv",na="")
# write.csv(pet,"/home/j_xiao/data/NCdemo/data/GIMMS_NDVI/MeanPET1981-2022.csv",na="")
# rm(rpet)

map <- read.csv("./data/NIRv/MeanMAP00-22.csv")[,2]
tmp <- read.csv("./data/NIRv/MeanTMP00-22.csv")[,2]
pet <- read.csv("./data/NIRv/MeanPET00-22.csv")[,2]
tp <- read.csv("./datasup/potential predictors/demo_dataframe_total_P.csv")[,2]
cn <- read.csv("./datasup/potential predictors/demo_dataframe_cn_topsoil_mean.csv")[,2]
biome <-  read.csv("./datasup/potential predictors/demo_dataframe_biome_type.csv")[,2]


biome[biome[]%in%c(1,10,20,30)] <- "tree"
biome[biome[]%in%c(37,38,39,40)] <- "shrubs"
biome[biome[]%in%c(41,42,43)] <- "grasses"
biome[biome[]%in%c(44)] <- "crops"

myc <- read.csv("./datasup/potential predictors/demo_dataframe_MycType.csv")[,2]
myc[which(myc==1)] <- "AM"
myc[which(myc==2)] <- "ECM"
myc[which(myc==4)] <- "NM"

mi <- map/(pet*365) # pet_units=mm/day,pre_units=mm/month,map_units=mm/year


##VPD##
rvpd <- raster("/home/j_xiao/data/TerraClimate/MeanVPD_8122.tif")
vpd <- as.numeric(as.matrix(rvpd))
vpd <- vpd*0.01

##SM##
rsm <- raster("/home/j_xiao/data/TerraClimate/MeanSM_8122.tif")
sm <- as.numeric(as.matrix(rsm))
sm <- sm*0.1

##AET##
raet <- raster("/home/j_xiao/data/TerraClimate/MeanAET_8122.tif")
aet <- as.numeric(as.matrix(raet))
aet <- aet*0.1

rm(rvpd,rsm,raet)

##PAR## 2000/03-2022/12
dirpar <- stack("CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_direct_clr_mon")
diffpar <- stack("CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_diff_clr_mon")
totalpar <- dirpar+diffpar
meantotalpar <- mean(totalpar)
meantotalpar <- disaggregate(meantotalpar,fact=2,)
par <-  as.numeric(as.matrix(meantotalpar))

rm(dirpar,diffpar,totalpar,meantotalpar,raet,rsm,rvpd)

####################
data <- data.frame(eco2=co2ef,MAP=map,TMP=tmp,PET=pet,MI=mi,AET=aet,VPD=vpd,SM=sm,P=tp,CNr=cn,Biome=biome,Myc=myc,PAR=par)
data <- na.omit(data)
data$Biome <- factor(data$Biome)
data$Myc <- factor(data$Myc)


#####H2O#####
library(h2o)
h2o.init(-1)
h2o.removeAll()
#h2o.no_progress()

df <- as.h2o(data)

predictors <- c("VPD", "SM", "MAP", "TMP", "PET","MI","AET","P","CNr","Biome","Myc","PAR")
response <- "eco2"

# split into train and validation sets
df_splits <- h2o.splitFrame(data =  df, ratios = 0.8, seed = 1234)
train <- df_splits[[1]]
test <- df_splits[[2]]

aml_full <- h2o.automl(x = predictors, y = response,
                     training_frame = df,
                     max_models = 10,
                     nfolds = 5,
                     exclude_algos = c("GLM", "DeepLearning"),
                     seed = 1234) 

aml_ndvi_80 <- h2o.automl(x = predictors, y = response,
                  training_frame = train,
                  max_models = 10,
                  nfolds = 5,
                  exclude_algos = c("GLM", "DeepLearning"),
                  seed = 1234) 

# View the AutoML Leaderboard
lb_gpp80 <- aml_gpp_80@leaderboard
lb_nirv80 <- aml_nirv_80@leaderboard
lb_ndvi80 <- aml_ndvi_80@leaderboard


print(lb_gpp80, n = nrow(lb_gpp80))  # Print all rows instead of default (6 rows)
print(lb_nirv80, n = nrow(lb_nirv80))  # Print all rows instead of default (6 rows)
print(lb_ndvi80, n = nrow(lb_ndvi80))  # Print all rows instead of default (6 rows)


bmodel <-  h2o.get_best_model(aml_gpp_80)

# get best model family
model_SE <- h2o.getModel("StackedEnsemble_AllModels_1_AutoML_4_20231123_143505")
model_DRF <- h2o.getModel("DRF_1_AutoML_3_20241010_163849")
model_XRT <- h2o.getModel("XRT_1_AutoML_4_20231123_143505")
model_GBM <- h2o.getModel("GBM_4_AutoML_4_20231123_143505")
model_XGBoost <- h2o.getModel("XGBoost_1_AutoML_4_20231123_143505")

h2o.r2(model_SE)
h2o.r2(model_DRF)
h2o.r2(model_XRT)
h2o.r2(model_GBM)
h2o.r2(model_XGBoost)

#PartialDependence: Partial dependency plot for vpd
pdvpd_modelSE <- h2o.partialPlot(object = model_SE,
                               data = test,
                               cols = "vpd",
                               plot_stddev = TRUE)
pdvpd_modelDRF <- h2o.partialPlot(object = model_DRF,
                                 data = test,
                                 cols = "vpd",
                                 plot_stddev = TRUE)
pdvpd_modelXRT <- h2o.partialPlot(object = model_XRT,
                                 data = test,
                                 cols = "vpd",
                                 plot_stddev = TRUE)
pdvpd_modelGBM <- h2o.partialPlot(object = model_GBM,
                                 data = test,
                                 cols = "vpd",
                                 plot_stddev = TRUE)
pdvpd_modelXGBoost <- h2o.partialPlot(object = model_XGBoost,
                                 data = test,
                                 cols = "vpd",
                                 plot_stddev = TRUE)

pdmatrix <- matrix(c(pdvpd_modelSE$mean_response,
pdvpd_modelDRF$mean_response,
pdvpd_modelXRT$mean_response,
pdvpd_modelGBM$mean_response,
pdvpd_modelXGBoost$mean_response),20,5)

model_median <- apply(pdmatrix, 1, median)

pddf <- data.frame(VPD=rep(pdvpd_modelSE$vpd,6),
                   Model=rep(c("Stacked Esemble","DRF","XRT","GBM","XGBoost","Median of models"),each=20),
                   mean_response=c(pdvpd_modelSE$mean_response,pdvpd_modelDRF$mean_response,pdvpd_modelXRT$mean_response,pdvpd_modelGBM$mean_response,pdvpd_modelXGBoost$mean_response,model_median),
                   std=c(pdvpd_modelSE$stddev_response,pdvpd_modelDRF$stddev_response,pdvpd_modelXRT$stddev_response,pdvpd_modelGBM$stddev_response,pdvpd_modelXGBoost$stddev_response,rep(NA,20)))
pddf$Model <- factor(pddf$Model, levels = c("DRF","XRT","GBM","XGBoost","Stacked Esemble","Median of models")) 
write.csv(pddf,"./results/FLUXCOM_GPP_vpd_pdplots_data.csv")

pddf[pddf$Model=="Median of models",]

ggplot(pddf,aes(x=VPD,y=mean_response,group=Model))+
  geom_ribbon(aes(ymin = mean_response - std, ymax = mean_response + std),fill="gray70",alpha=0.1)+
  geom_line(aes(color=Model,linewidth=Model))+
  geom_vline(xintercept = 0.27116959,linetype="dashed")+
  scale_color_manual(values=c("#B0A8B9","#4B4453","#FF8066","#C34A36","red","black"))+
  scale_linewidth_manual(values = c(1,1,1,1,1,1.8))+
  xlab("Multiyear mean VPD")+
  ylab("Mean response")+
  theme_few(base_size=15)+
  theme(#text=element_text(family="arial"),
        axis.title =element_text(face="bold"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black"),
        legend.background = element_rect(fill = NA),
        legend.position = c(.8, .8),
        legend.text = element_text(size = 9),
        legend.title = element_blank())



pred <- h2o.predict(object = rfmodel, newdata = test)
valid <- data.frame(prediction=as.vector(pred),observation=as.vector(test[,1]))

library(ggplot2)
ggplot(valid,aes(x=prediction,y=observation))+
  geom_point(alpha=.25)+
  geom_abline(intercept = 0,slope = 1,linetype="dashed",col="red")+
  theme_classic()

pdvpd <- h2o.partialPlot(object = rfmodel,
                      data = test,
                      cols = "vpd",
                      plot_stddev = FALSE)

va_plot <- h2o.varimp_heatmap(aml_gpp_80)

pd_plot <- h2o.pd_multi_plot(aml_gpp_80, test, "vpd")
pd_plot <- h2o.pd_multi_plot(aml_ndvi_80, test, "pet")

lbb <- as.data.frame(lb_gpp80)
write.csv(lbb,"/home/j_xiao/data/NCdemo/results/FLUXCOM_GPP/GPP_leaderboard.csv")

####RF####
library(randomForest)
library(pacman)
library(caret)
library(skimr)
library(pdp)
#data02 <- data02[sample(nrow(data02),size=100000,replace=F),]

set.seed(1234)

trains <- createDataPartition( #caret??
  y = data02$eco2,
  p = 0.8,
  list = F
)
data_train <- data02[trains,]
data_test <- data02[-trains,]

#function
rf <- randomForest(eco2~ ., data = data_train,importance = TRUE,ntree = 100)
#formula1 <- eco2 ~ map + tmp + pet + p + cnr + biome + myc
#rf1 <- randomForest(formula1, data = data_train,importance = TRUE,ntree = 100)
#formula2 <- eco2 ~ mi + tmp + p + cnr + biome + myc
#rf2 <- randomForest(formula2, data = data_train,importance = TRUE,ntree = 100)

varImpPlot(rf)
imp <- importance(rf,type=2)


#SHAP
# build the standardized coefficient magnitudes plot:
h2o.std_coef_plot(pros_glm)


# Predict the contributions using the model and test data:
contributions <- h2o.predict_contributions(aml, test)

# Plot SHAP summary plot:
h2o.shap_summary_plot(model_DRF, test)

# Plot SHAP contributions for one instance (e.g., row 5):
h2o.shap_explain_row_plot(model, prostate_test, row_index = 5)

library(ggplot2)
library(gridExtra)
colnames(imp) <- "Importance"
rownames(imp) <- c("Mean Annual Precipitation",
                   "Mean Annual Temperature",
                   "Potential Evapotranspiration",
                   "Moisture index",
                   "Actual Evapotranspiration",
                   "Vapour pressure deficit",
                   "Soil moisture",
                   "Total Phosphorus",
                   "Soil C:N ratio",
                   "Biome type",
                   "Mycorrhizal type")

imp <- data.frame(predictor=rownames(imp)[order(imp,decreasing = T)],
                  importance=imp[order(imp,decreasing = T)])





p <- ggplot(imp, aes(x=reorder(predictor,importance), y=importance)) +
  geom_bar(stat='identity',color="black",size=0.2, fill="#12639C") + #D9E6F6 #FDDB7D
  xlab("") + ylab ("Relative Importance of predictors") +
  coord_flip() + 
  xlab(NULL) + 
  #geom_hline(yintercept=cutoff,col="black",linetype="dashed") +
  #theme_few() + 
  #theme(legend.title = element_blank(),legend.position = c(0.73, 0.16),
  #      legend.box.background = element_rect(fill = "transparent"),
  #      legend.background = element_rect(fill = "transparent"),
  #      plot.margin = unit(rep(1,4),"lines"))

predictions <- predict(rf, newdata = data_test)
plot(predictions, data_test$eco2)


partialPlot(x = rf,
            pred.data = data_train,
            x.var = map)
plot(eco2 ~ vpd, data = data02)

ggplot(data02, aes(x=vpd, y=eco2))+
  
  geom_point() +
  geom_smooth(color="red",se=T,level = 0.95) +
  theme_bw()+
  ylab(expression(paste(CO[2],"effect on NDVI", sep="")))+
  xlab("VPD")


pdp1 <- partial(rf, pred.var = "vpd") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Vapour pressure deficit")
pdp2 <- partial(rf, pred.var = "pet") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Potential Evapotranspiration ")
pdp3 <- partial(rf, pred.var = "tmp") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Mean Annual Temperature")
pdp4 <- partial(rf, pred.var = "mi") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Moisture index")
pdp5 <- partial(rf, pred.var = "biome") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Biome type")
pdp6 <- partial(rf, pred.var = "aet") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Actual Evapotranspiration")
pdp7 <- partial(rf, pred.var = "sm") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Soil moisture")
pdp8 <- partial(rf, pred.var = "map") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Mean Annual Precipitation")
pdp9 <- partial(rf, pred.var = "cnr") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Soil C:N ratio")
pdp10 <- partial(rf, pred.var = "myc") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Mycorrhizal type")
pdp11 <- partial(rf, pred.var = "p") %>% plotPartial(smooth = F, lwd = 2, ylab = "eCO2 effects",xlab="Total Phosphorus")

pdp_vpd <- partial(rf, pred.var = "vpd")
grid.arrange(pdp1, pdp2, pdp3, pdp4,
             pdp5, pdp6, pdp7, pdp8,
             pdp9, pdp10,pdp11,
             ncol = 4, nrow = 3)




