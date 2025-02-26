###########Data preparing#########
library(raster)
library(dplyr)

setwd("/home/j_xiao/data/NCdemo")
data <- read.csv(".data/demo_ML_analysis_data.csv")
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

#SHAP
# build the standardized coefficient magnitudes plot:
h2o.std_coef_plot(pros_glm)

# Predict the contributions using the model and test data:
contributions <- h2o.predict_contributions(aml, test)

# Plot SHAP summary plot:
h2o.shap_summary_plot(model_DRF, test)

# Plot SHAP contributions for one instance (e.g., row 5):
h2o.shap_explain_row_plot(model, prostate_test, row_index = 5)




