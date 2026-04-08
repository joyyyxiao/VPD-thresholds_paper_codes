###########Data preparing#########
library(raster)
library(dplyr)
library(ggplot2)
library(ggthemes)


data <- read.csv("demo_data/demo_data_multiyearmean.csv")
clean_outliers <- function(x, probs = c(0.01, 0.99)) {
 
  x[!is.finite(x)] <- NA
  qs <- quantile(x, probs = probs, na.rm = TRUE)
    x[x < qs[1] | x > qs[2]] <- NA
  
  return(x)
}

#load eCO2 effect data
total_eco2 <- read.csv("results_mr1/FLUXCOM_gpp_total_CO2_effect.csv")
total_eco2 <- read.csv("results_mr1/GLASS_gpp_total_CO2_effect.csv")
total_eco2 <- read.csv("results_mr1/MODIS_gpp_total_CO2_effect.csv")

data$eco2_re <- clean_outliers(total_eco2$RelativeChange)
data$eco2_ab <- clean_outliers(total_eco2$AbsoluteChange)

#set factors
data$Biome <- factor(data$Biome)
data$Myc <- factor(data$Myc)
#data00 <- na.omit(data)

# Load land-use change masks at different percentage thresholds 
# (the percentage indicates the proportion of each pixel that experienced land-use change during 2000–2022)
mask_20 <- raster("demo_data/lulcc_0022_20perc.tif")
mask_20[mask_20 == 1] <- NA
mask_30 <- raster("demo_data/lulcc_0022_30perc.tif")
mask_30[mask_30 == 1] <- NA
mask_40 <- raster("demo_data/lulcc_0022_40perc.tif")
mask_40[mask_40 == 1] <- NA
mask_50 <- raster("demo_data/lulcc_0022_50perc.tif")
mask_50[mask_50 == 1] <- NA


data$lccmask <- as.numeric(as.matrix(mask_20))
data20 <- na.omit(data)
data$lccmask <- as.numeric(as.matrix(mask_30))
data30 <- na.omit(data)
data$lccmask <- as.numeric(as.matrix(mask_40))
data40 <- na.omit(data)
data$lccmask <- as.numeric(as.matrix(mask_50))
data50 <- na.omit(data)


#####H2O#####
Sys.setenv(JAVA_HOME = system("/usr/libexec/java_home -v 17", intern = TRUE))
system("java -version")
library(h2o)
h2o.init()
#h2o.no_progress()

predictors <- c("VPD", "SM", "MAP", "TMP", "PET","MI","AET","P","CNr","Biome","Myc","PAR")
response <- "eco2_re"
response <- "eco2_ab"
h2o.removeAll()
df <- as.h2o(data20)
aml <- h2o.automl(x = predictors, y = response,
                  training_frame = df,
                  max_models = 10,
                  nfolds = 5,
                  include_algos = c("GBM", "DRF","StackedEnsemble"),
                  seed = 1234) 

pd_obs <- as.data.frame(
  h2o.partialPlot(
    object = aml@leader,
    data = df,
    cols = "VPD",
    nbins = 20,
    plot = FALSE
  )
)
varimp_m <- h2o.varimp_heatmap(aml)
varimp_g <- h2o.varimp_heatmap(aml)
varimp_f <- h2o.varimp_heatmap(aml)
combined_plot <- gridExtra::grid.arrange(varimp_g, varimp_m, varimp_f, ncol = 1)
ggsave("plots/varimp_combined_ab.pdf", combined_plot, width = 5, height = 15, dpi = 300)
ggsave("plots/varimp_combined_re.pdf", combined_plot, width = 5, height = 15, dpi = 300)

###############
# View the AutoML Leaderboard
lb_df <- as.data.frame(aml@leaderboard)

# get best model family

model_SE <- h2o.getModel(lb_df$model_id[grep("^StackedEnsemble", lb_df$model_id)][1])
model_DRF <- h2o.getModel(lb_df$model_id[grep("^DRF", lb_df$model_id)][1])
model_XRT <- h2o.getModel(lb_df$model_id[grep("^XRT", lb_df$model_id)][1])
model_GBM <- h2o.getModel(lb_df$model_id[grep("^GBM", lb_df$model_id)][1])

lb_df$R2 <- NA_real_
lb_df$R2[grep("^StackedEnsemble", lb_df$model_id)[1]] <- as.numeric(h2o.r2(model_SE))
lb_df$R2[grep("^DRF", lb_df$model_id)[1]] <- as.numeric(h2o.r2(model_DRF))
lb_df$R2[grep("^XRT", lb_df$model_id)[1]] <- as.numeric(h2o.r2(model_XRT))
lb_df$R2[grep("^GBM", lb_df$model_id)[1]] <- as.numeric(h2o.r2(model_GBM))

lb_df
write.csv(lb_df,"results_mr1/reFLUX_aml_leaderboard_mask20.csv")
write.csv(lb_df,"results_mr1/abFLUX_aml_leaderboard_mask20.csv")
write.csv(lb_df,"results_mr1/reGLASS_aml_leaderboard_mask20.csv")
write.csv(lb_df,"results_mr1/abGLASS_aml_leaderboard_mask20.csv")
write.csv(lb_df,"results_mr1/reMODIS_aml_leaderboard_mask20.csv")
write.csv(lb_df,"results_mr1/abMODIS_aml_leaderboard_mask20.csv")

#PartialDependence: Partial dependency plot for vpd
pdvpd_modelSE <- h2o.partialPlot(object = model_SE,
                               data = df,
                               cols = "VPD",
                               plot_stddev = TRUE)
pdvpd_modelDRF <- h2o.partialPlot(object = model_DRF,
                                 data = df,
                                 cols = "VPD",
                                 plot_stddev = TRUE)
pdvpd_modelXRT <- h2o.partialPlot(object = model_XRT,
                                 data = df,
                                 cols = "VPD",
                                 plot_stddev = TRUE)
pdvpd_modelGBM <- h2o.partialPlot(object = model_GBM,
                                 data = df,
                                 cols = "VPD",
                                 plot_stddev = TRUE)

model_median <- apply(matrix(c(pdvpd_modelSE$mean_response,pdvpd_modelDRF$mean_response,pdvpd_modelXRT$mean_response,pdvpd_modelGBM$mean_response),20,4), 1, median)

pddf <- data.frame(VPD=rep(pdvpd_modelSE$VPD,5),
                   Model=rep(c("Stacked Esemble","DRF","XRT","GBM","Median of models"),each=20),
                   mean_response=c(pdvpd_modelSE$mean_response,pdvpd_modelDRF$mean_response,pdvpd_modelXRT$mean_response,pdvpd_modelGBM$mean_response,model_median),
                   std=c(pdvpd_modelSE$stddev_response,pdvpd_modelDRF$stddev_response,pdvpd_modelXRT$stddev_response,pdvpd_modelGBM$stddev_response,rep(NA,20)))
pddf$Model <- factor(pddf$Model, levels = c("DRF","XRT","GBM","Stacked Esemble","Median of models")) 


#Fig. 2
p <- ggplot(pddf,aes(x=VPD,y=mean_response,group=Model))+
  geom_ribbon(aes(ymin = mean_response - std, ymax = mean_response + std),fill="gray70",alpha=0.1)+
  geom_line(aes(color=Model,linewidth=Model))+
#  geom_vline(xintercept = 0.4658511,linetype="dashed")+
  scale_color_manual(values=c("#b5a5c7","#7c7484","#C34A36","#cfa095","black"))+
  scale_linewidth_manual(values = c(1,1,1,1,1,1.8))+
  xlab("Multiyear mean VPD (kPa)")+
  ylab("Absolute response")+
#  ylab("Relative response")+
  theme_few(base_size=15)+
  theme(#text=element_text(family="arial"),
        axis.title =element_text(face="bold"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black"),
        legend.background = element_rect(fill = NA),
        legend.key = element_blank(),
        legend.box.background = element_blank(),
        legend.position = c(.8, .8),
        legend.text = element_text(size = 9),
        legend.title = element_blank())

p

ggsave("plots/glass_ab_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)
ggsave("plots/glass_re_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)
ggsave("plots/flux_re_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)
ggsave("plots/flux_ab_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)
ggsave("plots/modis_re_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)
ggsave("plots/modis_ab_pdplot.pdf",p,width = 5,height = 4.5, dpi=300)



#######random null model#######
response <- "eco2_re"
response <- "eco2_ab"
B <- 10
pd_null_list <- vector("list", B)

for (b in 1:B) {
  cat("Null run", b, "\n")
  data_null <- data20
  data_null[[response]] <- sample(data_null[[response]])
  
  df_null <- as.h2o(data_null)
  
  aml_null <- h2o.automl(
    x = predictors, y = response,
    training_frame = df_null,
    max_models = 5,
    nfolds = 5,
    include_algos = c("GBM", "DRF"),
    seed = 1234 + b
  )
  
  pd_b <- as.data.frame(
    h2o.partialPlot(
      object = aml_null@leader,
      data = df_null,
      cols = "VPD",
      nbins = 20,
      plot = FALSE
    )
  )
  
  pd_b$boot <- b
  pd_null_list[[b]] <- pd_b
  h2o.removeAll()
}

pd_null <- do.call(rbind, pd_null_list)

vpd_vals <- sort(unique(pd_null$VPD))

pd_null_sum <- data.frame(
  VPD = vpd_vals,
  mean_null = sapply(vpd_vals, function(x) mean(pd_null$mean_response[pd_null$VPD == x], na.rm = TRUE)),
  low_null  = sapply(vpd_vals, function(x) quantile(pd_null$mean_response[pd_null$VPD == x], 0.025, na.rm = TRUE)),
  high_null = sapply(vpd_vals, function(x) quantile(pd_null$mean_response[pd_null$VPD == x], 0.975, na.rm = TRUE))
)

#Fig. S4
p_null_m_ab <- ggplot() +
  geom_ribbon(data = pd_null_sum,
              aes(x = VPD, ymin = low_null, ymax = high_null),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = pd_null_sum,
            aes(x = VPD, y = mean_null, color = "Random null model (n=10)"), 
            linewidth = 1) +
  geom_line(data = pd_obs,
            aes(x = VPD, y = mean_response, color = "Leader model"), 
            linewidth = 1.2) +
  scale_color_manual(name = NULL, 
                     values = c("Random null model (n=10)" = "grey40", 
                                "Leader model" = "red")) +
  labs(x = "VPD",
       y = "Absolute response") +
  theme_classic() +
  theme(legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"), 
        legend.background = element_blank()) 

p_null_m_re <- ggplot() +
  geom_ribbon(data = pd_null_sum,
              aes(x = VPD, ymin = low_null, ymax = high_null),
              fill = "grey80", alpha = 0.6) +
  geom_line(data = pd_null_sum,
            aes(x = VPD, y = mean_null, color = "Random null model (n=10)"), 
            linewidth = 1) +
  geom_line(data = pd_obs,
            aes(x = VPD, y = mean_response, color = "Leader model"), 
            linewidth = 1.2) +
  scale_color_manual(name = NULL, 
                     values = c("Random null model (n=10)" = "grey40", 
                                "Leader model" = "red")) +
  labs(x = "VPD",
       y = "Relative response") +
  theme_classic() +
  theme(legend.position = c(0.05, 0.95), 
        legend.justification = c("left", "top"), 
        legend.background = element_blank()) 

p <- gridExtra::grid.arrange(
  p_null_g_ab, p_null_m_ab, p_null_f_ab,
  p_null_g_re, p_null_m_re, p_null_f_re,
  nrow = 2, 
  ncol = 3
)
ggsave("plots/pdd_null_models.pdf", p, width = 9, height = 6, dpi = 300)
