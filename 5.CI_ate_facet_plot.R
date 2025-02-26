library(raster)
library(ggplot2)
library(ggthemes)
library(dplyr)
setwd("/home/j_xiao/data/NCdemo/")

ndvi_cf_ate <-  raster("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_ate_all_masked.tif")
gpp_cf_ate <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_all_masked.tif")
nirv_cf_ate <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_cf_ate_all_masked.tif")


# ndvi_rdd_ate <-  raster("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_D_estimate_masked.tif")
# gpp_rdd_ate <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_D_estimate_masked.tif")
# nirv_rdd_ate <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_D_estimate_masked.tif")

map <- read.csv("/home/j_xiao/data/NCdemo/causal_inference_demo1/MeanMAP00-22.csv")[,2]
tmp <- read.csv("/home/j_xiao/data/NCdemo/causal_inference_demo1/MeanTMP00-22.csv")[,2]
pet <- read.csv("/home/j_xiao/data/NCdemo/causal_inference_demo1/MeanPET00-22.csv")[,2]
tp <- read.csv("/home/j_xiao/data/NCdemo/causal_inference_demo1/demo_dataframe_total_P.csv")[,2]
cn <- read.csv("/home/j_xiao/data/NCdemo/causal_inference_demo1/demo_dataframe_cn_topsoil_mean.csv")[,2]
mi <- map/(pet*365) # pet_units=mm/day,pre_units=mm/month,map_units=mm/year

##SM##
rsm <- raster("./causal_inference_demo1/MeanSM_0022.tif")
sm <- as.numeric(as.matrix(rsm))
sm <- sm*0.1

##VPD##
rvpd <- raster("./causal_inference_demo1/MeanVPD_0022.tif")
vpd <- as.numeric(as.matrix(rvpd))
vpd <- vpd*0.01


##AET##
raet <- raster("./causal_inference_demo1/MeanAET_0022.tif")
aet <- as.numeric(as.matrix(raet))
aet <- aet*0.1

rm(rvpd,rsm,raet)

##PAR## 2000/03-2022/12
dirpar <- stack("./data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_direct_clr_mon")
diffpar <- stack("./data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_diff_clr_mon")
totalpar <- dirpar+diffpar
meantotalpar <- mean(totalpar)
meantotalpar <- disaggregate(meantotalpar,fact=2,)
par <-  as.numeric(as.matrix(meantotalpar))

rm(dirpar,diffpar,totalpar,meantotalpar,raet,rsm,rvpd)

##########

df <- data.frame(tmp=tmp,sm=sm,map=map,mi=mi,vpd=vpd,tp=tp,pet=pet,par=par,cn=cn,aet=aet, ATE=as.numeric(as_array(ndvi_cf_ate)))
df <- data.frame(tmp=tmp,sm=sm,map=map,mi=mi,vpd=vpd,tp=tp,pet=pet,par=par,cn=cn,aet=aet, ATE=as.numeric(as_array(gpp_cf_ate)))
df <- data.frame(tmp=tmp,sm=sm,map=map,mi=mi,vpd=vpd,tp=tp,pet=pet,par=par,cn=cn,aet=aet, ATE=as.numeric(as_array(nirv_cf_ate)))

df <- na.omit(df)
df <- round(df,2)
cleaned_df <- remove_outliers_df(df, "ATE")
# cleaned_df <- remove_outliers_df(cleaned_df, "X")
# cleaned_df <- na.omit(cleaned_df)
# ggplot(cleaned_df, aes(x = X, y = ATE)) +
#   geom_point(alpha = 0.5) +  # 散点图
#   geom_smooth(method = "loess") 


summary(cleaned_df$ATE) 
bins <- c(summary(cleaned_df$ATE)["Min."], summary(cleaned_df$ATE)["1st Qu."], 0 ,summary(cleaned_df$ATE)["Max."])
#bins <- c(summary(cleaned_df$ATE)["Min."], summary(cleaned_df$ATE)["Median"], 0 ,summary(cleaned_df$ATE)["Max."])

#labels <- c('Strongly Negative', 'Moderately Negative', 'Positive')


cleaned_df$ATE_group <- cut(cleaned_df$ATE, breaks = bins, labels = labels, include.lowest = TRUE)


table(cleaned_df$ATE_group)

#library(reshape2)
cleaned_df_melted <- melt(cleaned_df, id.vars = c("ATE_group"), 
                          measure.vars = c("tmp", "sm", "map", "mi", "vpd", "tp", "pet", "par", "cn", "aet"))
cleaned_df_melted$value_adjusted <- ifelse(cleaned_df_melted$variable == "tmp", pmax(pmin(cleaned_df_melted$value, 30), 15),
                                           ifelse(cleaned_df_melted$variable == "sm", pmax(pmin(cleaned_df_melted$value, 200), 0),
                                                  ifelse(cleaned_df_melted$variable == "map", pmax(pmin(cleaned_df_melted$value, 2500), 0),
                                                         ifelse(cleaned_df_melted$variable == "mi", pmax(pmin(cleaned_df_melted$value, 2), 0),
                                                                ifelse(cleaned_df_melted$variable == "vpd", pmax(pmin(cleaned_df_melted$value, 2), 0.5),
                                                                       ifelse(cleaned_df_melted$variable == "tp", pmax(pmin(cleaned_df_melted$value, 600), 0),
                                                                              ifelse(cleaned_df_melted$variable == "pet", pmax(pmin(cleaned_df_melted$value, 5), 3),
                                                                                     ifelse(cleaned_df_melted$variable == "par", pmax(pmin(cleaned_df_melted$value, 135), 115),
                                                                                            ifelse(cleaned_df_melted$variable == "cn", pmax(pmin(cleaned_df_melted$value, 1.5), 0.75),
                                                                                                   ifelse(cleaned_df_melted$variable == "aet", pmax(pmin(cleaned_df_melted$value, 100), 30),
                                                                                                          cleaned_df_melted$value))))))))))


library(ggthemes)
library(dplyr)

# single dataset plot
filtered_df <- cleaned_df_melted %>%
  filter(variable %in% c("tmp", "sm", "map", "vpd", "par", "aet"))

ggplot(filtered_df, aes(x = ATE_group, y = value_adjusted, fill = ATE_group)) +
  geom_boxplot(outlier.shape = NA, coef = 0) + 
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # 按不同的环境变量分面展示
  scale_fill_manual(values = c("Positive" = "#FF7F50", "Moderately Negative" = "#BDE4F4", "Strongly Negative" = "#404969")) +  # 自定义颜色
  theme_classic()+

  theme(axis.text.x = element_text(angle = 90, hjust = 1),  # 旋转x轴标签
        strip.text = element_text(size = 8),
        legend.position ="none") +  # 控制分面标签大小
  labs(x = "Average treatment effects", y = "Value")



# conbined dataset plot
combined_df <- dplyr::bind_rows(
  dplyr::mutate(cleaned_df_melted_gpp, dataset = "GPP"),
  dplyr::mutate(cleaned_df_melted_ndvi, dataset = "NDVI"),
  dplyr::mutate(cleaned_df_melted_nirv, dataset = "NIRv")
)
filtered_df <- combined_df %>%
  filter(variable %in% c("tmp", "sm", "map", "vpd", "par", "aet"))


filtered_df$variable <- toupper(filtered_df$variable )
filtered_df$dataset <- factor(filtered_df$dataset, levels = c("NDVI","NIRv","GPP"))
filtered_df$variable <- factor(filtered_df$variable,levels = c("TMP","VPD","SM","MAP","AET","PAR") )

p <- ggplot(filtered_df, aes(x = ATE_group, y = value_adjusted, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Facet by environmental variables
  scale_fill_manual(values = c("GPP" = "#FF7750", "NDVI" = "#BDE4F4", "NIRv" = "#284B63")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "top"  # Place the legend at the top
  ) +
  labs(
    y = "Value",
    x = "ATE Group",
    fill = "Dataset",
    title = "Comparison of NDVI, NIRv, and GPP Across ATE Groups"
  )

ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/ATE_cf_combined_data_facet.pdf",p, width = 10,height = 8, dpi=300)
