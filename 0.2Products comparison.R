library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
library(ggthemes)
#load data
GPP_ann_df  <- read.csv("GPPdatasets/FLUXCOM_GPP_annual_dataframe.csv",row.names=1)
GPP2_ann_df  <- read.csv("GPPdatasets/GLASS_GPP_annual_dataframe.csv",row.names=1)
GPP3_ann_df  <- read.csv("GPPdatasets/MODIS_GPP_annual_dataframe.csv",row.names=1)

years1 <- 2000:2018
years2 <- 2000:2025

colnames(GPP_ann_df) <- years1
colnames(GPP2_ann_df) <- years1
colnames(GPP3_ann_df) <- years2


#region data(north/tropics/south)
region3 <- raster("suppl_data/transcom_reg3.tif")
reg3 <- as.numeric(as.matrix(region3))
GPP_ann_df$region <- reg3
GPP2_ann_df$region <- reg3
GPP3_ann_df$region <- reg3

##########by_pixel

##by latitude
years <- 2000:2018
years <- 2000:2025

iav_res <- apply(GPP_ann_df, 1, function(y) {
  
  if (sum(!is.na(y)) < 18) {
    return(c(IAV = NA))
  }
  
  fit <- lm(y ~ years)
  resid <- residuals(fit)
  c(IAV = sd(resid, na.rm = TRUE))
})

#iav_map_gpp <- as_raster(iav_res)
GPP1_iav_res <- data.frame(iav=iav_res,region=reg3)
GPP2_iav_res <- data.frame(iav=iav_res,region=reg3)
GPP3_iav_res <- data.frame(iav=iav_res,region=reg3)

iav_res_gpp_bylat3 <- GPP_iav_res3 %>% filter(!is.na(region)) %>% 
                     mutate(coord = rownames(.),
                     latitude = as.numeric(sub(".*?,\\s*([^)]+)\\)", "\\1", coord))) %>%
                     group_by(latitude) %>%
  summarise(
    mean_slope = mean(iav, na.rm = TRUE),
    sd_slope   = sd(iav, na.rm = TRUE),
    n_pixel = dplyr::n(),
    se_slope   = sd_slope / sqrt(n_pixel),
    ci_lower   = mean_slope - 1.96 * se_slope,
    ci_upper   = mean_slope + 1.96 * se_slope
  ) %>%
  arrange(latitude)

gpp1_plot <- iav_res_gpp_bylat1 %>%
             mutate(product = "FLUXCOM")
gpp2_plot <- iav_res_gpp_bylat2 %>%
  mutate(product = "GLASS")
gpp3_plot <- iav_res_gpp_bylat3 %>%
  mutate(product = "MODIS")

lat_all <- bind_rows(gpp1_plot, gpp2_plot,gpp3_plot)

#Fig.S8
p <- ggplot(lat_all,
            aes(x = latitude,
                y = mean_slope,
                color = product,
                fill = product)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.2,
              colour = NA) +
  geom_line(linewidth = 1) +
  scale_fill_manual(values = c("#845EC2","#4e8397","#d5cabd"))+
  scale_color_manual(values = c("#845EC2","#4e8397","#d5cabd"))+
  theme_few() +
  labs(
    x = "Latitude (0.5° grid)",
    y = "Mean IAV",
    color = "Product",
    fill = "Product")+
  theme(
    legend.title = element_blank(),          
    legend.position = c(0.85, 0.75),         
    legend.background = element_blank()
  )

ggsave("plots/iav_bylat_3p.pdf",p,width = 8,height = 3, dpi=300)



#################BY REGION###################
biome <- raster("suppl_data/SYNMAP_biomass_4type_2000.tif")
climate <- raster("suppl_data/climateClass_KTC_fromCRU_2002to2021.tif")

GPP1_iav_res$biome <- as.numeric(as.matrix(biome))
GPP1_iav_res$climate <- as.numeric(as.matrix(climate))
GPP2_iav_res$biome <- as.numeric(as.matrix(biome))
GPP2_iav_res$climate <- as.numeric(as.matrix(climate))
GPP3_iav_res$biome <- as.numeric(as.matrix(biome))
GPP3_iav_res$climate <- as.numeric(as.matrix(climate))

GPP3_by_reg <- GPP3_iav_res %>%
  filter(!is.na(region)) %>%                     
  group_by(region) %>%                    
  summarise(
    mean_IAV = mean(iav, na.rm = TRUE),           
    n_pixel  = sum(!is.na(iav)), 
    sd = sd(iav, na.rm = TRUE),
    se = sd / sqrt(n_pixel),
    .groups = "drop"
  )%>%ungroup()

GPP3_by_cli<- GPP3_iav_res %>%
  filter(!is.na(climate) &climate != -1&climate != 0&climate != 60) %>%                     
  group_by(climate) %>%                    
  summarise(
    mean_IAV = mean(iav, na.rm = TRUE),           
    n_pixel  = sum(!is.na(iav)), 
    sd = sd(iav, na.rm = TRUE),
    se = sd / sqrt(n_pixel),
    .groups = "drop"
  )%>%ungroup()

GPP3_by_bio <- GPP3_iav_res %>%
  filter(!is.na(biome)) %>%                     
  group_by(biome) %>%                    
  summarise(
    mean_IAV = mean(iav, na.rm = TRUE),           
    n_pixel  = sum(!is.na(iav)), 
    sd = sd(iav, na.rm = TRUE),
    se = sd / sqrt(n_pixel),
    .groups = "drop"
  )%>%ungroup()


#plot
GPP1_by_cli$product ="FLUXCOM"
GPP1_by_bio$product ="FLUXCOM"
GPP1_by_reg$product ="FLUXCOM"
GPP2_by_cli$product ="GLASS"
GPP2_by_bio$product ="GLASS"
GPP2_by_reg$product ="GLASS"
GPP3_by_cli$product ="MODIS"
GPP3_by_bio$product ="MODIS"
GPP3_by_reg$product ="MODIS"

all_by_cli <- rbind(GPP1_by_cli,GPP2_by_cli,GPP3_by_cli)
all_by_bio <- rbind(GPP1_by_bio,GPP2_by_bio,GPP3_by_bio)
all_by_reg <- rbind(GPP1_by_reg,GPP2_by_reg,GPP3_by_reg)


#Fig.S10
p <- ggplot(all_by_cli,
            aes(x = factor(climate),
                y = mean_IAV,
                fill = factor(product))) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.7) +
  geom_errorbar(aes(ymin = mean_IAV - se,
                    ymax = mean_IAV + se),
                width = 0.2,
                linewidth = 0.4,
                position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#845EC2","#4e8397","#d5cabd"))+
  theme_few() +
  labs(
    x = "Climate Region",
    y = "Mean IAV"
  )+
  scale_x_discrete(
    labels = c(
      "10" = "Tropical humid climates",
      "20" = "Dry climates",
      "30" = "Subtropical climates",
      "40" = "Temperate climates",
      "50" = "Boreal climates"
    ))+
  theme(
    legend.title = element_blank(),          
    legend.position = c(0.85, 0.85),         
    legend.background = element_blank()
  )

ggsave("plots/iav_bycli_3p.pdf",p,width = 8,height = 4.5, dpi=300)

#Fig.S11
p <- ggplot(all_by_bio,
            aes(x = factor(biome),
                y = mean_IAV,
                fill = product)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.7) +
  geom_errorbar(aes(ymin = mean_IAV - se,
                    ymax = mean_IAV + se),
                width = 0.2,
                linewidth = 0.4,
                position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#845EC2","#4e8397","#d5cabd"))+
  theme_few() +
  labs(
    x = "Biome",
    y = "Mean IAV"
  )+
  scale_x_discrete(
    labels = c(
      "10" = "Trees",
      "20" = "Shrubs",
      "30" = "Grasses",
      "40" = "Crops"
    ))+
  theme(
    legend.title = element_blank(),          
    legend.position = c(0.35, 0.85),         
    legend.background = element_blank()
  )

ggsave("plots/iav_bybio_3p.pdf",p,width = 8,height = 4.5, dpi=300)

#Fig. S9
p <- ggplot(all_by_reg,
            aes(x = factor(region),
                y = mean_IAV,
                fill = product)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.7) +
  geom_errorbar(aes(ymin = mean_IAV - se,
                    ymax = mean_IAV + se),
                width = 0.2,
                linewidth = 0.4,
                position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("#845EC2","#4e8397","#d5cabd"))+
  theme_few() +
  labs(
    x = "TransCom Region",
    y = "Mean IAV"
  )+
  scale_x_discrete(
    labels = c(
      "33" = "North",
      "66" = "Tropics",
      "99" = "South"
    ))+
  theme(
    legend.title = element_blank(),          
    legend.position = c(0.15, 0.85),         
    legend.background = element_blank()
  )

ggsave("plots/iav_byreg_3p.pdf",p,width = 8,height = 4.5, dpi=300)
