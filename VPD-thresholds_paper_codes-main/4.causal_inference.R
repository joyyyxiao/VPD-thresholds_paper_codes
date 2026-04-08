library(grf)
library(raster)
library(dplyr)

####load data
load("suppl_data/data_monthly.RData")

##base functions
zscore <- function(data){
  z_scores <- (data - mean(data, na.rm = TRUE)) / sd(data, na.rm = TRUE)
  return(z_scores)
}

##load threshold data
f_th <- raster("results_mr1/FLUXCOM_threshold2_masked.tif")
g_th <- raster("results_mr1/GLASS_threshold2_masked.tif")
m_th <- raster("results_mr1/MODIS_threshold2_masked.tif")
f_th[f_th == 99] <- NA
g_th[g_th == 99] <- NA
m_th[m_th == 99] <- NA

thres1 <- data.frame(threshold=as.numeric(as.matrix(f_th)))
thres2 <- data.frame(threshold=as.numeric(as.matrix(g_th)))
thres3 <- data.frame(threshold=as.numeric(as.matrix(m_th)))

rownames(thres1) <- lonlat
rownames(thres2) <- lonlat
rownames(thres3) <- lonlat

#load monthly eCO2
meco2 <- read.csv("results_mr1/MODIS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- as.matrix(meco2)
meco2 <- meco2[,1:276]

meco2 <- read.csv("results_mr1/GLASS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- as.matrix(meco2)

meco2 <- read.csv("results_mr1/FLUXCOM_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- as.matrix(meco2)


###build function
###NOTE:check time period

cf_joy_track <- function(line, threshold){
    df <- na.omit(data.frame(month = rep(1:12, ncol(meco2)/12), 
                             eco2 = zscore(meco2[line,]),
                             par = mpar[line,],
                             pre = mpre[line,],
                             sm = msm[line,],
                             tmp = mtmp[line,],
                             vpd = mvpd[line,]))
    df$month <- factor(df$month, levels = 1:12, labels = month.name)
    
    out_null <- data.frame(ATE_all=NA,
                           lower_ci_all=NA,
                           upper_ci_all=NA,
                           ATE_treat=NA,
                           lower_ci_treat=NA,
                           upper_ci_treat=NA,
                           ATE_control=NA,
                           lower_ci_control=NA,
                           upper_ci_control=NA)
    
    if (is.null(threshold[line,]) || any(is.na(threshold[line,]))) {
      return(out_null)
    }
    
    test.vpdth <- threshold[line,]
    
    if (!is.numeric(df$vpd) || !is.numeric(test.vpdth)) {
      return(out_null)
    }
    
    df$vpd <- ifelse(df$vpd < test.vpdth, 0, 1)
    
    if (any(is.na(df)))
      return(out_null)
    
    set.seed(123)
    X <- model.matrix(~ month + par + pre + sm + tmp, data = df)
#   X <- df[,2:5]
    W <- df$vpd
    forest.W <- regression_forest(X, W, tune.parameters = "all")
    W.hat <- predict(forest.W)$predictions
    
    ### Filter
    overlap0 <- df %>%
      mutate(w.hat = W.hat) %>%
      mutate(keep = case_when(w.hat < 0.05 | w.hat > 0.95 ~ 0,
                              w.hat >= 0.05 & w.hat <= 0.95  ~ 1))
    overlap <- overlap0 %>% filter(keep == 1)
    
    if (nrow(overlap) < 80)
      return(out_null)
    
    X <- model.matrix(~ month + par + pre + sm + tmp, data = overlap)
  #X <- as.matrix(overlap[,2:5])
    Y <- overlap$eco2
    W <- overlap$vpd
    
    set.seed(123)
    forest.Y <- regression_forest(X, Y, tune.parameters = "all")
    Y.hat <- predict(forest.Y)$predictions
    W.hat <- overlap$w.hat
    
    tau.forest <- causal_forest(X, Y, W, W.hat = W.hat, Y.hat = Y.hat, tune.parameters = "all")
    ate_cf_all <- average_treatment_effect(tau.forest, target.sample = "all")
    ate_cf_treat <- average_treatment_effect(tau.forest, target.sample = "treat")
    ate_cf_control <- average_treatment_effect(tau.forest, target.sample = "control")
    
    out <- data.frame(ATE_all = ate_cf_all["estimate"],
                      lower_ci_all = ate_cf_all["estimate"] - 1.96 * ate_cf_all["std.err"],
                      upper_ci_all = ate_cf_all["estimate"] + 1.96 * ate_cf_all["std.err"],
                      ATE_treat = ate_cf_treat["estimate"],
                      lower_ci_treat = ate_cf_treat["estimate"] - 1.96 * ate_cf_treat["std.err"],
                      upper_ci_treat = ate_cf_treat["estimate"] + 1.96 * ate_cf_treat["std.err"],
                      ATE_control = ate_cf_control["estimate"],
                      lower_ci_control = ate_cf_control["estimate"] - 1.96 * ate_cf_control["std.err"],
                      upper_ci_control = ate_cf_control["estimate"] + 1.96 * ate_cf_control["std.err"])
    
    return(out)

}
cf_joy_track2 <- function(line, threshold) {
  log_entries <- character()  # Buffer for logging
  # Initialize result as out_null at the beginning
  result <- data.frame(ATE_all = NA, lower_ci_all = NA, upper_ci_all = NA,
                       ATE_treat = NA, lower_ci_treat = NA, upper_ci_treat = NA,
                       ATE_control = NA, lower_ci_control = NA, upper_ci_control = NA)
  result <- tryCatch({
    # Start processing log entry
    log_entry <- paste("Start processing line:", line, "at", Sys.time(), "\n")
    log_entries <- c(log_entries, log_entry)
    
    df <- na.omit(data.frame(
      month = rep(1:12, ncol(meco2)/12),  #2000-2022 23years or 2000-2018 19years
      eco2 = zscore(as.numeric(meco2[line,])),
      par = mpar[line,],
      pre = mpre[line,],
      sm = msm[line,],
      tmp = mtmp[line,],
      vpd = mvpd[line,]
    ))
    df$month <- factor(df$month, levels = 1:12, labels = month.name)
    
    out_null <- data.frame(ATE_all = NA, lower_ci_all = NA, upper_ci_all = NA,
                           ATE_treat = NA, lower_ci_treat = NA, upper_ci_treat = NA,
                           ATE_control = NA, lower_ci_control = NA, upper_ci_control = NA)
    
    # Custom warning checks
    if (is.null(threshold[line,]) || any(is.na(threshold[line,]))) {
      warning("Threshold is NULL or contains NA")
      return(out_null)
    }
    
    test.vpdth <- threshold[line,]
    
    if (!is.numeric(df$vpd) || !is.numeric(test.vpdth)) {
      warning("Non-numeric value found in df$vpd or test.vpdth")
      return(out_null)
    }
    
    df$vpd <- ifelse(df$vpd < test.vpdth, 0, 1)
    
    if (any(is.na(df))) {
      warning("NA values found in df after processing")
      return(out_null)
    }
    
    # Process overlap
    X <- model.matrix(~ month + par + pre + sm + tmp, data = df)
    W <- df$vpd
    forest.W <- regression_forest(X, W, tune.parameters = "all")
    W.hat <- predict(forest.W)$predictions
    
    overlap <- df %>%
      mutate(w.hat = W.hat) %>%
      mutate(keep = case_when(w.hat < 0.05 | w.hat > 0.95 ~ 0,
                              w.hat >= 0.05 & w.hat <= 0.95 ~ 1)) %>%
      filter(keep == 1)
    
    if (nrow(overlap) < 80) {
      warning("Not enough overlap found")
      return(out_null)
    }
    
    # Prediction
    X <- model.matrix(~ month + par + pre + sm + tmp, data = overlap)
    Y <- overlap$eco2
    W <- overlap$vpd
    
    forest.Y <- regression_forest(X, Y, tune.parameters = "all")
    Y.hat <- predict(forest.Y)$predictions
    W.hat <- overlap$w.hat
    
    tau.forest <- causal_forest(X, Y, W, W.hat = W.hat, Y.hat = Y.hat, tune.parameters = "all")
    ate_cf_all <- average_treatment_effect(tau.forest, target.sample = "all")
    ate_cf_treat <- average_treatment_effect(tau.forest, target.sample = "treat")
    ate_cf_control <- average_treatment_effect(tau.forest, target.sample = "control")
    
    out <- data.frame(ATE_all = ate_cf_all["estimate"],
                      lower_ci_all = ate_cf_all["estimate"] - 1.96 * ate_cf_all["std.err"],
                      upper_ci_all = ate_cf_all["estimate"] + 1.96 * ate_cf_all["std.err"],
                      ATE_treat = ate_cf_treat["estimate"],
                      lower_ci_treat = ate_cf_treat["estimate"] - 1.96 * ate_cf_treat["std.err"],
                      upper_ci_treat = ate_cf_treat["estimate"] + 1.96 * ate_cf_treat["std.err"],
                      ATE_control = ate_cf_control["estimate"],
                      lower_ci_control = ate_cf_control["estimate"] - 1.96 * ate_cf_control["std.err"],
                      upper_ci_control = ate_cf_control["estimate"] + 1.96 * ate_cf_control["std.err"])
    
    return(out)  # Successfully return the result
    
  }, warning = function(w) {
    # Check if the warning is custom by checking the message
    custom_warnings <- c("Threshold is NULL or contains NA",
                         "Non-numeric value found in df$vpd or test.vpdth",
                         "NA values found in df after processing",
                         "Not enough overlap found")
    
    if (conditionMessage(w) %in% custom_warnings) {
      log_entries <<- c(log_entries, paste("Custom Warning in line", line, ":", w$message, "\n"))
      return(out_null)  # Return out_null if a custom warning occurs
    } else {
      log_entries <<- c(log_entries, paste("Non-Custom Warning in line", line, ":", w$message, "\n"))
      return(result)  # Return result if it's a non-custom warning
    }
  }, error = function(e) {
    log_entries <<- c(log_entries, paste("Error in line", line, ":", e$message, "\n"))
    return(out_null)  # Return NA if an error occurs
  }, finally = {
    # Log everything at the end of the process
    if (length(log_entries) > 0) {
      cat(log_entries, file = "process_log_re.txt", append = TRUE)
    }
  })
  
  return(result)
}



###loop

library(foreach)
library(doSNOW)
library(doRNG)
library(parallel)

cl <- makeCluster(4)
registerDoSNOW(cl)
registerDoRNG(seed = 123)

# Setup the progress bar
pb <- txtProgressBar(max = length(lonlat), style = 3)
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress) 

print(paste("Starting parellel processing at:", Sys.time()))
ptime  <- system.time(
  cf_result_fluxcom <- foreach(
    line = lonlat,
    .options.snow = opts,
#   .errorhandling = "pass", 
    .combine = rbind,
    .packages = c('grf', 'dplyr')
  ) %dopar% {
    cf_joy_track2(line,threshold=thres1)
  }
)

print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)

####save data
write.csv(cf_result_modis,"results_mr1/MODIS_cf_result.csv")
write.csv(cf_result_glass,"results_mr1/GLASS_cf_result.csv")
write.csv(cf_result_fluxcom,"results_mr1/FLUXCOM_cf_result.csv")


########Plot
library(ggplot2)
library(ggthemes)
library(reshape2)

as_raster <- function(x){
  mat <- matrix(x,nrow=360,ncol=720)
  r <- raster(mat, xmn=-180,xmx=180,ymn=-90,ymx=90,crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  return(r)
}
remove_except <- function(keep) {
  
  all_objects <- ls(envir = .GlobalEnv)
  to_remove <- setdiff(all_objects, keep)
  rm(list = to_remove, envir = .GlobalEnv)
  cat("Kept objects:", keep, "\n")
}
remove_outliers_df <- function(data, column) {
  
  Q1 <- quantile(data[[column]], 0.25, na.rm = TRUE)
  Q3 <- quantile(data[[column]], 0.75, na.rm = TRUE)
  
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  cleaned_data <- data[data[[column]] >= lower_bound & data[[column]] <= upper_bound, ]
  
  return(cleaned_data)
}
remove_outliers <- function(data) {
  if (!is.numeric(data)) {
    stop("Please provide a numeric vector.")
  }
  
  Q1 <- quantile(data, 0.25, na.rm = TRUE)
  Q3 <- quantile(data, 0.75, na.rm = TRUE)
  
  IQR <- Q3 - Q1
  
  lower_bound <- Q1 - 1.5 * IQR
  upper_bound <- Q3 + 1.5 * IQR
  
  data[data < lower_bound | data > upper_bound] <- NA
  
  return(data)
}


ate_all_f <- as_raster(remove_outliers(cf_result_fluxcom$ATE_all))
ate_all_g <- as_raster(remove_outliers(cf_result_glass$ATE_all))
ate_all_m <- as_raster(remove_outliers(cf_result_modis$ATE_all))

mean_data <- read.csv("suppl_data/data_multiyearmean.csv")

df <- data.frame(tmp=mean_data$TMP,sm=mean_data$SM,map=mean_data$MAP,mi=mean_data$MI,vpd=mean_data$VPD,tp=mean_data$P,pet=mean_data$PET,par=mean_data$PAR,cn=mean_data$CNr,aet=mean_data$AET, ATE=remove_outliers(cf_result_fluxcom$ATE_all))
df <- data.frame(tmp=mean_data$TMP,sm=mean_data$SM,map=mean_data$MAP,mi=mean_data$MI,vpd=mean_data$VPD,tp=mean_data$P,pet=mean_data$PET,par=mean_data$PAR,cn=mean_data$CNr,aet=mean_data$AET, ATE=remove_outliers(cf_result_glass$ATE_all))
df <- data.frame(tmp=mean_data$TMP,sm=mean_data$SM,map=mean_data$MAP,mi=mean_data$MI,vpd=mean_data$VPD,tp=mean_data$P,pet=mean_data$PET,par=mean_data$PAR,cn=mean_data$CNr,aet=mean_data$AET, ATE=remove_outliers(cf_result_modis$ATE_all))

df <- na.omit(df)
df <- round(df,2)

summary(df$ATE) 
bins <- c(summary(df$ATE)["Min."], summary(df$ATE)["1st Qu."], 0 ,summary(df$ATE)["Max."])
labels <- c('Strongly Negative', 'Moderately Negative', 'Positive')
df$ATE_group <- cut(df$ATE, breaks = bins, labels = labels, include.lowest = TRUE)
table(df$ATE_group)


df_melted <- melt(df, id.vars = c("ATE_group"), 
                          measure.vars = c("tmp", "sm", "map", "mi", "vpd", "tp", "pet", "par", "cn", "aet"))
df_melted$value_adjusted <- ifelse(df_melted$variable == "tmp", pmax(pmin(df_melted$value, 30), 15),
                                           ifelse(df_melted$variable == "sm", pmax(pmin(df_melted$value, 200), 0),
                                                  ifelse(df_melted$variable == "map", pmax(pmin(df_melted$value, 2500), 0),
                                                         ifelse(df_melted$variable == "mi", pmax(pmin(df_melted$value, 2), 0),
                                                                ifelse(df_melted$variable == "vpd", pmax(pmin(df_melted$value, 2), 0.5),
                                                                       ifelse(df_melted$variable == "tp", pmax(pmin(df_melted$value, 600), 0),
                                                                              ifelse(df_melted$variable == "pet", pmax(pmin(df_melted$value, 5), 3),
                                                                                     ifelse(df_melted$variable == "par", pmax(pmin(df_melted$value, 135), 115),
                                                                                            ifelse(df_melted$variable == "cn", pmax(pmin(df_melted$value, 1.5), 0.75),
                                                                                                   ifelse(df_melted$variable == "aet", pmax(pmin(df_melted$value, 100), 30),
                                                                                                          df_melted$value))))))))))
df_melted_f <- df_melted
df_melted_g <- df_melted
df_melted_m <- df_melted

# conbined dataset plot
combined_df <- dplyr::bind_rows(
  dplyr::mutate(df_melted_f, dataset = "FLUXCOM GPP"),
  dplyr::mutate(df_melted_g, dataset = "GLASS GPP"),
  dplyr::mutate(df_melted_m, dataset = "MODIS GPP")
)
filtered_df <- combined_df %>%
  filter(variable %in% c("tmp", "sm", "map", "vpd", "par", "aet"))

filtered_df$variable <- toupper(filtered_df$variable )
filtered_df$dataset <- factor(filtered_df$dataset, levels = c("FLUXCOM GPP","GLASS GPP","MODIS GPP"))
filtered_df$variable <- factor(filtered_df$variable,levels = c("TMP","VPD","SM","MAP","AET","PAR") )

#Fig.4
p <- ggplot(filtered_df, aes(x = ATE_group, y = value_adjusted, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +  # Facet by environmental variables
  scale_fill_manual(values = c("FLUXCOM GPP" = "#FF7750", "GLASS GPP" = "#BDE4F4", "MODIS GPP" = "#284B63")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(size = 8),
    legend.position = "top"  # Place the legend at the top
  ) +
  labs(
    y = "Value",
    x = NULL,
    fill = "Products",
  )

ggsave("plots/ATE_cf_combined_data_facet.pdf",p, width = 10,height = 8, dpi=300)
