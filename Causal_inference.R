library(grf)
library(raster)
library(dplyr)

####load data
setwd("/home/j_xiao/data/NCdemo")
lonlat <- read.csv("/home/j_xiao/data/NCdemo/line.csv")[,2]


##VPD##
file <- paste0("/home/j_xiao/data/TerraClimate/TerraVPD/","vpd",rep(2000:2022,each=12),"_",1:12,".tif")
file <- paste0("/home/j_xiao/data/TerraClimate/TerraVPD/","vpd",rep(2000:2018,each=12),"_",1:12,".tif")
rvpd <- stack(file)
mvpd <- as_array(rvpd)
rm(rvpd,file)
mvpd <- mvpd*0.01

##TMP##
rtmp <- stack("/home/j_xiao/data/CRUclimate/cru_ts4.06.1981.1990.tmp.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.1991.2000.tmp.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.2001.2010.tmp.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.2011.2020.tmp.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.07.2021.2022.tmp.dat.nc")
mtmp <- as_array(rtmp)
rm(rtmp)
mtmp <- mtmp[,229:504] #00-22

##PRE##
rpre <- stack("/home/j_xiao/data/CRUclimate/cru_ts4.06.1981.1990.pre.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.1991.2000.pre.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.2001.2010.pre.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.06.2011.2020.pre.dat.nc",
              "/home/j_xiao/data/CRUclimate/cru_ts4.07.2021.2022.pre.dat.nc")
mpre <- as_array(rpre)
rm(rpre)
mpre <- mpre[,229:504] #00-22

##ERA5_SM##
rsm <- stack("/home/j_xiao/data/ERA5sm0_100_global/sm_aggr_720*360_00-22.tif")
msm <- as_array(rsm)
rm(rsm)

##CERES_SYN1deg_PAR##
dirpar <- stack("./data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_direct_clr_mon")
diffpar <- stack("./data/CERES_SYN1deg-Month_Terra-Aqua-MODIS_Ed4.1_Subset_200003-202212.nc",varname="adj_sfc_par_diff_clr_mon")
totalpar <- dirpar+diffpar
totalpar <- disaggregate(totalpar,fact=2)
mpar <- as_array(totalpar)
rm(dirpar,diffpar,totalpar)
mpar <- cbind(rep(NA,259200),rep(NA,259200),mpar)


##eCO2##
meco2  <- read.csv("./results/GIMMS_NDVI/GIMMS_NDVI_monthly_ppm_CO2_effect_nainterp.csv",row.names = 1) #1982-2022 meco2 <- meco2[,217:492] #from2000
meco2 <- read.csv("./results/MODIS_NIRv/NIRv_monthly_ppm_CO2_effect_nainterp.csv",row.names = 1) #meco2 <- cbind(rep(NA,259200),meco2)
meco2 <- read.csv("./results/FLUXCOM_GPP/FLUXCOM_GPP_monthly_ppm_CO2_effect_nainterp.csv",row.names = 1) #1980-2018 meco2 <- meco2[,241:468] #from2000
# mpar <- mpar[, 1:228]
# mtmp <- mtmp[, 1:228]
# msm <- msm[, 1:228]
# mvpd <- mvpd[, 1:228]
# mpre <- mpre[, 1:228]

rownames(meco2) <- lonlat
colnames(meco2) <- 1:276
meco2 <- as.matrix(meco2)

##load threshold data
ndvi_thres <- read.csv("/home/j_xiao/data/NCdemo/results/GIMMS_NDVI/NDVI_threshold2_h2o_gbm.csv")[-1,]
ndvi_thres <- ndvi_thres$threshold
ndvi_thres <- data.frame(thres=ndvi_thres)
rownames(ndvi_thres) <- lonlat

nirv_thres <- read.csv("/home/j_xiao/data/NCdemo/results/MODIS_NIRv/NIRv_threshold2_h2o_gbm.csv")[-1,]
nirv_thres <- nirv_thres$threshold
nirv_thres <- data.frame(thres=nirv_thres)
rownames(nirv_thres) <- lonlat

gpp_thres <- read.csv("/home/j_xiao/data/NCdemo/results/FLUXCOM_GPP/GPP_threshold2_h2o_gbm.csv")[-1,]
gpp_thres <- gpp_thres$threshold
gpp_thres <- data.frame(thres=gpp_thres)
rownames(gpp_thres) <- lonlat



##test case
line <- "(80.75, 21.25)"
X <- data.frame(par=mpar[line,],
                pre=mpre[line,],
                sm=msm[line,],
                tmp=mtmp[line,])

X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)

test.vpd <- mvpd[line,]
test.vpdth <- ndvi_thres[line,]
W <- ifelse(test.vpd < test.vpdth, 0, 1)

Y <- zscore(meco2[line,])

tau.forest <- causal_forest(X, Y, W)
tau.forest
#> GRF forest object of type causal_forest 
#> Number of trees: 2000 
#> Number of training samples: 2000 
#> Variable importance: 
#>     1     2     3     4     5     6     7     8     9    10 
#> 0.691 0.032 0.040 0.048 0.030 0.029 0.032 0.030 0.029 0.041
#> 
tau.hat.oob <- predict(tau.forest)
hist(tau.hat.oob$predictions)

#Estimate treatment effects for the test sample
tau.hat <- predict(tau.forest, X.test)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 2)

#Estimate the conditional average treatment effect on the full sample (CATE)
average_treatment_effect(tau.forest, target.sample = "all")
#>   estimate    std.err 
#> 0.36741777 0.04986737

# Estimate the conditional average treatment effect on the treated sample (ATET).
average_treatment_effect(tau.forest, target.sample = "treated")

# Estimate the conditional average treatment effect on the treated sample (ATEC).
average_treatment_effect(tau.forest, target.sample = "control")

# Add confidence intervals for heterogeneous treatment effects; growing more trees is now recommended.
tau.forest <- causal_forest(X, Y, W, num.trees = 4000)
tau.hat <- predict(tau.forest, X.test, estimate.variance = TRUE)
sigma.hat <- sqrt(tau.hat$variance.estimates)
plot(X.test[, 1], tau.hat$predictions, ylim = range(tau.hat$predictions + 1.96 * sigma.hat, tau.hat$predictions - 1.96 * sigma.hat, 0, 2), xlab = "x", ylab = "tau", type = "l")
lines(X.test[, 1], tau.hat$predictions + 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], tau.hat$predictions - 1.96 * sigma.hat, col = 1, lty = 2)
lines(X.test[, 1], pmax(0, X.test[, 1]), col = 2, lty = 1)

# In some examples, pre-fitting models for Y and W separately may
# be helpful (e.g., if different models use different covariates).
# In some applications, one may even want to get Y.hat and W.hat
# using a completely different method (e.g., boosting).

# Generate new data.
n <- 4000
p <- 20
X <- matrix(rnorm(n * p), n, p)
TAU <- 1 / (1 + exp(-X[, 3]))
W <- rbinom(n, 1, 1 / (1 + exp(-X[, 1] - X[, 2])))
Y <- pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + W * TAU + rnorm(n)

forest.W <- regression_forest(X, W, tune.parameters = "all")
W.hat <- predict(forest.W)$predictions

forest.Y <- regression_forest(X, Y, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

forest.Y.varimp <- variable_importance(forest.Y)

# Note: Forests may have a hard time when trained on very few variables
# (e.g., ncol(X) = 1, 2, or 3). We recommend not being too aggressive
# in selection.
selected.vars <- which(forest.Y.varimp / mean(forest.Y.varimp) > 0.2)

tau.forest <- causal_forest(X[, selected.vars], Y, W,
                            W.hat = W.hat, Y.hat = Y.hat,
                            tune.parameters = "all")

tau.forest <- causal_forest(X, Y, W,
                            W.hat = W.hat, Y.hat = Y.hat,
                            tune.parameters = "all")

# See if a causal forest succeeded in capturing heterogeneity by plotting
# the TOC and calculating a 95% CI for the AUTOC.
train <- sample(1:n, n / 2)
train.forest <- causal_forest(X[train, ], Y[train], W[train])
eval.forest <- causal_forest(X[-train, ], Y[-train], W[-train])
rate <- rank_average_treatment_effect(eval.forest,
                                      predict(train.forest, X[-train, ])$predictions)
plot(rate)
paste("AUTOC:", round(rate$estimate, 2), "+/", round(1.96 * rate$std.err, 2))


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

cl <- makeCluster(8)
registerDoSNOW(cl)
registerDoRNG(seed = 123)

# Setup the progress bar
pb <- txtProgressBar(max = length(lonlat), style = 3)
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress) 

print(paste("Starting parellel processing at:", Sys.time()))
ptime  <- system.time(
  cf_result_ndvi <- foreach(
    line = lonlat,
    .options.snow = opts,
#   .errorhandling = "pass", 
    .combine = rbind,
    .packages = c('grf', 'dplyr')
  ) %dopar% {
    cf_joy_track2(line,threshold=ndvi_thres)
  }
)

print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)


write.csv(cf_result_ndvi,"/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_result.csv")


####save data
library(raster)
as_array <- function(raster){
  matrix <- matrix(NA,259200,nlayers(raster))
  for(i in 1:nlayers(raster)){
    matrix[,i] <-   as.numeric(as.matrix(subset(raster,i)))
  }
  lonlat <- read.csv("/home/j_xiao/data/NCdemo/line.csv")[,2]
  rownames(matrix) <- lonlat
  colnames(matrix) <- rep(NA,nlayers(raster))
  return(matrix)
}
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

cf_result_gpp <- read.csv("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_result.csv")[,-1]
thres_mask <- raster("/home/j_xiao/data/NCdemo/results/FLUXCOM_GPP/vpd_threshold_masked.tif")

ate_all <- remove_outliers(cf_result_gpp$ATE_all)
ate_all <- as_raster(ate_all)
ate_all <- mask(ate_all,thres_mask)

ate_control <- remove_outliers(cf_result_gpp$ATE_control)
ate_control <- as_raster(ate_control)
ate_control <- mask(ate_control,thres_mask)

ate_treat <- remove_outliers(cf_result_gpp$ATE_treat)
ate_treat <- as_raster(ate_treat)
ate_treat <- mask(ate_treat,thres_mask)

plot(ate_treat)

writeRaster(ate_all,"/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_all_masked.tif",overwrite=TRUE)
writeRaster(ate_control,"/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_control_masked.tif",overwrite=TRUE)
writeRaster(ate_treat,"/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_treat_masked.tif",overwrite=TRUE)


######
cf_result_ndvi <- read.csv("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_result.csv")[,-1]
thres_mask <- raster("/home/j_xiao/data/NCdemo/results/GIMMS_NDVI/vpd_threshold_masked.tif")

ate_all <- remove_outliers(cf_result_ndvi$ATE_all)
ate_all <- as_raster(ate_all)
ate_all <- mask(ate_all,thres_mask)

ate_control <- remove_outliers(cf_result_ndvi$ATE_control)
ate_control <- as_raster(ate_control)
ate_control <- mask(ate_control,thres_mask)

ate_treat <- remove_outliers(cf_result_ndvi$ATE_treat)
ate_treat <- as_raster(ate_treat)
ate_treat <- mask(ate_treat,thres_mask)

plot(ate_all)
plot(ate_treat)
plot(ate_control)

writeRaster(ate_all,"/home/j_xiao/data/NCdemo/results/causal_inference/ndvi_cf_ate_all_masked.tif",overwrite=TRUE)
writeRaster(ate_control,"/home/j_xiao/data/NCdemo/results/causal_inference/ndvi_cf_ate_control_masked.tif",overwrite=TRUE)
writeRaster(ate_treat,"/home/j_xiao/data/NCdemo/results/causal_inference/ndvi_cf_ate_treat_masked.tif",overwrite=TRUE)


######
cf_result_nirv <- read.csv("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_cf_result.csv")[,-1]
thres_mask <- raster("/home/j_xiao/data/NCdemo/results/MODIS_NIRv/vpd_threshold_masked.tif")

ate_all <- remove_outliers(cf_result_nirv$ATE_all)
ate_all <- as_raster(ate_all)
ate_all <- mask(ate_all,thres_mask)

ate_control <- remove_outliers(cf_result_nirv$ATE_control)
ate_control <- as_raster(ate_control)
ate_control <- mask(ate_control,thres_mask)

ate_treat <- remove_outliers(cf_result_nirv$ATE_treat)
ate_treat <- as_raster(ate_treat)
ate_treat <- mask(ate_treat,thres_mask)

plot(ate_all)
plot(ate_treat)
plot(ate_control)

writeRaster(ate_all,"/home/j_xiao/data/NCdemo/results/causal_inference/nirv_cf_ate_all_masked.tif",overwrite=TRUE)
writeRaster(ate_control,"/home/j_xiao/data/NCdemo/results/causal_inference/nirv_cf_ate_control_masked.tif",overwrite=TRUE)
writeRaster(ate_treat,"/home/j_xiao/data/NCdemo/results/causal_inference/nirv_cf_ate_treat_masked.tif",overwrite=TRUE)


####Plot
library(raster)
library(ggplot2)
library(ggthemes)
library(dplyr)
#GPP
ate_all <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_all_masked.tif")
ate_control <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_control_masked.tif")
ate_treat <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_cf_ate_treat_masked.tif")

#NDVI
ate_all <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_ate_all_masked.tif")
ate_control <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_ate_control_masked.tif")
ate_treat <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_cf_ate_treat_masked.tif")

#NIRv
ate_all <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_cf_ate_all_masked.tif")
ate_control <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_cf_ate_control_masked.tif")
ate_treat <- raster("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_cf_ate_treat_masked.tif")





df <- data.frame(ATE=as.numeric(as_array(ate_all)),
                 ATEC=as.numeric(as_array(ate_control)),
                 ATET=as.numeric(as_array(ate_treat)))


data <- tidyr::pivot_longer(df, 
                            cols = everything(),  
                            names_to = "Category", 
                            values_to = "Value")   

p1 <- 
  ggplot(data,aes(x=Category,y=Value,fill=Category)) +
  geom_boxplot(width=0.5,size=0.4,alpha=1,outlier.alpha = 0.25)+
  geom_hline(yintercept = 0,linetype = "dashed",size=0.3)+
  scale_fill_manual(values=c("#B8001F","#FCFAEE","#507687"))+
  theme_few()+
  theme(legend.position = "none")+
  ylab("Estimate")+
  xlab(element_blank())

ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_ate_allcontreat.pdf",p1, width = 4,height = 5, dpi=300)
ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_ate_allcontreat.pdf",p1, width = 4,height = 5, dpi=300)
ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_ate_allcontreat.pdf",p1, width = 4,height = 5, dpi=300)


summary_data <- na.omit(data) %>%
  group_by(Category) %>%
  summarize(
    mean_Value = mean(Value),
    ci_lower = mean(Value) - qt(0.975, df=n()-1) * sd(Value) / sqrt(n()),  
    ci_upper = mean(Value) + qt(0.975, df=n()-1) * sd(Value) / sqrt(n())   
  )


fill_colors <- c("#B8001F", "black", "#507687")


p2 <- 
  ggplot(summary_data, aes(x=Category, y=mean_Value,color=Category)) +
  geom_point(size=5) +      
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), 
                 size=1) +  
  scale_color_manual(values=fill_colors) +
  theme_few()+
  theme(legend.position = "none")+
  ylab("Estimate")+
  xlab(element_blank())

ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/GPP_ate_mean_allcontreat.pdf",p2, width = 3,height = 5, dpi=300)
ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/NDVI_ate_mean_allcontreat.pdf",p2, width = 3,height = 5, dpi=300)
ggsave("/home/j_xiao/data/NCdemo/results/causal_inference/NIRv_ate_mean_allcontreat.pdf",p2, width = 3,height = 5, dpi=300)



