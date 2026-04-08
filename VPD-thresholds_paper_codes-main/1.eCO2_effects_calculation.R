
library(xts)
library(zoo)
library(lubridate)
library(TSS.RESTREND)
library(jsonlite)
library(forecast)
library(mblm)
library(raster)


VIdf <- read.csv("GPPdatasets/FLUXCOM_GPP_monthly_dataframe.csv",row.names=1, check.names=FALSE)
VIdf <- read.csv("GPPdatasets/GLASS_GPP_monthly_dataframe.csv",row.names=1, check.names=FALSE)
VIdf <- read.csv("GPPdatasets/MODIS_GPP_monthly_dataframe.csv",row.names=1, check.names=FALSE)


C4df <- read.csv("suppl_data/dataframe_C4frac.csv", row.names=1, check.names=FALSE) # C4 frac
C4df[C4df < 0] = 0.
C4df[C4df > 1] = 1.


CO2ts <- read.csv("suppl_data/co2_mm_mlo.csv") # CO2 time series
CO2ts <-  ts(CO2ts[,4], start=c(1958, 3), end=c(2026,2), frequency = 12) #monthly col4=average 

# ==========frank CO2 function==========
franksCO2m <- function(CTSR.VI, C4frac, CO2, refyear=1980){
  # ========== Check the data ==========
  if (class(CO2) != "ts"){
    stop("CTSR.VI Not a time series object. Please check the data")
  }else if ((C4frac > 1)||(C4frac<0)){
    print(C4frac)
    stop("C4frac must be between 0 and 1")
  }
  
  
  franks_FvC <- function(Ca){
    theta = 0.7         # shape of the light response curve
    gamma_star = 40.0   # CO2 compensation point
    return ((theta * Ca - gamma_star) / (theta * Ca + 2.0 * gamma_star))
  }
  
  # ========== Calculate the baseline ==========
  baseline_Ca = window(CO2, refyear, refyear+1)[[1]]
  baseAnet    = franks_FvC(baseline_Ca)
  
  # +++++ Get the start year and start month of the data +++++
  yst <- start(CTSR.VI)[1]
  mst <- start(CTSR.VI)[2]
  yfn <- start(tail(CTSR.VI, n=1))[1]
  mfn <- start(tail(CTSR.VI, n=1))[2]
  
  # # ========== get the CO2 concentration of the dataset years ==========
  # +++++ build a matching CO2 series +++++
  Ca <- as.numeric(window(CO2, c(yst,mst), c(yfn,mfn)))
  # +++++ Check the length +++++
  if (length(CTSR.VI) == length(Ca)){
    Ca <- Ca
  }else{
    stop("Lengh problem in the CO2 adjustment, Ts must be either annual or monthly")
  }
  
  Anet           = franks_FvC(Ca)
  model_response = Anet / baseAnet
  
  # ========== Adjust the CTSRVI ==========
  VIadjC3 <- CTSR.VI / model_response
  
  # ========== Account for C4  ==========
  CTSR.VIadj <- (VIadjC3 * (1-C4frac)) + (CTSR.VI *C4frac)
  
  return(CTSR.VIadj)
  
}

# eCO2 total effects function
tssr.attr.total <- function(line, VI, C4frac, par){
  # =========== Function is applied to one pixel at a time ===========
  # ========== Perfrom the data checks and if anything fails skipp processing ==========
  # There is a data check for NANs in the TSSRattribution function, If SkipError is True
  # It then returns an opject of the same structure as actual results but filled with NaN
  # Usefull stacking using the foreac::do command.
  if (any(is.na(VI)) | (sd(VI)==0)){
    results = rep(NA,3)
    
  }else{
    if (!par){print(line)}
    
    # ========== Deal with Dates ==========
    # +++++ COnvert the VI to a TS object +++++

#     CTSR.VI <- ts(as.numeric(VI), start=c(1980, 1), end=c(2018,12), frequency = 12) #FLUXCOMGPP
#     CTSR.VI <- ts(as.numeric(VI), start=c(2000, 1), end=c(2018,12), frequency = 12) #GLASSGPP
      CTSR.VI <- ts(as.numeric(VI), start=c(2000, 1), end=c(2025,12), frequency = 12) #MODISGPP
     
     CTSR.VI <- window(CTSR.VI, start=c(2000, 1), end=c(2022,12))
    # ========== Get the annual  VI values ==========
    AnSum <- function(x, na.rm = TRUE) {
       
       if (!is.ts(x)) stop("Input must be a ts object.")
       if (frequency(x) != 12) stop("Must be monthly (frequency = 12).")
       
       yrs <- floor(time(x))
       mos <- cycle(x)
       vals <- as.numeric(x)
       
       mdays_common <- c(31,28,31,30,31,30,31,31,30,31,30,31)
       mdays <- mdays_common[mos]
       
       
       annual_sum <- tapply(vals * mdays, yrs, sum, na.rm = na.rm)
       
       out <- data.frame(
         Annual.sum = as.numeric(annual_sum),
         Annual.mean.rate = as.numeric(annual_sum) / 365
       )
       
       rownames(out) <- names(annual_sum)
       
       return(out)
     }
    rawmax.df  <- AnSum(CTSR.VI)
    CTSR.VIadj = franksCO2m(
      CTSR.VI, CO2= CO2ts, C4frac=round(C4frac, digits = 4), refyear=1980)
    adjmax.df  <- AnSum(CTSR.VIadj)
    
    # ========== Calculate Observed Change ==========
    ti       <- time(rawmax.df$Annual.sum)
    t0 <- min(ti)
    t1 <- max(ti)
    span <- t1 - t0
    obsVI    <- rawmax.df$Annual.sum
    OBStheil <- mblm(obsVI~ti, repeated=FALSE)
    ObservedChange = as.numeric(OBStheil[[1]][2]) *span
    
    # ========== Calculate CO2-induced Change ==========
    if (C4frac < 1){
      
      CO2change <- rawmax.df$Annual.sum - adjmax.df$Annual.sum
      CO2theil  <- mblm(CO2change~ti, repeated=FALSE)
      AbsoluteChange = as.numeric(CO2theil[[1]][2]) *span
      
      b0 <- coef(CO2theil)[1]
      b1 <- coef(CO2theil)[2]
      E1_hat  <- b0 + b1 * t0
      RelativeChange = AbsoluteChange/E1_hat 
      
    }else{
      # Avoid cases where the CO2 values will all be zero
      AbsoluteChange = 0.0
      RelativeChange = 0.0
    }

    results = c(ObservedChange,AbsoluteChange,RelativeChange)
  }
  # ========== return the results ==========
  ret <- t(results)
  rownames(ret) <- line         # add the row name back
  colnames(ret) <- c("ObservedChange","AbsoluteChange","RelativeChange")
  return(ret)
}

# eCO2 monthly effects function
# 
tssr.attr.monthly.ppm <- function(line, VI, C4frac,par){
  # =========== Function is applied to one pixel at a time ===========
  # ========== Perfrom the data checks and if anything fails skipp processing ==========
  # There is a data check for NANs in the TSSRattribution function, If SkipError is True
  # It then returns an opject of the same structure as actual results but filled with NaN
  # Usefull stacking using the foreac::do command.
  year <- 26
  if (all(is.na(VI))){
    results = rep(NA,year*12)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- 1:(year*12)
    return(ret)
  }
  if (!par){print(line)}
  # ========== Deal with Dates ==========
  # +++++ COnvert the VI to a TS object +++++

  # CTSR.VI <- ts(as.numeric(VI), start=c(1980, 1), end=c(2018,12), frequency = 12) #FLUXCOMGPP
  # CTSR.VI <- ts(as.numeric(VI), start=c(2000, 1), end=c(2018,12), frequency = 12) #GLASSGPP
    CTSR.VI <- ts(as.numeric(VI), start=c(2000, 1), end=c(2025,12), frequency = 12) #MODISGPP

    CTSR.VI <- window(CTSR.VI, start=c(2000, 1), end=c(2022,12))
  
  MonTotal <- function(x) {
    if (!is.ts(x)) stop("Input must be a ts object.")
    if (frequency(x) != 12) stop("Must be monthly (frequency = 12).")
    
    mos  <- cycle(x)
    vals <- as.numeric(x)
    
    mdays_common <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    mdays <- mdays_common[mos]
    
    out_vals <- vals * mdays
    
    out <- ts(out_vals, start = start(x), frequency = 12)
    return(out)
  }
  CTSR.VI <- MonTotal(CTSR.VI)
  
  if (sum(is.na(CTSR.VI))>75){
    results = rep(NA,year*12)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- 1:(year*12)
    return(ret)
  }
  if (sd(VI,na.rm=T)==0){
    results = rep(NA,year*12)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- 1:(year*12)
    return(ret)
  }
  CTSR.VI <- forecast::na.interp(CTSR.VI)
  # ========== CO2_dif calculate ==========
  CO2_dif <- function(CO2,start,end,refyear=1980){
    if (class(CO2) != "ts"){
      stop("CTSR.VI Not a time series object. Please check the data")
    }
    # ========== Calculate the baseline ==========
    baseline_Ca = window(CO2, refyear, refyear+1)[[1]]
    # # ========== get the CO2 concentration of the dataset years ==========
    # +++++ build a matching CO2 series +++++
    Ca <- as.numeric(window(CO2, start, end+1))[1:(length(as.numeric(window(CO2, start, end+1)))-1)]
    dif <- Ca-baseline_Ca
    return(dif)
  }
  # ========== get the results ==========
  CTSR.VIadj = franksCO2m(
    CTSR.VI, C4frac=round(C4frac, digits = 4),
    CO2=CO2ts, refyear=1980)
  co2dif <- CO2_dif(CO2=CO2ts,start=2000,end=2018,refyear=1980)
  CO2change <- (CTSR.VI - CTSR.VIadj)/co2dif
  results <- as.vector(CO2change)
  if (is.null(results)){
    browser()
  }
  # ========== return the results ==========
  ret <- t(results)
  rownames(ret) <- line         # add the row name back
  colnames(ret) <- 1:(year*12)
  return(ret)
}


# ========== Loop over the rows in parallel ==========
library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(4)
registerDoSNOW(cl)

# Setup the progress bar
print(paste("Starting parellel processing at:", Sys.time()))

pb <- txtProgressBar(max =dim(VIdf)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

ptime  <- system.time(
  tss.atdf <- foreach(
    line=rownames(VIdf), .combine = rbind, .options.snow = opts,
    .packages = c('TSS.RESTREND', "xts", "lubridate","mblm")) %dopar% {
      tssr.attr.total(line, VIdf[line, ],C4df[line, ], par=TRUE)
    })                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)

###save data###
write.csv(tss.atdf,"results_mr1/GLASS_gpp_total_CO2_effect.csv")
write.csv(tss.atdf,"results_mr1/FLUXCOM_gpp_total_CO2_effect.csv")
write.csv(tss.atdf,"results_mr1/MODIS_gpp_total_CO2_effect.csv")
#tss.atdf <- read.csv("results_mr1/FLUXCOM_gpp_total_CO2_effect.csv")
#tss.atdf <- read.csv("results_mr1/GLASS_gpp_total_CO2_effect.csv")
#tss.atdf <- read.csv("results_mr1/MODIS_gpp_total_CO2_effect.csv")

biome <- raster("suppl_data/SYNMAP_biomass_4type_2000.tif")
r1 <- as_raster(tss.atdf$AbsoluteChange)
r1 <- mask(r1, biome)
r2 <- as_raster(tss.atdf$RelativeChange)
r2 <- mask(r2, biome)
writeRaster(r1,"results_mr1/fluxcom_abtotalchange_map.tif", overwrite = T)
writeRaster(r2,"results_mr1/fluxcom_retotalchange_map.tif", overwrite = T)

writeRaster(r1,"results_mr1/glass_abtotalchange_map.tif", overwrite = T)
writeRaster(r2,"results_mr1/glass_retotalchange_map.tif", overwrite = T)

writeRaster(r1,"results_mr1/modis_abtotalchange_map.tif", overwrite = T)
writeRaster(r2,"results_mr1/modis_retotalchange_map.tif", overwrite = T)


########monthly change per ppm
write.csv(tss.atdf,"results_mr1/FLUXCOM_gpp_ppm_CO2_effect.csv")
write.csv(tss.atdf,"results_mr1/GLASS_gpp_ppm_CO2_effect.csv")
write.csv(tss.atdf,"results_mr1/MODIS_gpp_ppm_CO2_effect.csv")

