
library(xts)
library(zoo)
library(lubridate)
library(TSS.RESTREND)
library(jsonlite)

library(mblm)
library(raster)

setwd("/home/j_xiao/data/NCdemo")

VIdf <- read.csv("./data/demo_dataframe_NIRv.csv",row.names=1, check.names=FALSE)
fnC4 <- "./data/demo_dataframe_C4frac.csv"
C4df <- read.csv(fnC4, row.names=1, check.names=FALSE) # C4 frac
C4df[C4df < 0] = 0.
C4df[C4df > 1] = 1.
C4frac <- C4df

CO2ts <- read.csv("./data/CO2_forts.csv") # CO2 time series
CO2ts <-  ts(CO2ts[23:779,5], start=c(1960, 1), end=c(2023,1), frequency = 12) #monthly 4=average 5=deseasonal

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
  # mst <- start(CTSR.VI)[2]
  yfn <- start(tail(CTSR.VI, n=1))[1]
  
  # # ========== get the CO2 concentration of the dataset years ==========
  # +++++ build a matching CO2 series +++++
  Ca <- as.numeric(window(CO2, yst, yfn+1))[2:(length(as.numeric(window(CO2, yst, yfn+1)))-1)]
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
    results = rep(NA,2)
    
  }else{
    if (!par){print(line)}
    
    # ========== Deal with Dates ==========
    # +++++ COnvert the VI to a TS object +++++
    # VIdates <- as.POSIXlt(colnames(VI))
    # VIys    <- VIdates[1]$year + 1900
    # VIms    <- month(VIdates[1])
    # VIyf    <- tail(VIdates, n=1)$year+ 1900
    # VImf    <- month(tail(VIdates, n=1))
    CTSR.VI <- ts(as.numeric(VI), start=c(2000, 2), end=c(2022,12), frequency = 12)
    
    # ========== Get the annual Max VI values ==========
    rawmax.df  <- AnMaxVI(CTSR.VI)
    CTSR.VIadj = franksCO2m(
      CTSR.VI, CO2= CO2ts, C4frac=round(C4frac, digits = 4), refyear=1980)
    adjmax.df  <- AnMaxVI(CTSR.VIadj)
    
    # ========== Calculate Observed Change ==========
    ti       <- time(rawmax.df$Max)
    obsVI    <- rawmax.df$Max
    OBStheil <- mblm(obsVI~ti, repeated=FALSE)
    ObservedChange = as.numeric(OBStheil[[1]][2]) *23
    
    # ========== Calculate CO2 Change ==========
    if (C4frac < 1){
      
      CO2change <- rawmax.df$Max - adjmax.df$Max
      CO2theil  <- mblm(CO2change~ti, repeated=FALSE)
      CO2 = as.numeric(CO2theil[[1]][2]) *23
      
    }else{
      # Avoid cases where the CO2 values will all be zero
      CO2 = 0.0
    }
    
    
    
    results = c(ObservedChange,CO2)
  }
  # ========== return the results ==========
  ret <- t(results)
  rownames(ret) <- line         # add the row name back
  colnames(ret) <- c("ObservedChange","TotalCO2change")
  return(ret)
}
# eCO2 monthly effects function
tssr.attr.monthly.ppm <- function(line, VI, C4frac, CO2,par){
  # =========== Function is applied to one pixel at a time ===========
  # ========== Perfrom the data checks and if anything fails skipp processing ==========
  # There is a data check for NANs in the TSSRattribution function, If SkipError is True
  # It then returns an opject of the same structure as actual results but filled with NaN
  # Usefull stacking using the foreac::do command.
  if (all(is.na(VI))){
    results = rep(NA,275)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- colnames(VIdf)
    return(ret)
  }
  if (!par){print(line)}
  # ========== Deal with Dates ==========
  # +++++ COnvert the VI to a TS object +++++
  # VIdates <- as.POSIXlt(colnames(VI))
  # VIys    <- VIdates[1]$year + 1900
  # VIms    <- month(VIdates[1])
  # VIyf    <- tail(VIdates, n=1)$year+ 1900
  # VImf    <- month(tail(VIdates, n=1))
  CTSR.VI <- ts(as.numeric(VI), start=c(2000, 2), end=c(2022,12), frequency = 12)
  if (sum(is.na(CTSR.VI))>75){
    results = rep(NA,275)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- colnames(VIdf)
    return(ret)
  }
  if (sd(VI,na.rm=T)==0){
    results = rep(NA,275)
    ret <- t(results)
    rownames(ret) <- line         # add the row name back
    colnames(ret) <- colnames(VIdf)
    return(ret)
  }
  CTSR.VI <- na.interp(CTSR.VI)
  # ========== CO2_dif calculate ==========
  CO2_dif <- function(CO2,start,end,refyear=1980){
    if (class(CO2) != "ts"){
      stop("CTSR.VI Not a time series object. Please check the data")
    }
    # ========== Calculate the baseline ==========
    baseline_Ca = window(CO2, refyear, refyear+1)[[1]]
    # # ========== get the CO2 concentration of the dataset years ==========
    # +++++ build a matching CO2 series +++++
    Ca <- as.numeric(window(CO2, start, end+1))[2:(length(as.numeric(window(CO2, start, end+1)))-1)]
    dif <- Ca-baseline_Ca
    return(dif)
  }
  # ========== get the results ==========
  CTSR.VIadj = franksCO2m(
    CTSR.VI, C4frac=round(C4frac, digits = 4),
    CO2=CO2ts, refyear=1980)
  co2dif <- CO2_dif(CO2=CO2ts,start=2000,end=2022,refyear=1980)
  CO2change <- (CTSR.VI - CTSR.VIadj)/co2dif
  results <- as.vector(CO2change)
  if (is.null(results)){
    browser()
  }
  # ========== return the results ==========
  ret <- t(results)
  rownames(ret) <- line         # add the row name back
  colnames(ret) <- colnames(VIdf)
  return(ret)
}


# ========== Loop over the rows in parallel ==========
library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(30)
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
      tssr.attr.monthly.ppm(line, VIdf[line, ],C4df[line, ], par=TRUE)
    })
print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)

###save data###
write.csv(tss.atdf,"./results/NIRv_monthly_ppm_CO2_effect.csv")

