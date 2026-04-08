library(library)
library(dplyr)
library(raster)

as_array <- function (raster, lonlat = NULL, nrows = 360, ncols = 720) {
  if (!is.null(lonlat)) {
    lonlat <- read.csv(lonlat)[, 2]
  }
  else {
    lonlat <- 1:(nrows * ncols)
  }
  matrix <- matrix(NA, nrows * ncols, nlayers(raster))
  for (i in 1:nlayers(raster)) {
    matrix[, i] <- as.numeric(as.matrix(raster[[i]]))
  }
  rownames(matrix) <- lonlat
  colnames(matrix) <- rep(NA, nlayers(raster))
  return(matrix)
}
as_matrix <- 
#GLASS GPP 2000-2018 8-day average gC m-2 day-1 to annual gC m-2 year-1
GLASS_GPP_ann_df <- matrix(
  NA,
  nrow = 360*720,
  ncol = 1
)

for (year in 2000:2018){
  filelist <-  list.files("/Volumes/JoyXiao/j_xiao/data/GLASS_AVHRR_GPP/",full.names = T,pattern = paste0("A",year, ".*\\.tif$"))
  if (length(filelist) != 46) {
    stop("Number of files does not match")
  } #original GLASS AVHRR GPP data downloaded from http://www.glass.umd.edu/
  
  stack <- stack(filelist)
  days <- c(rep(8, 45), 5)
  ann_GLASSGPP <- calc(stack,function(x){  x[x == 65535] <- NA
                                           sum((x*0.01) * days, na.rm = TRUE)})#scale =0.01
  values_year <- as.numeric(as.matrix(ann_GLASSGPP))
  GLASS_GPP_ann_df <- cbind(GLASS_GPP_ann_df, values_year)
  
  print(year)
}
GLASS_GPP_ann_df <- GLASS_GPP_ann_df [,-1]
colnames(GLASS_GPP_ann_df) <- 2000:2018
rownames(GLASS_GPP_ann_df) <- read.csv("line.csv")$x
write.csv(GLASS_GPP_ann_df,"GPPdatasets/GLASS_GPP_annual_dataframe.csv")

#GLASS GPP 2000-2018 8-day average gC m-2 day-1 to monthly gC m-2 day-1
stackall <- list()
for (year in 2000:2018){
  filelist <-  list.files("/Volumes/JoyXiao/j_xiao/data/GLASS_AVHRR_GPP/",full.names = T,pattern = paste0("A",year, ".*\\.tif$"))
  if (length(filelist) != 46) {
    stop("Number of files does not match")
  }#original GLASS AVHRR GPP data downloaded from http://www.glass.umd.edu/
  
  stack_8d <- stack(filelist)
  stack_8d <- calc(stack_8d,function(x){x[x == 65535] <- NA
                                  x <- x*0.01; return(x)})#scale =0.01
  
  days <- c(rep(8, 45), 5)
  stack_365 <- stack(
    lapply(1:46, function(i) {
      replicate(days[i], stack_8d[[i]], simplify = FALSE)
    }) |> unlist(recursive = FALSE)
  )
  
  days_in_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  month_index <- rep(1:12, times = days_in_month)
  
  stack_month_mean <- stack(
    lapply(1:12, function(m) {
      calc(stack_365[[which(month_index == m)]], mean)
    })
  )
  
  stackall[[as.character(year)]] <- stack_month_mean
  print(year)
}
stack_mm <- stack(stackall)
dataframe_GLASSGPP <- as_array(stack_mm,"line.csv")
write.csv(dataframe_GLASSGPP,"GPPdatasets/GLASS_GPP_monthly_dataframe.csv")

#FLUXCOM GPP 2000-2018 monthly average gC m-2 day-1 to annual gC m-2 year-1
VIdf <- read.csv("GPPdatasets/FLUXCOM_GPP_monthly_dataframe.csv",row.names=1, check.names=FALSE) #1980-2018 #FLUXCOM GPP data downloaded from https://www.fluxcom.org/CF-Products/
monthtoannual_GPP <- function(GPP_mat, start_year, end_year, return_annual = FALSE) {
  nyear <- end_year - start_year + 1
  
  if (ncol(GPP_mat) != nyear * 12) {
    stop("Number of columns does not match nyear * 12")
  }
  
  
  years  <- rep(seq(start_year, end_year), each = 12)
  months <- rep(1:12, times = nyear)
  
  
  first_next_month <- as.Date(paste(years, months, "01", sep="-")) + 32
  first_next_month <- as.Date(format(first_next_month, "%Y-%m-01"))
  first_this_month <- as.Date(paste(years, months, "01", sep="-"))
  days_in_month <- as.integer(first_next_month - first_this_month)
  
  # gC/day * days -> gC/month
  GPP_monthly_total <- sweep(GPP_mat, 2, days_in_month, `*`)
  
  
  year_index <- rep(seq_len(nyear), each = 12)
  
  # gC/yr
  GPP_annual <- sapply(seq_len(nyear), function(i) {
    cols <- which(year_index == i)
    rowSums(GPP_monthly_total[, cols, drop = FALSE], na.rm = TRUE)
  })
  # dim:n_pixel × nyear
  
  # mean (gC/yr)
  GPP_mean_annual <- rowMeans(GPP_annual, na.rm = TRUE)
  
  if (return_annual) {
    return(list(mean_annual = GPP_mean_annual, annual = GPP_annual))
  } else {
    return(GPP_mean_annual)
  }
} #function for converting to annual data
FLUXCOM_GPP_ann_df <- as.data.frame(monthtoannual_GPP(VIdf,1980,2018,T)$annual)

colnames(FLUXCOM_GPP_ann_df) <- c(1980:2018)
rownames(FLUXCOM_GPP_ann_df) <- read.csv("line.csv")$x

FLUXCOM_GPP_ann_df <- FLUXCOM_GPP_ann_df[,as.character(2000:2018)]
write.csv(FLUXCOM_GPP_ann_df,"GPPdatasets/FLUXCOM_GPP_annual_dataframe.csv")

#MODIS GPP 2000-2022 monthly and annual
list.files("/Volumes/JoyXiao/j_xiao/data/MODIS_GPP_05deg") #original MODIS GPP data downloaded from https://lpdaac.usgs.gov
filelist <- paste0("/Volumes/JoyXiao/j_xiao/data/MODIS_GPP_05deg/MODIS_GPP_resample",rep(2000:2025,each=12),"_",c(1:12),".tif")
stack <- stack(filelist)
VIdf <- as_array(stack,"line.csv")
write.csv(VIdf,"GPP_datasets/MODIS_GPP_monthly_dataframe.csv")

calc_sum_annual <- function(data,start_year,end_year) {
  
  nyear <- end_year - start_year + 1
  
  if (ncol(data) != nyear * 12) {
    stop("Number of columns does not match nyear * 12")
  }
  
  year_index <- rep(seq_len(nyear), each = 12)
  days_in_month <- rep(c(31,28,31,30,31,30,31,31,30,31,30,31),nyear)
  
  annual <- sapply(seq_len(nyear), function(i) {
    cols <- which(year_index == i)
    days <- days_in_month[cols]
    apply(data[, cols, drop = FALSE], 1, function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        sum(x * days, na.rm = F)
      }
    })
  })
}

GPP3_ann_df <- calc_sum_annual(VIdf,2000,2025)
colnames(GPP3_ann_df) <- 2000:2025
write.csv(GPP3_ann_df,"GPPdatasets/MODIS_GPP_annual_dataframe.csv")

