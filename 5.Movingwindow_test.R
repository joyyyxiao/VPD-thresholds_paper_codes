#load data
load("suppl_data/data_monthly.RData")
meco2 <- read.csv("results_mr1/MODIS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- read.csv("results_mr1/GLASS_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- read.csv("results_mr1/FLUXCOM_gpp_ppm_CO2_effect.csv",row.names = 1)
meco2 <- as.matrix(meco2)
meco2 <- meco2[,1:276]
mpar <- mpar [,1:228]
msm <- msm [,1:228]
mtmp <- mtmp [,1:228]
mvpd <- mvpd [,1:228]
mpre <- mpre [,1:228]
#######
zscore <- function(data){
  z_scores <- (data - mean(data, na.rm = TRUE)) / sd(data, na.rm = TRUE)
  return(z_scores)
}
##############

sensi_mon_alpha0 <- function(eco2,tmp,pre,sm,vpd,par,line,win_size){
  
  if(all(is.na(eco2))){
    
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)  
    rownames(temp) <- line
    return(temp)
  }
  
  
  mat <- matrix(data = c(eco2,tmp,pre,sm,vpd,par),nrow=length(eco2),ncol=6,
                dimnames = list(1:length(eco2),c("y","tmp","pre","sm","vpd","par")))
  mat <- na.omit(mat)
  
  if(nrow(mat) < (win_size-20)){
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)
    rownames(temp) <- line
    return(temp)
  }
  x <- mat[,-1]
  y <- mat[,1]
  
  if(sd(y)==0){
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)
    rownames(temp) <- line
    return(temp)
  }
  if (any(apply(x, 2, sd) == 0)) {
    temp <- data.frame(lambda.min=NA, mse=NA, intercept=NA, 
                       coef.tmp=NA, coef.pre=NA, coef.sm=NA, 
                       coef.vpd=NA, coef.par=NA)
    rownames(temp) <- line
    return(temp)
  }
  alpha0.fit <- try(cv.glmnet(x, y, type.measure = "mse", alpha=0, family="gaussian"), silent = TRUE)
  if (inherits(alpha0.fit, "try-error")) {
    temp <- data.frame(lambda.min=NA, mse=NA, intercept=NA, 
                       coef.tmp=NA, coef.pre=NA, coef.sm=NA, 
                       coef.vpd=NA, coef.par=NA)
    rownames(temp) <- line
    return(temp)
  }
  lambda.min <- alpha0.fit$lambda.min
  mse <- alpha0.fit$cvm[[alpha0.fit$index[1]]]
  coef <- coef(alpha0.fit, s = "lambda.min")
  temp <- data.frame(lambda.min=lambda.min,mse=mse,intercept=coef["(Intercept)",],coef.tmp=coef["tmp",],coef.pre=coef["pre",],coef.sm=coef["sm",],coef.vpd=coef["vpd",],coef.par=coef["par",])
  rownames(temp) <- line
  return(temp)
}


##############
## Moving window sensitivity analysis
## Explore different window sizes

library(raster)
library(glmnet)
library(trend)  
library(foreach)
library(parallel)
library(doParallel)

mask <- raster("suppl_data/SYNMAP_biomass_4type_2000.tif")

cl <- makeCluster(4)
#registerDoSNOW(cl)
registerDoParallel(cl)
#pb <- txtProgressBar(max =dim(meco2)[1], style = 3)
#progress <- function(n) setTxtProgressBar(pb, n)
#opts <- list(progress = progress)



window_sizes <- c(60, 84, 120, 180)  # 5, 7, 10, 15years
sens_results <- list()  

for (win_size in window_sizes) {
  mean_sens_values <- c()  

  for (i in seq(1, 228, by = 12)) {  #
    if (i + win_size - 1 > 228) break  
    
    wineco <- meco2[, i:(i + win_size - 1)]
    wintmp <- mtmp[, i:(i + win_size - 1)]
    winpre <- mpre[, i:(i + win_size - 1)]
    winsm  <- msm[, i:(i + win_size - 1)]
    winvpd <- mvpd[, i:(i + win_size - 1)]
    winpar <- mpar[, i:(i + win_size - 1)]
    
    lines <- lonlat  
    nlines <- length(lonlat)
    
    chunk_size <- 5000
    all_idx    <- seq(1, nlines, by=chunk_size)
    
    sensi.df <- foreach(start = all_idx, .combine = rbind, 
      .packages = c("glmnet")
    ) %dopar% {
        end <- min(start + chunk_size - 1, nlines)
        sub_lines <- lines[start:end]
      
        partial_res <- lapply(sub_lines, function(line) {
          sensi_mon_alpha0(
            eco2 = zscore(wineco[line, ]),
            tmp  = zscore(wintmp[line, ]),
            pre  = zscore(winpre[line, ]),
            sm   = zscore(winsm[line, ]),
            vpd  = zscore(winvpd[line, ]),
            par  = zscore(winpar[line, ]),
            line = line,
            win_size = win_size
          )
        })
        do.call(rbind, partial_res)
    }
    
    if(!identical(rownames(sensi.df), lines)){
      sensi.df <- sensi.df[match(lines, rownames(sensi.df)), ]
    }
    
    vpd_sens_matrix <- matrix(sensi.df$coef.vpd, nrow = 360, ncol = 720)
    
    r <- raster(vpd_sens_matrix, 
                xmn = -180, xmx = 180, 
                ymn = -90, ymx = 90, 
                crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
    r <- mask(r,mask)
    r <- abs(r)
    
    
    mean_val <- cellStats(r, "mean", na.rm = TRUE)
    sd_val   <- cellStats(r, "sd", na.rm = TRUE)
    
    lower_limit  <- mean_val - 1.5 * sd_val
    upper_limit  <- mean_val + 1.5 * sd_val
    r_cleaned    <- r
    r_cleaned[r_cleaned < lower_limit | r_cleaned > upper_limit] <- NA
    
    mean_cleaned <- cellStats(r_cleaned, "mean", na.rm = TRUE)
    
    mean_sens_values <- c(mean_sens_values, mean_cleaned)
    
    message(paste0("Window size ", win_size, 
                   ", time index = ", i, 
                   " finished, cleaned_mean = ", mean_cleaned))
  }
  
  if (length(mean_sens_values) > 1) {
    scaled_mean <- (mean_sens_values - min(mean_sens_values, na.rm=TRUE)) /
      (max(mean_sens_values, na.rm=TRUE) - min(mean_sens_values, na.rm=TRUE))
  } else {
    scaled_mean <- mean_sens_values  
  }
  
  sens_results[[as.character(win_size)]] <- data.frame(
    window      = seq_along(mean_sens_values),
    mean        = mean_sens_values,
    scaled_mean = scaled_mean,
    vi          = "FLUXCOM GPP"
  )
}

###stop cluster###
registerDoSEQ()
stopCluster(cl)

sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(sens_results,"results_mr1/MODIS_GPP_sens_globalmean.csv")

sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(sens_results,"results_mr1/GLASS_GPP_sens_globalmean.csv")

sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(sens_results,"results_mr1/FLUXCOM_GPP_sens_globalmean.csv")


##### Mann-Kendall (MK) test
#read data
t_m <- read.csv("results_mr1/MODIS_GPP_sens_globalmean.csv",row.names = 1)
t_g <- read.csv("results_mr1/GLASS_GPP_sens_globalmean.csv",row.names = 1)
t_f <- read.csv("results_mr1/FLUXCOM_GPP_sens_globalmean.csv",row.names = 1)
sens_results <- split(t_m, t_m$window_length)
sens_results <- split(t_g, t_g$window_length)
sens_results <- split(t_f, t_f$window_length)

global_mean_sens_mk <- data.frame(VI=rep(c("MODIS GPP","GLASS GPP","FLUXCOM GPP"),each=4),
                                  win_size=rep(c(5,7,10,15)*12,3),
                                  z=rep(NA,12),
                                  p=rep(NA,12))



for (win_size in names(sens_results)) {
  df <- sens_results[[win_size]]
  
  mk_test <- trend::mk.test(df$scaled_mean)  
  
  
  index <- which(global_mean_sens_mk$VI == "MODIS GPP" & global_mean_sens_mk$win_size == as.numeric(win_size))
  

  global_mean_sens_mk$z[index] <- mk_test$statistic
  global_mean_sens_mk$p[index] <- mk_test$p.value
  
  print(paste0("Window size: ", win_size, 
               " | MK test Z-value: ", mk_test$statistic, 
               " | p-value: ", mk_test$p.value))
}

write.csv(global_mean_sens_mk,"results_mr1/global_mean_vpdsens_acrosspds_mk.csv")

##plot
library(ggplot2)
library(ggthemes)
library(dplyr)
global_mean_sens_mk <- global_mean_sens_mk %>%
  mutate(z = round(z, 2))  

#Fig. S6
p <- ggplot(global_mean_sens_mk, aes(x = factor(win_size), y = VI, fill = p)) +
  geom_tile(color = "white") +  
  geom_text(aes(label = z), color = "white", size = 5) +  
  scale_fill_gradient(low = "lightblue", high = "darkblue", limits= c(0,0.05),name = "p value") +  
  scale_x_discrete(labels = c("5", "7", "10", "15")) +
  theme_minimal() +
  labs(
    x = "Length of moving window (years)",
    y = "Products",
    title = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("Trends_movingwindow_test.pdf",p, width = 6,height = 4, dpi=300)