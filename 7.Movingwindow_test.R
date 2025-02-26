##############

sensi_mon_alpha0 <- function(eco2,tmp,pre,sm,vpd,par,line,win_size){
  
  if(all(is.na(eco2))){
    
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)  
    rownames(temp) <- line
    return(temp)
  }
  
  
  mat <- matrix(data = c(eco2,zscore(tmp),zscore(pre),zscore(sm),zscore(vpd),zscore(par)),nrow=length(eco2),ncol=6,
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
  
  alpha0.fit <- cv.glmnet(x,y,type.measure = "mse",alpha=0,family="gaussian")  
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
library(trend)  
mask <- shapefile("/home/j_xiao/data/Terrestrial Ecoregions of the World_WWF/mask_shpfile.shp")



library(foreach)
#library(doSNOW)
library(parallel)
library(doParallel)

cl <- makeCluster(8)
#registerDoSNOW(cl)
registerDoParallel(cl)
#pb <- txtProgressBar(max =dim(meco2)[1], style = 3)
#progress <- function(n) setTxtProgressBar(pb, n)
#opts <- list(progress = progress)



window_sizes <- c(60, 84, 120, 180)  # 5, 7, 10, 15years
sens_results <- list()  

for (win_size in window_sizes) {
  mean_sens_values <- c()  

  for (i in seq(1, 276, by = 12)) {  #
    if (i + win_size - 1 > 276) break  
    
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
    r <- mask(r,mask,inverse=T)
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
    vi          = "NDVI"
  )
}

###stop cluster###
stopCluster(cl)

GPP_sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(GPP_sens_results,"/home/j_xiao/data/NCdemo/results/FLUXCOM_GPP/GPP_sens_globalmean.csv")

NIRv_sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(NIRv_sens_results,"/home/j_xiao/data/NCdemo/results/MODIS_NIRv/NIRv_sens_globalmean.csv")

NDVI_sens_results <- dplyr::bind_rows(sens_results, .id = "window_length")
write.csv(NDVI_sens_results,"/home/j_xiao/data/NCdemo/results/GIMMS_NDVI/NDVI_sens_globalmean.csv")





##### Mann-Kendall (MK) test
#read data
t1 <- read.csv("/home/j_xiao/data/NCdemo/results/FLUXCOM_GPP/GPP_sens_globalmean.csv")[,-1]
t2 <- read.csv("/home/j_xiao/data/NCdemo/results/MODIS_NIRv/NIRv_sens_globalmean.csv")[,-1]
t3 <- read.csv("/home/j_xiao/data/NCdemo/results/GIMMS_NDVI/NDVI_sens_globalmean.csv")[,-1]
sens_results <- split(t1, t1$window_length)
sens_results <- split(t2, t2$window_length)
sens_results <- split(t3, t3$window_length)

global_mean_sens_mk <- data.frame(VI=rep(c("GPP","NIRv","NDVI"),each=4),
                                  win_size=rep(c(5,7,10,15)*12,3),
                                  z=rep(NA,12),
                                  p=rep(NA,12)
)


for (win_size in names(sens_results)) {
  df <- sens_results[[win_size]]
  
  mk_test <- trend::mk.test(df$scaled_mean)  
  
  
  index <- which(global_mean_sens_mk$VI == "NDVI" & global_mean_sens_mk$win_size == as.numeric(win_size))
  

  global_mean_sens_mk$z[index] <- mk_test$statistic
  global_mean_sens_mk$p[index] <- mk_test$p.value
  
  print(paste0("Window size: ", win_size, 
               " | MK test Z-value: ", mk_test$statistic, 
               " | p-value: ", mk_test$p.value))
}

write.csv(global_mean_sens_mk,"./results/global_mean_vpdsens_acrossVIs_mk.csv")

##plot
global_mean_sens_mk <- global_mean_sens_mk %>%
  mutate(z = round(z, 2))  

ggplot(global_mean_sens_mk, aes(x = factor(win_size), y = VI, fill = p)) +
  geom_tile(color = "white") +  
  geom_text(aes(label = z), color = "white", size = 5) +  
  scale_fill_gradient(low = "lightblue", high = "darkblue", limits= c(0,0.001),name = "p value") +  # 颜色梯度
  theme_minimal() +
  labs(
    x = "Length of moving window (y)",
    y = "Variable (VI)",
    title = "P values of the declining trends in VPD sensitivity"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )
