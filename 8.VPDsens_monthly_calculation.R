library(glmnet)
load("./data/demo_data_monthly.RData")
x <- matrix(data = c(zscore(mtmp[line,]),zscore(mpre[line,]),zscore(mvpd[line,]),zscore(mpar[line,]),zscore(msm[line,])),nrow=276,ncol=5,
            dimnames = list(1:276,c("tmp","pre","vpd","par","sm")))

##split data

train_rows <- sample(1:length(y),0.8*10)
x.train <- x[train_rows,]
x.test <- x[-train_rows,] 

y.train <- y[train_rows]
y.test <- y[-train_rows] 

alpha0.fit <- cv.glmnet(x.train,y.train,type.measure = "mse",alpha=0,family="gaussian")  #set the type.measure to "deviance" when applying Elastic-Net Regression to Logistic Regression; alpha=0 --Ridge, alpha=1, Lasso, ortherwise mixture; family="gaussian" indicates we are doing linear regression, set it to "binomial" when we doing logistic regression
alpha0.predicted <- predict(alpha0.fit,s=alpha0.fit$lambda.min, newx=x.test) #s indicates the size of the penalty; lambda.1se -- within 1 se of the lambda.min & simplest model 
plot(alpha0.fit) 
plot(alpha0.predicted,y.test)
abline(0,1)
coef(alpha0.fit, s = "lambda.min")


##############
sensi_mon_alpha0 <- function(eco2,tmp,pre,sm,vpd,par,line){

  if(all(is.na(eco2))){
    
    temp <- data.frame(lambda.min=NA,mse=NA,intercept=NA,coef.tmp=NA,coef.pre=NA,coef.sm=NA,coef.vpd=NA,coef.par=NA)  
    rownames(temp) <- line
    return(temp)
  }
  
  
      mat <- matrix(data = c(eco2,zscore(tmp),zscore(pre),zscore(sm),zscore(vpd),zscore(par)),nrow=length(eco2),ncol=6,
                  dimnames = list(1:length(eco2),c("y","tmp","pre","sm","vpd","par")))
      mat <- na.omit(mat)
      
  if(nrow(mat)<100){
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

sensi_mon_alpha0(eco2=zscore(meco2[line,]),tmp=mtmp[line,],pre=mpre[line,],sm=msm[line,],vpd=mvpd[line,],par=mpar[line,],line)

################
library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(30)
registerDoSNOW(cl)

# Setup the progress bar


pb <- txtProgressBar(max =dim(meco2)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

print(paste("Starting parellel processing at:", Sys.time()))
ptime  <- system.time(
  sensi.df <- foreach(
    line=rownames(meco2), .combine = rbind, .options.snow = opts,
    .packages = c('glmnet')) %dopar% {
      
      sensi_mon_alpha0(eco2=zscore(meco2[line,]),tmp=mtmp[line,],pre=mpre[line,],sm=msm[line,],vpd=mvpd[line,],par=mpar[line,],line)
      
      
    })
print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)

###save data###
write.csv(sensi.df,"./results/FLUXCOM_GPP/FLUXCOM_GPP_Ridge_sensitivity_alpha0_ERA5_monthly.csv")


###make raster###
library(raster)
vpd_sens_map <- as_raster(sensi.df$coef.vpd)
plot(vpd_sens_map,col=rainbow(20))
mask <- shapefile("/home/j_xiao/data/Terrestrial Ecoregions of the World_WWF/mask_shpfile.shp")
vpd_sens_map <- mask(vpd_sens_map,mask,inverse=T)
writeRaster(vpd_sens_map,"./results/FLUXCOM_GPP/FLUXCOM_GPP_vpd_monthly_sensmap_alpha0_ERA5.tif")

############

##Moving window
#window == 10 years ==120 months

for(i in seq(1,228,by=12)[1:10]){
  
  wineco2 <- meco2[,i:(i+119)]
  wintmp <- mtmp[,i:(i+119)]
  winpre <- mpre[,i:(i+119)]
  winsm <- msm[,i:(i+119)]
  winvpd <- mvpd[,i:(i+119)]
  winpar <- mpar[,i:(i+119)]
  
  sensi.df <- foreach(
    line=rownames(wineco2), .combine = rbind, .options.snow = opts,
    .packages = c('glmnet')) %dopar% {
      sensi_mon_alpha0(eco2=zscore(wineco2[line,]),tmp=wintmp[line,],pre=winpre[line,],sm=winsm[line,],vpd=winvpd[line,],par=winpar[line,],line)
      
    }
  
  vpd_sens_matrix <- matrix(sensi.df$coef.vpd,nrow=360,ncol=720)
  vpd_sens_map <- raster(vpd_sens_matrix, xmn=-180,xmx=180,ymn=-90,ymx=90,crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  writeRaster(vpd_sens_map,paste0("./results/FLUXCOM_GPP/FLUXCOM_GPP_vpd_sensmap_alpha0_ERA5_monthly","_window",i,".tif"),overwrite=T)
  print(paste0("window",i,"finished"))           
  
}


registerDoSEQ()
stopCluster(cl)



