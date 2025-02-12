library(chngpt)
library(h2o)
h2o.init()
h2o.no_progress()
h2o.removeAll()
h2o.ls()
line <- rownames(eco2_m)
which(line== "(-166.75, 66.25)")
#h2o.shutdown()


#######
h2o.init(nthreads = -1)
makeDF <- function(pd = pd, varimp = varimp, threshold = th, modelid = modelid){
  df <- data.frame(pd.vpd1=pd[1,"vpd"],
                   pd.vpd2=pd[2,"vpd"],
                   pd.vpd3=pd[3,"vpd"], 
                   pd.vpd4=pd[4,"vpd"], 
                   pd.vpd5=pd[5,"vpd"], 
                   pd.vpd6=pd[6,"vpd"], 
                   pd.vpd7=pd[7,"vpd"], 
                   pd.vpd8=pd[8,"vpd"], 
                   pd.vpd9=pd[9,"vpd"], 
                   pd.vpd10=pd[10,"vpd"], 
                   pd.vpd11=pd[11,"vpd"], 
                   pd.vpd12=pd[12,"vpd"], 
                   pd.vpd13=pd[13,"vpd"], 
                   pd.vpd14=pd[14,"vpd"], 
                   pd.vpd15=pd[15,"vpd"], 
                   pd.vpd16=pd[16,"vpd"], 
                   pd.vpd17=pd[17,"vpd"], 
                   pd.vpd18=pd[18,"vpd"], 
                   pd.vpd19=pd[19,"vpd"], 
                   pd.vpd20=pd[20,"vpd"], 
                   pd.mean_response1=pd[1,"mean_response"], 
                   pd.mean_response2=pd[2,"mean_response"], 
                   pd.mean_response3=pd[3,"mean_response"], 
                   pd.mean_response4=pd[4,"mean_response"], 
                   pd.mean_response5=pd[5,"mean_response"], 
                   pd.mean_response6=pd[6,"mean_response"], 
                   pd.mean_response7=pd[7,"mean_response"], 
                   pd.mean_response8=pd[8,"mean_response"], 
                   pd.mean_response9=pd[9,"mean_response"], 
                   pd.mean_response10=pd[10,"mean_response"], 
                   pd.mean_response11=pd[11,"mean_response"], 
                   pd.mean_response12=pd[12,"mean_response"], 
                   pd.mean_response13=pd[13,"mean_response"], 
                   pd.mean_response14=pd[14,"mean_response"], 
                   pd.mean_response15=pd[15,"mean_response"], 
                   pd.mean_response16=pd[16,"mean_response"], 
                   pd.mean_response17=pd[17,"mean_response"], 
                   pd.mean_response18=pd[18,"mean_response"], 
                   pd.mean_response19=pd[19,"mean_response"], 
                   pd.mean_response20=pd[20,"mean_response"],
                   pd.stddev_response1 = pd[1, "stddev_response"], 
                   pd.stddev_response2 = pd[2, "stddev_response"], 
                   pd.stddev_response3 = pd[3, "stddev_response"], 
                   pd.stddev_response4 = pd[4, "stddev_response"], 
                   pd.stddev_response5 = pd[5, "stddev_response"], 
                   pd.stddev_response6 = pd[6, "stddev_response"], 
                   pd.stddev_response7 = pd[7, "stddev_response"], 
                   pd.stddev_response8 = pd[8, "stddev_response"], 
                   pd.stddev_response9 = pd[9, "stddev_response"], 
                   pd.stddev_response10 = pd[10, "stddev_response"], 
                   pd.stddev_response11 = pd[11, "stddev_response"], 
                   pd.stddev_response12 = pd[12, "stddev_response"], 
                   pd.stddev_response13 = pd[13, "stddev_response"], 
                   pd.stddev_response14 = pd[14, "stddev_response"], 
                   pd.stddev_response15 = pd[15, "stddev_response"], 
                   pd.stddev_response16 = pd[16, "stddev_response"], 
                   pd.stddev_response17 = pd[17, "stddev_response"], 
                   pd.stddev_response18 = pd[18, "stddev_response"], 
                   pd.stddev_response19 = pd[19, "stddev_response"], 
                   pd.stddev_response20 = pd[20, "stddev_response"],
                   pd.std_error_mean_response1 = pd[1, "std_error_mean_response"], 
                   pd.std_error_mean_response2 = pd[2, "std_error_mean_response"], 
                   pd.std_error_mean_response3 = pd[3, "std_error_mean_response"], 
                   pd.std_error_mean_response4 = pd[4, "std_error_mean_response"], 
                   pd.std_error_mean_response5 = pd[5, "std_error_mean_response"], 
                   pd.std_error_mean_response6 = pd[6, "std_error_mean_response"], 
                   pd.std_error_mean_response7 = pd[7, "std_error_mean_response"], 
                   pd.std_error_mean_response8 = pd[8, "std_error_mean_response"], 
                   pd.std_error_mean_response9 = pd[9, "std_error_mean_response"], 
                   pd.std_error_mean_response10 = pd[10, "std_error_mean_response"], 
                   pd.std_error_mean_response11 = pd[11, "std_error_mean_response"], 
                   pd.std_error_mean_response12 = pd[12, "std_error_mean_response"], 
                   pd.std_error_mean_response13 = pd[13, "std_error_mean_response"], 
                   pd.std_error_mean_response14 = pd[14, "std_error_mean_response"], 
                   pd.std_error_mean_response15 = pd[15, "std_error_mean_response"], 
                   pd.std_error_mean_response16 = pd[16, "std_error_mean_response"], 
                   pd.std_error_mean_response17 = pd[17, "std_error_mean_response"], 
                   pd.std_error_mean_response18 = pd[18, "std_error_mean_response"], 
                   pd.std_error_mean_response19 = pd[19, "std_error_mean_response"], 
                   pd.std_error_mean_response20 = pd[20, "std_error_mean_response"],
                   varimp.tmp.relative_importance = varimp[varimp[, "variable"] == "tmp", "relative_importance"],
                   varimp.tmp.scaled_importance = varimp[varimp[, "variable"] == "tmp", "scaled_importance"],
                   varimp.tmp.percentage = varimp[varimp[, "variable"] == "tmp", "percentage"],
                   varimp.sm.relative_importance = varimp[varimp[, "variable"] == "sm", "relative_importance"],
                   varimp.sm.scaled_importance = varimp[varimp[, "variable"] == "sm", "scaled_importance"],
                   varimp.sm.percentage = varimp[varimp[, "variable"] == "sm", "percentage"],
                   varimp.par.relative_importance = varimp[varimp[, "variable"] == "par", "relative_importance"],
                   varimp.par.scaled_importance = varimp[varimp[, "variable"] == "par", "scaled_importance"],
                   varimp.par.percentage = varimp[varimp[, "variable"] == "par", "percentage"],
                   varimp.pre.relative_importance = varimp[varimp[, "variable"] == "pre", "relative_importance"],
                   varimp.pre.scaled_importance = varimp[varimp[, "variable"] == "pre", "scaled_importance"],
                   varimp.pre.percentage = varimp[varimp[, "variable"] == "pre", "percentage"],
                   varimp.vpd.relative_importance = varimp[varimp[, "variable"] == "vpd", "relative_importance"],
                   varimp.vpd.scaled_importance = varimp[varimp[, "variable"] == "vpd", "scaled_importance"],
                   varimp.vpd.percentage = varimp[varimp[, "variable"] == "vpd", "percentage"],
                   threshold = threshold, modelid = modelid
                   
                   )
  return(df)
}
emptydf <- makeDF(makeDF(pd = pd, varimp = x1$varimp, threshold = x1$threshold, modelid = x1$modelid))[1,] <- NA
thres <- function(Loc){

  df <- data.frame(response=zcore(meco2[Loc,]),vpd=mvpd[Loc,],sm=msm[Loc,],pre=mpre[Loc,],tmp=mtmp[Loc,],par=mpar[Loc,])
  df <- na.omit(df)
#check data
  if(nrow(df) < 200 |any(apply(df,2,sd)==0))
    return(emptydf)
 # h2o.init()
  h2o.no_progress()
  df <- as.h2o(df)

  predictors <- c("vpd", "sm", "pre", "tmp", "par")
  response <- "response"
  
  # split into train and validation sets
  splits <- h2o.splitFrame(data =  df, ratios = 0.8, seed = 1234)
  train <- splits[[1]]
  test <- splits[[2]]
  # aml <- h2o.automl(x = predictors, y = response,
  #                   training_frame = train,
  #                   max_models = 10,
  #                   nfolds = 5,
  #                   exclude_algos = c("GLM", "DeepLearning"),
  #                   seed = 1234)
  # model <-  h2o.get_best_model(aml)
  
  model <- h2o.gbm(
    x = predictors,
    y = response,
    nfolds = 5,
    seed = 1234,
    training_frame = train
  )
  
  pd <- h2o.partialPlot(object = model,
                        data = test,
                        cols = "vpd",
                        plot=FALSE)
  varimp <- h2o.varimp(model)
  fitst <- chngptm(formula.1=mean_response~1, formula.2=~vpd,pd,type="stegmented", family="gaussian")
  th <- fitst$coefficients[[5]]
  modelid <-  model@model_id[1]
  
#  mylist <- list(pd = pd, varimp = varimp, threshold = th, modelid = modelid)
  
  mydata <- makeDF(pd = pd, varimp = varimp, threshold = th, modelid = modelid)
  h2o.removeAll()
  gc()
  return(mydata)
}
rr_gbm <- emptydf

for (Loc in line){

  rr <- thres(Loc)
  rownames(rr) <- Loc
  rr_gbm <- rbind(rr_gbm,rr)
  print(Loc)
 
}


######################Parallel with 1-core H2O cluster at a different port#############


library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(8)
registerDoSNOW(cl)

# Setup the progress bar
pb <- txtProgressBar(max = length(line), style = 3)
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress)

print(paste("Starting parellel processing at:", Sys.time()))
ptime  <- system.time(
  rr_gbm <- foreach(
    Loc = line,
    .options.snow = opts,
    .combine = rbind,
    .packages = c('h2o', 'chngpt')
  ) %dopar% {
    port <- 54321 + 3*length(Loc)
    print(paste0("http://localhost:", port))
    h2o.init(nthreads = 10, max_mem_size = "10G", port = port)
    thres(Loc)
  }
)

print(paste("\n Parellel processing complete at:", Sys.time()))

###stop cluster###
registerDoSEQ()
stopCluster(cl)


write.csv(rr_gbm,"/home/j_xiao/data/NCdemo/results/MODIS_NIRv/threshold2_h2o_gbm2.csv")
