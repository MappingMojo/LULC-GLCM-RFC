#This implementation of the Random Forest algorithm has been adapted from 
# Millard and Richardson, Remote Sens. 2015, 7(7), 8489-8515; doi:10.3390/rs70708489
##source: https://github.com/robingenuer/VSURF/commit/7c8d1243d39a83696690779668ea4cc217caca37

  #install packages
  install.packages('sp')
  install.packages('caret')
  install.packages('cbind')
  install.packages('magick')
  install.packages('snow')
  #install.packages("Rcpp")
  
  #load libraries
  library(raster)
  library(randomForest)
  library(sp)
  library(rgdal)
  library(ggplot2)
  library (VSURF)
  library(lattice)
  library(dplyr)
  library(parallel)
  library(caret)
  library(magick)
  library(e1071)
  library(MASS)
  #library(varSelRF) 

#===========================================================================================
#set working directory to data folder where your data is stored
#===========================================================================================
  
  setwd("D:/R_Data Analysis/LUC")
  file.edit(file.path("~", ".Rprofile"))  #invisible(utils::memory.limit(size = 60000))
  
  # Start the clock!
  ptm <- proc.time()

#===========================================================================================
#Read input data (image and label data) source: #https://stackoverflow.com/questions/34513734/how-do-i-stack-raster-files-in-r
#===========================================================================================
  files <- list.files(path="D:/R_Data Analysis/LUC/Dataset1", pattern="*.tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
  allInfo = image_info(image_read(files))   
  allInfo$fileName = files         #Attach the file names
  files
  infile <-  stack(files) 

#==========================================================================================
#normalize bands of  dataset
#==========================================================================================  
  
  norm <- function(x){(x-min)/(max-min)}
  
  for(j in 1:nlayers(infile)){
    
    cat(paste("Currently processing layer:", j,"/",nlayers(infile), "\n"))
    
    min <- cellStats(infile[[j]],'min')
    max <- cellStats(infile[[j]],'max')
    
    #initialize cluster
    #number of cores to use for clusterR function (max recommended: ncores - 1)
    beginCluster(9)
    
    #normalize
    infile[[j]] <- clusterR(infile[[j]], calc, args=list(fun=norm), export=c('min',"max"))
    
    #end cluster
    endCluster()
  }
  
  infileVS.1 <-stack(infile)
  infileVS.1
  
  ##reset normalized raster layer names to original names
  infile <-  stack(files) 
  infile
  names(infileVS.1) = names(infile) 
  names(infileVS.1)
  
  
  ##Coerce raster stack to spatial grid dataframe (first convert it to a spatial pixels dataframe)
  r <- as(as(infileVS.1, 'SpatialPixelsDataFrame'), 'SpatialGridDataFrame')
  infileVS <- readGDAL(r)
  infileVS
  
  
  ##Write/save the file 
  # outpath <- ("D:/R_Data Analysis/LUC/")
  # outraster = paste(outpath, r, "1.tif")
  # writeRaster(r, filename = outraster, format='GTiff',
  #             options='COMPRESS=LZW', datatype='FLT4S', overwrite=TRUE)
  
  #==========================================================================================
#Input the training and validation shapefiles 
#==========================================================================================
  
  Path = ("D:/R_Data Analysis/LUC/TrainTest")
  TR.shape.file="Train"
  VS.shape.file="Validation"
  
  #===========================================================================================
  #Feature selection using VSURF
  #===========================================================================================
  
  #selection of bands in image you want to use 
  selection<-c(1:nlayers(infileVS)) 
  
  #===========================================================================================
  #Read polygons from shape file (training set)
  #===========================================================================================
  
  VTS <- readOGR(dsn=Path, layer=TR.shape.file,  stringsAsFactors = FALSE)
  np <- length(VTS@polygons)
  
  names(VTS@data)[11] <- "lulc"  #column with class values
  VTS@data$lulc <- as.numeric(VTS@data$lulc)
  
  Ncl <- length(table(VTS@data$lulc))
  
  #==================================================================================
  # Subset image by polygons
  #==================================================================================
  
  xy <- coordinates(infileVS)[,]
  
  lulc <- array(0,0)
  pix_id <- array(0,0)
  
  
  for(id in 1:np)
  {
    xy_pol <- VTS@polygons[[id]]@Polygons[[1]]@coords
    tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_pol[,1],xy_pol[,2])
    ind <- which(tmp.in>0)
    lulc <- c(lulc,rep(VTS$lulc[id],length(ind)))
    pix_id <- c(pix_id,ind)
  }
  
  Labels <- data.frame(lulc=lulc)
  

  TR <- cbind(infileVS@data[pix_id,],lulc=lulc)
  TR
  
  training_response <- as.factor(TR[,"lulc"])
  training_response

  training_predictors = as.factor(TR[,selection])
  training_predictors

#==========================================================================================
# Create the VSURF model
#==========================================================================================
  # memory.limit(size=56000) # increase RAM memory capacity
  
  x=training_predictors
  y=training_response
  thesis.vsurf<- VSURF(x, y, ntree = 2000, verbose=TRUE, RFimplem = "randomForest", parallel = FALSE)
                       #,nfor.thres=100,nfor.interp=40,nfor.pred=18) #,ncores=20)  
  ##OR
  # thesis.vsurf<- VSURF(x=TR[,selection], y=TR$lulc, ntree = 2000, verbose=TRUE, RFimplem = "randomForest", parallel = FALSE,
  #                 nfor.thres=100,nfor.interp=40,nfor.pred=18) #,ncores=20)  
  
  thesis.vsurf
  summary(thesis.vsurf)
  plot(thesis.vsurf, step = "thres", imp.sd = FALSE, var.names = TRUE, ylim = c(0, 2e-4))
  thesis.vsurf$varselect.pred
  names(infile[, c(thesis.vsurf$varselect.pred)])
  
  #==========================================================================================
  #### Plot VSURF results
  plot(thesis.vsurf, step = "all", var.names = TRUE,
       imp.mean = TRUE, imp.sd = TRUE,
       nvar.imp.mean = length(thesis.vsurf$imp.mean.dec),
       nvar.imp.sd = length(thesis.vsurf$imp.sd.dec),
       nvar.interp = length(thesis.vsurf$varselect.thres),
       nvar.pred = length(thesis.vsurf$varselect.pred))
  
  ## S3 method for class 'VSURF_thres'
  plot(thesis.vsurf, var.names = TRUE, imp.mean = TRUE,
       imp.sd = TRUE, nvar.imp.mean = length(thesis.vsurf$imp.mean.dec),
       nvar.imp.sd = length(thesis.vsurf$imp.sd.dec))
  
  ## S3 method for class 'VSURF_interp'
  plot(thesis.vsurf, var.names = TRUE,
       nvar.interp = length(thesis.vsurf$varselect.thres))
  
  ## S3 method for class 'VSURF_pred'
  plot(thesis.vsurf, var.names = TRUE,
       nvar.pred = length(thesis.vsurf$varselect.pred))
  
  #==========================================================================================    
  #==========================================================================================    
  #variable importance plot   
  
  fdir <- ("D:/R_Data Analysis/LCC/VSURF_Plots/")
  scaleFUN = function(x) sprintf("%.3f", x)
  
  #plot variables selected for prediction (stage III)
  fileout  = "vsurf_pred_lcc"
  indent3 = data.frame(cbind(colnames(training_predictors)[thesis.vsurf$varselect.pred]),thesis.vsurf$imp.mean.dec[1:length(thesis.vsurf$varselect.pred)])
  indent3
  
  colnames(indent3)= c("Var","Imp")
  indent3$Imp = sapply(indent3$Imp,as.numeric)
  barP = ggplot(data=indent3,aes(reorder(Var,Imp),Imp)) + geom_bar(fill="gold",stat="identity",position="dodge") + coord_flip() #horizontal bar graph
  barP = barP + scale_y_continuous(expand=c(0,0),limits=c(0,0.07),breaks=seq(0,0.07,by=0.01),labels=scaleFUN(seq(0,0.07,by=0.01)),name=expression(paste("VI Mean Decrease")))
  barP = barP + scale_x_discrete(name="Predictor Variables",expand=c(0,0))
  barP = barP + theme(axis.title.x=element_text(family="mono",color="black",size=11,vjust=0.1),axis.title.y=element_text(family="mono",color="black",size=14,vjust=1)) #modified axies lables
  barP = barP + theme(axis.ticks=element_blank(),axis.text.x=element_text(family="mono",color="black",size=11),axis.text.y=element_text(family="mono",color="black",size=14)) #modifies tick labels
  barP = barP + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,size=0.5,colour="gray65")) #removes standard ggplot background
  barP = barP + theme(plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm")) #formats area
  
  ggsave(filename=paste(fdir,fileout,".importance.RF",".tiff",sep=""),plot=barP,dpi=600,limitsize=T)
  plot(barP)
  
  #plot variables selected for interpretation (stage II)
  fileout  = "vsurf_interp_lcc"
  ident = data.frame(cbind(colnames(training_predictors)[thesis.vsurf$varselect.interp]),thesis.vsurf$imp.mean.dec[1:length(thesis.vsurf$varselect.interp)])
  ident
  
  colnames(ident)= c("Var","Imp")
  ident$Imp = sapply(ident$Imp,as.numeric)
  barP = ggplot(data=ident,aes(reorder(Var,Imp),Imp)) + geom_bar(fill="orange",stat="identity",position="dodge") + coord_flip() #horizontal bar graph
  barP = barP + scale_y_continuous(expand=c(0,0),limits=c(0,0.07),breaks=seq(0,0.07,by=0.01),labels=scaleFUN(seq(0,0.07,by=0.01)),name=expression(paste("VI Mean Decrease")))
  barP = barP + scale_x_discrete(name="Interpret Variables",expand=c(0,0))
  barP = barP + theme(axis.title.x=element_text(family="mono",color="black",size=11,vjust=0.1),axis.title.y=element_text(family="mono",color="black",size=14,vjust=1)) #modified axies lables
  barP = barP + theme(axis.ticks=element_blank(),axis.text.x=element_text(family="mono",color="black",size=11),axis.text.y=element_text(family="mono",color="black",size=14)) #modifies tick labels
  barP = barP + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,size=0.5,colour="gray65")) #removes standard ggplot background
  barP = barP + theme(plot.margin=unit(c(0.8,0.8,0.8,0.8),"cm")) #formats area
  
  
  ggsave(filename=paste(fdir,fileout,".importance.RF",".tiff",sep=""),plot=barP,dpi=600,limitsize=T)
  plot(barP)
  
  #===============================================================================================================================================================
  #plot all variables showing the number selected at each stage (threshold, interpretation and prediction)
  
  fileout  = "vsurf_err_lcc"
  ident.MSE                             = data.frame(seq(1:length(thesis.vsurf$err.interp)),thesis.vsurf$err.interp)
  #ident.MSE                             = data.frame(seq(1:length(humid.RFmodel150new[["err.interp"]])),humid.RFmodel150new[["err.interp"]])
  colnames(ident.MSE)                   = c("VarNum","OOB")
  lineP                                 = ggplot(ident.MSE,aes(x=VarNum,y=OOB)) + geom_line(color="black",size=0.7)
  lineP                                 = lineP + geom_vline(xintercept=length(thesis.vsurf$varselect.pred),color="gold",size=0.5)
  lineP                                 = lineP + geom_vline(xintercept=length(thesis.vsurf$varselect.interp),linetype="dotted", color="orange",size=0.5)
  lineP                                 = lineP + geom_vline(xintercept=length(thesis.vsurf$varselect.thres),linetype="dotted", color="red",size=0.5)
  lineP                                 = lineP + annotate("text",x=(xintercept=length(thesis.vsurf$varselect.pred)+1),y=0.50,label="Prediction",vjust=0.7,family="sans",fontface=1,size=5,color="gold",angle=90)
  lineP                                 = lineP + annotate("text",x=(xintercept=length(thesis.vsurf$varselect.interp)+1),y=0.50,label="Interpretation",vjust=0.7,family="sans",fontface=1,size=5,color="orange",angle=90)
  lineP                                 = lineP + annotate("text",x=(xintercept=length(thesis.vsurf$varselect.thres)+1),y=0.70,label="Threshold",vjust=0.7,family="sans",fontface=1,size=5,color="red",angle=90)
  lineP                                 = lineP + scale_x_continuous(expand=c(0,0),limits=c(0,nrow(ident.MSE)),breaks=seq(0,nrow(ident.MSE),by=3),labels=seq(0,nrow(ident.MSE),by=3),name=("No. of Variables"))
  lineP                                 = lineP + scale_y_continuous(expand=c(0,0),limits=c(0,0.99),breaks=seq(0,0.99,by=0.10),labels=scaleFUN(seq(0,0.99,by=0.10)), name=expression("Mean OOB Error Rate")) ### CHANGE THE 0.10 by 0.005
  lineP                                 = lineP + theme(axis.ticks=element_line(size=0.5,color="black"),axis.text.x=element_text(margin=ggplot2::margin(0.5,0.5,1,1,"cm")),axis.text.y=element_text(margin=ggplot2::margin(0.5,0.5,1,1,"cm")),axis.ticks.length=unit(-0.25,"cm"),legend.position="none")
  lineP                                 = lineP + theme(axis.title.x=element_text(family="sans",color="black",size=9,vjust=0.05),axis.title.y=element_text(family="sans",color="black",size=12,vjust=0.2)) #modified axies lables
  lineP                                 = lineP + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),panel.background=element_blank(),panel.border=element_rect(fill=NA,size=0.5,colour="black")) #removes standard ggplot background
  lineP                                 = lineP + theme(plot.margin=unit(c(0.6,0.6,0.6,0.6),"cm")) #formats area
  ggsave(filename=paste(fdir,fileout,".error.RF",".tiff",sep=""),plot=lineP,dpi=600,limitsize=T)
  plot(lineP)
  
  
#===========================================================================================
#Join other data (Land-Cover band)
#===========================================================================================
  
  #subset selected VSURF variables
  infileR <- infileVS[[thesis.vsurf$varselect.pred]]
  infileR
  
 
  #Read land cover data and join with normalized selected bands
  files <- list.files(path="D:/R_Data Analysis/Output/Multi-class1", pattern="rf_iteration_7_NumberForTraining_850_.pix", all.files=FALSE, full.names=TRUE,recursive=TRUE)
  allInfo = image_info(image_read(files))
  #Attach the file names
  allInfo$fileName = files
  dim(files)
  LCC <- stack(files)

  
  #join the two datasets
  infileRF<- stack(other.data, infileR)
  dim(infileRF)
  names(infileRF)
  infileRF
  
#===========================================================================================
#Set up RFC model including splitting training and validation data
#===========================================================================================
  
    #selection of bands in image you want to use 
    selection<-c(1:nlayers(infileRF)) 
    
    NumIterations<-1
    #number of times to run the classification
  
    cumulative_errorTable = NULL    
    #this variable will hold the error table populated in subsequent iterations of the loop

  for (i in 1:NumIterations){
    cat("Running classification iteration ", i)
  
#===========================================================================================
#Read polygons from shape file (training set)
#===========================================================================================
    
    VTS <- readOGR(dsn=Path, layer=TR.shape.file,  stringsAsFactors = FALSE)
    np <- length(VTS@polygons)
    
    names(VTS@data)[11] <- "lulc"
    VTS@data$lulc <- as.numeric(VTS@data$lulc)
    
    Ncl <- length(table(VTS@data$lulc))
    
 
    # Subset image by polygons
    xy <- coordinates(infileRF)[,]
    
    lulc <- array(0,0)
    pix_id <- array(0,0)
    
    
    for(id in 1:np)
    {
      xy_pol <- VTS@polygons[[id]]@Polygons[[1]]@coords
      tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_pol[,1],xy_pol[,2])
      ind <- which(tmp.in>0)
      lulc <- c(lulc,rep(VTS$lulc[id],length(ind)))
      pix_id <- c(pix_id,ind)
    }
    
    Labels <- data.frame(lulc=lulc)
    
    # Extract the training data 
    TR <- cbind(infileRF@data[pix_id,],Labels)
    TR
    
    #training_response = as.factor(TR.shape.file$lulc) 
    training_responseRF = as.factor(TR$lulc)
    training_responseRF
    #here, our csv file has a field called "lulc" for names/code of each of the classes 
    
    training_predictorsRF = TR[,selection] 
    training_predictorsRF
    
    # Create and save the forest
    ntree <- 100    #number of trees to produce per iteration
    #mtry = 5
    r_forest = randomForest(x=training_predictorsRF, y=training_responseRF, ntree = ntree, keep.forest=TRUE, 
                            importance = TRUE, proximity=TRUE) 
    
    #===========================================================================================
    #Investigate the OOB (Out-Of-the bag) error
    #===========================================================================================
    r_forest
    
    #===========================================================================================
    # Assessment of variable importance
    #===========================================================================================
    imp<-importance(r_forest)  #for ALL classes individually
    imp   #display importance output in console
    varImpPlot(r_forest)
    
    #=========================================================================================
    #Evaluate the impact of the ntree on the accuracy
    #========================================================================================
    tree_nr_evaluation <-data.frame(
      Trees=rep(1:nrow(r_forest$err.rate), times=7),
      Type=rep(c("OOB", "1", "2", "3", "7", "8", "9"),
               each=nrow(r_forest$err.rate)),
      Error = c(r_forest$err.rate[, "OOB"],
                r_forest$err.rate[, "1"],
                r_forest$err.rate[, "2"],
                r_forest$err.rate[, "3"],
                r_forest$err.rate[, "7"],
                r_forest$err.rate[, "8"],
                r_forest$err.rate[, "9"]))
    
    ggplot(data=tree_nr_evaluation, aes(x=Trees, y=Error)) + geom_line(aes(color=Type))
    
    #=======================================================================================
    #Evaluate the impact of the mtry on the accuracy
    #========================================================================================
    mtry_evaluation <-vector (length=nlayers(infileRF))
    for(i in 1:nlayers(infileRF)){
      test_model <- randomForest(training_predictorsRF, y=training_responseRF, mtry=i, ntree = ntree, keep.forest=TRUE, importance = TRUE, proximity=TRUE) 
      mtry_evaluation[i] <-test_model$err.rate[nrow(test_model$err.rate),1]
    }
    mtry_evaluation
    #==========================================================================================
    
    #set up the file name for the raster to print
    outraster = paste("rf_iterationLUC", i, "NumberForTraining", NumberForTraining, ".pix", sep = "_")
    
    #Classify the whole image
    #gets the raster data to use for classification
    inraster_rf = stack(infileRF)
    predictor_data = subset(inraster_rf, selection)	
    
    #change this to your output directory if different 
    setwd("D:/R_Data Analysis/Output") 
    
    #write out the predicted classification to a raster file
    predictions = predict(predictor_data, r_forest, filename=outraster, format="PCIDSK", overwrite=TRUE, progress="text", type="response")
    remove(predictor_data) 
    
    #do independent validation
    proj4string(outraster)<-CRS("+proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs") 
    
    #reads polygons from the classification you just produced 
    VVS <- readOGR(dsn=Path, layer=VS.shape.file,  stringsAsFactors = FALSE)
    VVS
    np_valid <- length(VVS@polygons)
    
    names(VVS@data)[11] <- "lulc"
    VVS@data$lulc <- as.numeric(VVS@data$lulc)
    
    Ncl.1 <- length(table(VVS@data$lulc))
    
    #==================================================================================
    # Subset image by polygons
    #==================================================================================
    
    xy <- coordinates(outraster)[,]
    
    lulc <- array(0,0)
    pix_id <- array(0,0)
    
    
    for(id in 1:np_valid)
    {
      xy_pol <- VVS@polygons[[id]]@Polygons[[1]]@coords
      tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_pol[,1],xy_pol[,2])
      ind <- which(tmp.in>0)
      lulc <- c(lulc,rep(VVS$lulc[id],length(ind)))
      pix_id <- c(pix_id,ind)
    }
    
    Labels.1 <- data.frame(lulc=lulc)
    
    # Extract the validation data 
    VS <- cbind(outraster@data[pix_id,],Labels.1)
    VS
    
    #extracts the value of the classified raster at the validation point locations  
    
    ValidationJoin<-cbind(VVS$lulc,VS) 
    dim(ValidationJoin)
    ValidationJoin
    #joins the predicted with observed
    
    totalcorrect<-(subset(ValidationJoin,ValidationJoin[,1]==ValidationJoin[,2])) 
    #determines which samples were the same for observed and predicted
    
    indep_errorrate<-100 - (length(totalcorrect[,1])/length(validation)*100)
    #determines the independent error rate (%)
    
    rfOOB<-r_forest$err.rate[ntree,1]*100 #obtains the rfOOB error rate 
    differror<-rfOOB - indep_errorrate
    current_error<-cbind(names(outraster), rfOOB, indep_errorrate, differror)
    #binds all rfOOB, indep_errorrate, and differror together
    
    cumulative_errorTable<-rbind(current_error, cumulative_errorTable)	
    # populates a variable of cumulatively collected validation data within the loop
    
  } #end of loop
  
  current_error

#==========================================================================================
#confusion matrix & Overall Accuracy
#==========================================================================================
  cm = table(prediction=VS,truth=VVS$lulc)
  oa =sum(diag(cm))/(length(VVS$lulc))
  cm
  oa
  tableName = paste("rf_iteration", i,"NumberForTraining", NumberForTraining, ".csv", sep = "_")
  write.csv(cumulative_errorTable, file = tableName)	
  #prints out a table with error from all classifications
  warning()

#==========================================================================================
####F1 Score using CARET
#==========================================================================================
  #F1 <- (2 * precision * recall) / (precision + recall)
  
  predicted <- as.factor(c(VS))
  realized  <- as.factor(c(VVS$lulc))
  
  # Compute the confusion matrix and all the statistics
  result <- confusionMatrix(predicted, realized, mode="prec_recall")
  
  result
  result$byClass["Precision"]
  result$byClass["Recall"]
  result$byClass["F1"]
  
  # Stop the clock
  proc.time() - ptm
#==========================================================================================
##################################### THE END ############################################
#==========================================================================================

