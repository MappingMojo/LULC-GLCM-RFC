##source: https://riptutorial.com/r/example/12880/calculating-glcm-texture
#install required packages
install.packages('glcm')
install.packages('sp')

#read libraries
library(glcm)
library(sp)
library(raster)

#inpath <- ("path")
#outpath <- ("path")

#read raster from path (only if a single raster file is in the folder)
r <- raster("path")

plot(r)

#Generate GLCM textures
rglcm <- glcm(r, 
              window = c(27,27), 
              shift = c(1,1), 
              statistics = c("mean", "variance", "homogeneity", "contrast", 
                             "dissimilarity", "entropy", "second_moment", "correlation"),
              na_val = 0
)

rglcm.1 <- stack(rglcm)
plot(rglcm.1)
nlayers(rglcm.1)

#name variables
names(rglcm.1) <- c("mean", "variance", "homogeneity", "contrast", 
                    "dissimilarity", "entropy", "second_moment", "correlation")

#set output directory
setwd("path")

#save each variable/band as a seperate file
for(i in 1:nlayers(rglcm.1)){
  band<-rglcm.1[[i]]
  #save raster in a separate file
  writeRaster(band,paste(names(rglcm.1[[i]]), "B8_27",  '.tif', sep=""))
}

#offload data
remove(rglcm) 

################################################################################
