##source: https://riptutorial.com/r/example/12880/calculating-glcm-texture

install.packages('glcm')
install.packages('sp')

library(glcm)
library(sp)
library(raster)

#inpath <- ("C:/Users/pnkay/Downloads/Thesis/Data/Final4/B8_Nairobi.tif")
#outpath <- ("C:/Users/pnkay/Downloads/Thesis/Data/Final4/Text21_nir")

r <- raster("C:/Users/pnkay/Downloads/Thesis/Data/Final4/B8_Nairobi.tif")

plot(r)

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

names(rglcm.1) <- c("mean", "variance", "homogeneity", "contrast", 
                    "dissimilarity", "entropy", "second_moment", "correlation")


setwd("C:/Users/pnkay/Downloads/Thesis/Data/Final4/B8/")

for(i in 1:nlayers(rglcm.1)){
  band<-rglcm.1[[i]]
  #save raster in a separate file
  writeRaster(band,paste(names(rglcm.1[[i]]), "B8_27",  '.tif', sep=""))
}


remove(rglcm) 

################################################################################