library(blockCV)
library(raster)
library(sf)
library(ggplot2)


# Load package data -------------------------------------------------------

# package data (rsaster covariates of Austrailian Wet Tropic region)
awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))

#presence-absence data includes 116 presence points and 138 absence points
#appropriate format of species data for blockCV packages is sf or SpatialPointsDataFrame
# import presence-absence data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))

# make a SpatialPointsDataFrame object from data.frame
pa_data <- st_as_sf(PA, coords = c("x", "y"), crs = crs(awt))

# plot points on map
plot(awt[[1]]) # plot raster data
plot(pa_data[which(pa_data$Species==1), ], pch = 16, col="red", add=TRUE) # add presence points
plot(pa_data[which(pa_data$Species==0), ], pch = 16, col="blue", add=TRUE) # add absence points
legend(x=500000, y=8250000, legend=c("Presence","Absence"), col=c(2, 4), pch=c(16,16), bty="n")

# import presence-background species data - 116 presence points, 10,000 random backround points as 0s

PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pb_data <- st_as_sf(PB, coords = c("x", "y"), crs = crs(awt))
# number of presence and background records
table(pb_data$Species)

# Spatial Blocking -----------------------------------------------------------
# by specified range with random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species",
                   rasterLayer = awt,
                   theRange = 70000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)

# spatial blocking by rows and columns with checkerboard assignment
sb2 <- spatialBlock(speciesData = pb_data, # presence-background data
                    species = "Species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 6,
                    k = 5,
                    selection = "systematic",
                    biomod2Format = TRUE)

# spatial blocking by rows with systematic assignment
sb3 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    rasterLayer = awt,
                    rows = 6,
                    selection = "checkerboard",
                    biomod2Format = TRUE)

# adding points on saptialBlock plot
sb$plots + geom_sf(data = pa_data, alpha = 0.5)

# Buffers -----------------------------------------------------------
# buffering with presence-absence data
bf1 <- buffering(speciesData = pa_data,
                  theRange = 70000,
                  species = "Species", # to count the number of presences and absences/backgrounds
                  spDataType = "PA", # presence-absence  data type
                  progress = TRUE)

 # buffering with presence-background data (takes a while)
bf2 <- buffering(speciesData = pb_data, # presence-background data
                  theRange = 70000,
                  species = "Species",
                  spDataType = "PB", # presence-background data type
                  addBG = TRUE, # add background data to testing folds
                  progress = TRUE)

 # environmental clustering
 eb <- envBlock(rasterLayer = awt,
                speciesData = pa_data,
                species = "Species",
                k = 5,
                standardization = "standard", # rescale variables between 0 and 1
                rasterBlock = FALSE,
                numLimit = 50)
 
 # Spatial Autocorr -----------------------------------------------------------
 # measure spatial autocorrelation in the predictor raster files
 sac <- spatialAutoRange(rasterLayer = awt,
                         sampleNumber = 5000,
                         doParallel = TRUE,
                         showPlots = TRUE)
 # class of the output result
 class(sac)
 
 # summary statistics of the output
 summary(sac)
 
 plot(sac$variograms[[1]]) # doesn't work
 
# Fold explorer Shiny App -----------------------------------------------------------
 
# explore generated folds
foldExplorer(blocks = sb,
               rasterLayer = awt,
               speciesData = pa_data)
 
# Range explorer Shiny App -----------------------------------------------------------
# explore the block size
rangeExplorer(rasterLayer = awt) # the only mandatory input

# add species data to add them on the map
rangeExplorer(rasterLayer = awt,
              speciesData = pa_data,
              species = "Species",
              rangeTable = NULL,
              minRange = 30000, # limit the search domain
              maxRange = 100000)

 
# Maxent modeling -----------------------------------------------------------
 # loading the libraries
 library(maxnet)
# install.packages("precrec")
 library(precrec)
 # library(ggplot2)

 # extract the raster values for the species points as a dataframe
 mydata <- raster::extract(awt, pb_data)
 mydata <- as.data.frame(mydata)
 # create a vector of 1 (for presence) and 0 (for background samples)
 pb <- pb_data$Species

 # extract the folds in spatialBlock object created
 # in the previous section (with presence-background data)
 # the foldID only works for spatialBlock and envBlock folds
 folds <- sb2$foldID

 # create an empty vector to store the AUC of each fold
 AUCs <- vector(mode = "numeric")
 for(k in seq_len(5)){
   # extracting the training and testing indices
   # this way only works with foldID
   trainSet <- which(folds != k) # training set indices
   testSet <- which(folds == k) # testing set indices
   # fitting a maxent model using linear, quadratic and hinge features
   mx <- maxnet(p = pb[trainSet],
                data = mydata[trainSet, ],
                maxnet.formula(p = pb[trainSet],
                               data = mydata[trainSet, ],
                               classes = "default"))
   testTable <- pb_data[testSet, ] # a table for testing predictions and reference data
   testTable$pred <- predict(mx, mydata[testSet, ], type = "cloglog") # predict the test set
   # calculate area under the ROC curve
   precrec_obj <- evalmod(scores = testTable$pred, labels = testTable$Species)
   AUCs[k] <- auc(precrec_obj)[1,4] # extract AUC-ROC
 }

 # print the mean of AUCs
 print(mean(AUCs))
 
 # The model fitting is not run to save the vignette generation time
 # this AUC is based on the actual run
 print(0.8664762)
 
