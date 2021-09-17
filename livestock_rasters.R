
# =============================================================================
# Name:           livestock_rasters.R
# Description:    takes an SDM prediction raster and an FAO livestock abundance raster and transforms them to the same resolution.  Also outlines how to filter livestock abundance by a threshold.
# 
# Inputs:         monitoring.csv
#                 point_check.csv

# Outputs:        labels_mk.csv
# 
# Notes:          - 
#                 - 
# =============================================================================


# load packages -----------------------------------------------------------
library(raster)
library(rgdal)

# load data ---------------------------------------------------------------
# loading an example prediction raster we got for elephants from Wallace (example data not mean to be used for analysis)
eleph_raster <- raster("livestock/example_wallace_raster/layer.grd")
plot(eleph_raster)
# checking the crs
crs(eleph_raster)  # +proj=longlat +datum=WGS84 +no_defs

# loading global 2010 cattle data
# data originally from: http://www.fao.org/livestock-systems/global-distributions/cattle/en/
cattle <- raster("livestock/fao_datasets/global_cattle/6_Ct_2010_Aw.tif")
plot(cattle)
crs(cattle) # +proj=longlat +datum=WGS84 +no_defs

# checking resolution of livestock data and prediction raster
res(cattle) # 5 arc minutes
res(eleph_raster) # 2.5 arc minutes

# projectRaster() to crop and resample cattle raster at the prediction raster resolution
cattle_crop <- projectRaster(cattle, eleph_raster, method = "bilinear")
res(cattle_crop) # 2.5 arc minutes

# they can now be plotted together
plot(cattle_crop, legend = FALSE)
plot(eleph_raster, add = TRUE, legend = FALSE)

# Might want to focus on areas of high density, so extracting pixels where suitability and cattle density are over a cretain threshold
cattle_filter <-  cattle_crop > 5000 # 5000 chosen at random, if going this route, you would want to pick an more informative value
cattle_mask <- mask(cattle_crop, cattle_filter, maskvalue = 0) # crop the raster to high density areas
plot(cattle_mask)

# from here, you could identify pixels where there are overlap, maybe combine with land use data information, look at changes between current vs. predicted ranges, etc.