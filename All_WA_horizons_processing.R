library(terra)
library(sf)
library(dplyr)
library(tidyterra)

setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")

##### basic dataframe####
wa_hor <- read.csv("SOIL CARBON/ANALYSIS/All_WA_horizons.csv")
##### #convert to spatial vector points##### 
wa_hor_pts <- terra::vect(wa_hor, geom = c("lon", "lat"), crs = "EPSG:4326") 

##### #The larger rasters need to stay in epsg:4326 because projecting takes too long##### 
All_NDYI <- merge(rast("SOIL CARBON/All_WA/WA_Spectral_NDYI-0000000000-0000000000.tif"), rast("SOIL CARBON/All_WA/WA_Spectral_NDYI-0000000000-0000023296.tif")) 
All_geomorph <- rast("WA_DNR/Project/hydrology/wa_geomorphons.tif")
large_rasts <- c(All_NDYI, All_geomorph)

#####  #extract large raster metrics ##### 
#simple extract 
wa_hor_pts_ext <- terra::extract(large_rasts, wa_hor_pts, ID = F, bind = T)

#### Climate
wa_clim <- c(rast("PRISM_Climate/PRISM_tmean_30yr_normal_800mM4_annual_asc/PRISM_tmean_30yr_normal_800mM4_annual_asc.asc"),rast("PRISM_Climate/PRISM_ppt_30yr_normal_800mM4_annual_asc/PRISM_ppt_30yr_normal_800mM4_annual_asc.asc") )
wa_hor_pts_ext <- terra::extract(wa_clim, wa_hor_pts_ext, ID = F, bind = T, method = "bilinear")



##### #bring in the study area raster data for covariates##### 
    #need to merge the GEE layers
hoh_spec <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/GEE20yr/hoh_spec_20yr_sea2.tif") |> select(NDVI, MNDWI, EVI, SCI) 
mas_spec <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/GEE20yr/mas_spec_20yr_sea2.tif")|> select(NDVI, MNDWI, EVI, SCI) 
col_spec <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/GEE20yr/col_spec_20yr_sea2.tif")|> select(NDVI, MNDWI, EVI, SCI)
gee_spec <- merge(hoh_spec, mas_spec, col_spec)

#####  #extract gee metrics ##### 
#simple extract 
wa_hor_pts_ext <- terra::extract(gee_spec, wa_hor_pts_ext, ID = F, bind = T, method = "bilinear")

#geology 
wa_geo <- vect("WA_Geo/WA_Geology_100K.gpkg") |> select(GEOLOGIC_AGE) |> terra::project("EPSG:4326")
#wa_hor_pts_ext <- terra::extract(wa_geo, wa_hor_pts_ext)
#wa_hor_pts_geo
wa_hor_pts_ext$GEOLOGIC_AGE <- wa_hor_pts_geo$GEOLOGIC_AGE

#WIP has been merged & in epsg:26910
All_WIP <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/ALL_WA_WIP_MERGE.tif")

wa_hor_pts_extwip <- terra::project(wa_hor_pts_ext, "EPSG:26910")
wa_hor_pts_extwip <- terra::extract(All_WIP, wa_hor_pts_extwip, ID = F, bind = T, method = "bilinear")




wa_hor_pts_extwip_df <- terra::as.data.frame(wa_hor_pts_extwip, geom = "XY") |> 
    dplyr::rename(NDYI = NDYI_median,
                  Geomorphon_Class = wa_geomorphons,
                  WIP = WET,
                  Temp = PRISM_tmean_30yr_normal_800mM4_annual_asc,
                  Precip = PRISM_ppt_30yr_normal_800mM4_annual_asc)


write.csv(wa_hor_pts_extwip_df, "SOIL CARBON/ANALYSIS/All_WA_horizons_spec_geoage.csv")
