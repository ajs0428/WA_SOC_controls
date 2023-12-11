library(terra)
library(tidyterra)

setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")


#### Points ####
wa_dat <- read.csv("SOIL CARBON/ANALYSIS/All_WA_horizons_spec_geoage.csv") |>
    dplyr::select(sample_ID, site, depth_cm, BD_g_cm3, carbon_perc, carbon_stock_g_cm2, x, y)
str(wa_dat)

wa_dat_pts <- vect(wa_dat, geom = c("x", "y"), crs = "EPSG:26910")
plot(wa_dat_pts)

#### Raster and Vector data ####
# 
# catconv <- data.frame(id = 1:10, geomorphons = letters[1:10])
# 
# #Colville
# col_hli <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/col_HLI.tif") |> project("EPSG:26910")
# col_geomorph <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/col_geomorphons_11_23.tif")
# col_geomorph_cat <- categories(col_geomorph, layer = 1, value = catconv) |> project("EPSG:26910")
# 
# 
# #Mashel
# mas_hli <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mas_HLI.tif") |> project("EPSG:26910")
# mas_geomorph <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mas_geomorphons_11_23.tif")|> 
#     project("EPSG:26910") 
# mas_geomorph_cat <- categories(mas_geomorph, layer = 1, value = catconv) |> project("EPSG:26910")
# 
# 
# #Hoh
# hoh_hli <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_HLI.tif") |> project("EPSG:26910")
# hoh_geomorph <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_geomorphons_11_23.tif")
# hoh_geomorph_cat <- categories(hoh_geomorph, layer = 1, value = catconv) |> project("EPSG:26910")

#Merge Layers
geomorph_stack <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_geomorph_merge.tif")#terra::merge(col_geomorph_cat, mas_geomorph_cat, hoh_geomorph_cat, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_geomorph_merge.tif")
hli_stack <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_HLI_merge.tif")#terra::merge(col_hli, mas_hli, hoh_hli, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_HLI_merge.tif")

merge_wip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_WIP_MERGE.tif")

#All WA layers
wa_state <- vect("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_State_Boundary/WA_State_Boundary.shp") |> project("EPSG:4326")

wa_precip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_Precip.tif") |> terra::project("EPSG:26910")
wa_temp <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/ALL_WA_Temp.tif") |> terra::project("EPSG:26910")

wa_geology <- vect("WA_Geo/ger_portal_surface_geology_100k/surface_geology_100k.gdb", layer = "geologic_unit_poly_100k") |> 
    terra::project("EPSG:26910") |> tidyterra::select(LITHOLOGY, GEOLOGIC_AGE)

wa_spec1 <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_median_region/WA_median_regions-0000000000-0000000000.tif") 
wa_spec2 <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_median_region/WA_median_regions-0000000000-0000019968.tif")#
wa_spec3 <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_median_region/WA_median_regions-0000006656-0000006656.tif")#

wa_spec_merge <- terra::merge(wa_spec1, wa_spec2, wa_spec3, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_median_GEE_spec_median_merge.tif", overwrite = T) |>
    project("EPSG:26910")

#### Extraction ####
wa_geo_ext <- read.csv("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_geo_extract.csv")

wa_dat_pts_ext <- wa_dat_pts |> terra::extract(x = geomorph_stack, bind = TRUE, method = "simple") |>
    terra::extract(x = hli_stack, bind = TRUE, method = "simple") |>
    terra::extract(x = wa_precip, bind = TRUE, method = "simple") |> 
    terra::extract(x = wa_temp, bind = TRUE, method = "simple") |>
    terra::extract(x = merge_wip, bind = TRUE, method = "simple") |>
    terra::extract(x = wa_spec_merge, bind = TRUE, method = "simple") |> 
    tidyterra::rename(HLI = lyr.1,
                      Precip = PRISM_ppt_30yr_normal_800mM4_annual_asc,
                      Temp = PRISM_tmean_30yr_normal_800mM4_annual_asc,
                      WIP = WET) |>
    tidyterra::select(-uncertainty) |>
    tidyterra::bind_spat_cols(wa_geo_ext[,c("LITHOLOGY", "GEOLOGIC_AGE")])
(values(wa_dat_pts_ext))
names(wa_dat_pts_ext)

write.csv(wa_dat_pts_ext, file = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/All_WA_horizons_spec_geoagelith.csv")
writeVector(wa_dat_pts_ext, file = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/All_WA_horizons_spec_geoagelith.gpkg", overwrite = T)





