library(terra)
library(aqp)
library(ithir)
library(mpspline2)

setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")

#### Data and spline data ####
wa_dat <- read.csv("SOIL CARBON/ANALYSIS/All_WA_horizons_spec_geoage.csv")
str(wa_dat)

wa_spl_dat_0_30_60_100 <- read.csv("SOIL CARBON/ANALYSIS/R/WA_SOC_controls/wa_spl_dat_0_30_60_100.csv") |>
    mutate(across(where(is.character),as.factor))
colnames(wa_spl_dat_0_30_60_100) <- gsub("_median","",colnames(wa_spl_dat_0_30_60_100))
str(wa_spl_dat_0_30_60_100)

#### Modeling from other Rmd ####
library(lme4)
library(lmerTest)
library(randomForest)

columns_to_exclude_spl <- c("SOC_stock_spline", "lower_depth") 

wa_spl_dat_scale <- wa_spl_dat_0_30_60_100 |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude_spl),
        ~dplyr::case_when(TRUE ~ scale(.))))

full <- wa_spl_dat_0_30_60_100 |>
    dplyr::select(-c(sample_ID, upper_depth, X, site, SCI, HLI, SAVI, EMBI, DSI, DSWI1, LSWI, LITHOLOGY, MNDWI, tree_canopy_cover))

mod2.1 <- lmer(log10(SOC_stock_spline) ~ 
                   WIP+Precip+Temp  + NDYI + GEOLOGIC_AGE +
                   lower_depth + (lower_depth||sample_ID),  
               data = wa_spl_dat_scale, REML = F)

dmod <- lmer(log10(SOC_stock_spline) ~ lower_depth + NDYI + Precip + Temp +  
                 WIP + (1 | sample_ID) + (0 + lower_depth | sample_ID) + NDYI:Precip +  
                 NDYI:Temp + NDYI:WIP + Precip:Temp + Precip:WIP + Temp:WIP +      NDYI:Precip:Temp + NDYI:Temp:WIP,
             data = wa_spl_dat_scale, REML = F)

rf_model <- randomForest((SOC_stock_spline) ~ .,
                         ntree = 1000, mtry = 4,
                         importance = TRUE, data = full)
#### spatial data ####

hoh_wip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/Hoh_WIP_Mask0_10_2022.tif")
test_box <- as.polygons(ext(c(411521.887913524, 420933.712237848, 5292301.97939307, 5298455.8645282)))
test_wip <- terra::crop(hoh_wip, test_box) |> terra::scale(center = 0.4429349, scale = 0.2924148)

precip <- rast("PRISM_Climate/PRISM_ppt_30yr_normal_800mM4_annual_asc/PRISM_ppt_30yr_normal_800mM4_annual_asc.asc") |> terra::project("EPSG:26910")
test_precip <- precip |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = 1870.354, scale = 1022.088)

temp <- rast("PRISM_Climate/PRISM_tmean_30yr_normal_800mM4_annual_asc/PRISM_tmean_30yr_normal_800mM4_annual_asc.asc") |> terra::project("EPSG:26910")
test_temp <- temp |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = 8.42412, scale = 1.616208)

ndyi <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/all_spec_GEE/hoh_Spectral_NDYI.tif")
test_ndyi <- ndyi |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = 0.4192793, scale = 0.08742682)

spec <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/All_WA/WA_GEE_spec_median.tif") |> terra::project("EPSG:26910")
test_ndvi <- spec["NDVI"] |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = 0.8315177, scale = 0.08425498)
test_andwi <- spec["ANDWI"] |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = -0.7033847, scale = 0.05644566)
test_evi <- spec["EVI"] |> terra::crop(test_box) |> terra::resample(test_wip) |> terra::scale(center = 0.4798586, scale = 0.09722779)

wa_geology <- vect("WA_Geo/ger_portal_surface_geology_100k/surface_geology_100k.gdb", layer = "geologic_unit_poly_100k") |> 
    terra::project("EPSG:26910") |> tidyterra::select(LITHOLOGY, GEOLOGIC_AGE)
test_geo <- wa_geology |> terra::crop(test_box) |> terra::rasterize(test_wip, field = "GEOLOGIC_AGE") |> 
    terra::subst("Present", "Quaternary")

catconv <- data.frame(id = 1:10, geomorphons = letters[1:10])
hoh_geomorph <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_geomorphons.tif")
test_geomorph <- terra::crop(hoh_geomorph, test_box) |> 
    categories(layer = 1, value = catconv) |> project("EPSG:26910")#|> terra::resample(hoh_geomorph)

m <- c(0, 1, 30)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
test_depth <- classify(test_wip, rclmat, include.lowest = T)

pred_stack <- c(test_wip, test_precip, test_temp, test_ndyi, test_geo)
names(pred_stack) <- c("WIP", "Precip", "Temp", "NDYI", "GEOLOGIC_AGE")
plot(pred_stack)

map_func <- function(stack){
    stack_names <- names(stack)
    print(stack_names)
    wipcolr <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))
    preccolr <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
    tempcolr <- colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
    ndyicolr <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrBr"))
    geocolr <- colorRampPalette(RColorBrewer::brewer.pal(4, "Accent"))
    
    wip_plot<- levelplot(stack[stack_names[1]], margin = F, par.settings=list(
        axis.line=list(col='transparent') # suppress axes and legend outline
    ), scales=list(draw=FALSE), maxpixels = 9e4, col.regions = wipcolr,
    main = "WIP")
    
    prec_plot <- levelplot(stack[stack_names[2]], margin = F, par.settings=list(
        axis.line=list(col='transparent') # suppress axes and legend outline
    ), scales=list(draw=FALSE), maxpixels = 9e3, col.regions = preccolr,
    main = "Precipitation")
    
    temp_plot<- levelplot(stack[stack_names[3]], margin = F, par.settings=list(
        axis.line=list(col='transparent') # suppress axes and legend outline
    ), scales=list(draw=FALSE), maxpixels = 9e3, col.regions = tempcolr,
    main = "Temperature")
    
    ndyi_plot <- levelplot(stack[stack_names[4]], margin = F, par.settings=list(
        axis.line=list(col='transparent') # suppress axes and legend outline
    ), scales=list(draw=FALSE), maxpixels = 9e4, col.regions = ndyicolr,
    main = "NDYI")
    
    geo_plot <- levelplot(stack[stack_names[5]], margin = F, par.settings=list(
        axis.line=list(col='transparent') # suppress axes and legend outline
    ), scales=list(draw=FALSE), maxpixels = 9e3, col.regions =geocolr,
    main = "Geologic Age")
    
    top_plots<- ggpubr::ggarrange(wip_plot, prec_plot, temp_plot, ndyi_plot, ncol = 2, nrow = 2, align = "v")
    all_plots <- ggpubr::ggarrange(top_plots, geo_plot, ncol = 1, nrow = 2, heights = c(2,1))
    return(all_plots)
}


map_func(pred_stack)

#### spatial prediction ####

test_depth30 <- data.frame(lower_depth = 30)
test_depth60 <- data.frame(lower_depth = 60)
test_depth100 <- data.frame(lower_depth = 100)

mod2.1_predict30 <- predict(pred_stack, mod2.1, re.form = NA, allow.new.levels = TRUE, const = test_depth30)
writeRaster(10**mod2.1_predict30, filename = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict30.tif",
            overwrite = T)
                            
mod2.1_predict60 <- predict(pred_stack, mod2.1, re.form = NA, allow.new.levels = TRUE, const = test_depth60)
writeRaster(10**mod2.1_predict60, filename = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict60.tif",
                            overwrite = T)
mod2.1_predict100 <- predict(pred_stack, mod2.1, re.form = NA, allow.new.levels = TRUE, const = test_depth100)
writeRaster(10**mod2.1_predict100, filename = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict100.tif",
                             overwrite = T)
plot(10**(mod2.1_predict30))
plot(10**(mod2.1_predict60))
plot(10**(mod2.1_predict100))

soc_sum <- 10**mod2.1_predict30 + 10**mod2.1_predict60 + 10**mod2.1_predict100
plot(soc_sum)


pred_stack_rf <- c(test_wip, test_precip, test_temp, test_ndyi, test_ndvi, test_evi, test_andwi, test_geo, test_geomorph)
names(pred_stack_rf) <- c("WIP", "Precip", "Temp", "NDYI", "NDVI", "EVI", "ANDWI", "GEOLOGIC_AGE", "geomorphons")
plot(pred_stack_rf)
rf_predic30 <- predict(pred_stack_rf, rf_model, re.form = NA, allow.new.levels = TRUE, const = test_depth30)
