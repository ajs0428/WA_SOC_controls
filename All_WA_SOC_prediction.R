## ----setup, include=FALSE------------------------------------------------------------------------
library(formatR)
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.show = "hold", time_it = TRUE, dpi = 100)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = T, collapse = TRUE)
knitr::opts_knit$set(root.dir = '/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')
library(rgl)
library(terra)
library(lme4)
library(MASS)
library(lmerTest)
library(MuMIn)
library(terra)
library(sf)
library(ggplot2)
library(ggeffects)
library(merTools)
library(glmnet)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(webshot)
library(kableExtra)
library(formatR)
library(dplyr)
library(stringr)


knitr::knit_hooks$set(webgl = hook_webgl)
rgl::setupKnitr(autoprint = TRUE)

terraOptions(
    memfrac = 0.1
)

setGDALconfig("GDAL_PAM_ENABLED", "FALSE")


## ------------------------------------------------------------------------------------------------
pred_path <- "SOIL CARBON/All_WA/data/Rasters/PredictorStacks/"


## ------------------------------------------------------------------------------------------------
hoh_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/hoh_pts_2855.gpkg"))
mas_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/mas_pts_2856.gpkg"))
col_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/col_pts_2855.gpkg"))

wa_dat <- rbind(hoh_dat, mas_dat, col_dat) |> 
  mutate(
  GEO = as.factor(GEO),
  geomorphons = as.factor(geomorphons),
  site = as.factor(site),
  site = forcats::fct_reorder(site, SOC_stock_spline, .fun = "median")) |>
  dplyr::rename_with(~gsub("_median", "", .x, fixed = TRUE))

columns_to_exclude <- c("SOC_stock_spline") 

wa_dat_scale <- wa_dat |> 
    dplyr::select(sample_ID, lower_depth, SOC_stock_spline, 
                  site, DTM, GEO, WIP,HLI, CHM) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
                  ~dplyr::case_when(TRUE ~ scale(.))),
        site = forcats::fct_reorder(site, SOC_stock_spline, .fun = "median")) 




## ------------------------------------------------------------------------------------------------
wa_dat_scale_params <- wa_dat |> 
    dplyr::select(sample_ID, lower_depth, SOC_stock_spline, 
                  site, DTM, GEO, WIP,HLI, CHM) |>
    dplyr::summarise(across(dplyr::where(is.numeric) & !all_of(columns_to_exclude),
                            list(mean = mean, sd = sd)))


## ------------------------------------------------------------------------------------------------
hoh_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Hoh_PredictorStack_Class.tif") |>
    terra::subset(c("DTM", "GEO", "HLI", "WIP", "CHM"))
mas_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Mas_PredictorStack_Class.tif") |>
    terra::subset(c("DTM", "GEO", "HLI", "WIP", "CHM"))
col_stack <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/Col_PredictorStack_Class.tif") |>
    terra::subset(c("DTM", "GEO", "HLI", "WIP", "CHM"))


## ------------------------------------------------------------------------------------------------
#| eval: false

## rast_scale_func <- function(stack, path) {
##     stack_scale <- stack |>
##         tidyterra::mutate(DTM = (DTM - wa_dat_scale_params$DTM_mean)/wa_dat_scale_params$DTM_sd,
##                           HLI = (HLI - wa_dat_scale_params$HLI_mean)/wa_dat_scale_params$HLI_sd,
##                           WIP = (WIP - wa_dat_scale_params$WIP_mean)/wa_dat_scale_params$WIP_sd,
##                           CHM = (CHM - wa_dat_scale_params$CHM_mean)/wa_dat_scale_params$CHM_sd)
##     writeRaster(stack_scale, paste0(getwd(), "/", path, deparse(substitute(stack)), "_scale", ".tif"),
##                 overwrite = TRUE)
##     return(stack_scale)
## }
## 
## hoh_stack_scale <- rast_scale_func(hoh_stack, pred_path)
## mas_stack_scale <- rast_scale_func(mas_stack, pred_path)
## col_stack_scale <- rast_scale_func(col_stack, pred_path)
## 
## 
## 


## ------------------------------------------------------------------------------------------------
hoh_stack_scale <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/hoh_stack_scale.tif")
mas_stack_scale <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/mas_stack_scale.tif")
col_stack_scale <- rast("SOIL CARBON/All_WA/data/Rasters/PredictorStacks/col_stack_scale.tif")


## ------------------------------------------------------------------------------------------------
hoh_consolidate_geo <- vect("SOIL CARBON/All_WA/data/Vectors/Hoh_geology_consolidate.shp")
mas_consolidate_geo <- vect("SOIL CARBON/All_WA/data/Vectors/Mas_geology_consolidate.shp")
col_consolidate_geo <- vect("SOIL CARBON/All_WA/data/Vectors/Col_geology_consolidate.shp") 
plot(hoh_consolidate_geo, "GEO")
plot(hoh_stack_scale$GEO)


## ------------------------------------------------------------------------------------------------
hoh_stack_scale$GEO <- subst(hoh_stack_scale$GEO, from = c(1,2,3), 
                       c("MioceneEocene", "Pleistocene", "Quaternary"))
mas_stack_scale$GEO <- subst(mas_stack_scale$GEO, from = c(1,2,3), 
                       c("MioceneEocene", "Pleistocene", "Quaternary"))
col_stack_scale$GEO <- subst(col_stack_scale$GEO, from = c(1,2,3,4), 
                       c("MioceneEocene", "Pleistocene", "PreTertiary", "Quaternary"))

plot(hoh_stack_scale$GEO)
plot(mas_stack_scale$GEO)
plot(col_stack_scale$GEO)


## ------------------------------------------------------------------------------------------------
model <- get(load("SOIL CARBON/All_WA/analysis/models/All_WA_Model7_Spline.RData"))
formula(model)


## ------------------------------------------------------------------------------------------------
SOC_pred_path <- "SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/"

SOC_pred_func <- function(stack, depths, model, site, path){
    for(i in 1:length(depths)){
      
      const_df_15 <- data.frame(site = site, 
                                lower_depth = ((15 - wa_dat_scale_params$lower_depth_mean)/wa_dat_scale_params$lower_depth_sd), 
                                stringsAsFactors = TRUE)
      const_df_30 <- data.frame(site = site, 
                                lower_depth = ((30 - wa_dat_scale_params$lower_depth_mean)/wa_dat_scale_params$lower_depth_sd), 
                                stringsAsFactors = TRUE)
      const_df_60 <- data.frame(site = site, 
                                lower_depth = ((60 - wa_dat_scale_params$lower_depth_mean)/wa_dat_scale_params$lower_depth_sd), 
                                stringsAsFactors = TRUE)
      const_df_100 <- data.frame(site = site, 
                                 lower_depth = ((100 - wa_dat_scale_params$lower_depth_mean)/wa_dat_scale_params$lower_depth_sd), 
                                 stringsAsFactors = TRUE)
      const_df_200 <- data.frame(site = site, 
                                 lower_depth = ((200 - wa_dat_scale_params$lower_depth_mean)/wa_dat_scale_params$lower_depth_sd), 
                                 stringsAsFactors = TRUE)
      
      const_df_list <- list(const_df_15, const_df_30, const_df_60, const_df_100, const_df_200)
      depth_num_list <- list(15, 30, 60, 100, 200)
      
          predict(stack, model, #se.fit = TRUE,
                             const = const_df_list[[i]],
                             na.rm = TRUE, 
                             re.form = NA, allow.new.levels = TRUE,
                             cores = (parallel::detectCores()-2),
                             filename = paste0(getwd(), 
                                               "/", 
                                               path, 
                                               stringr::str_remove(deparse(substitute(stack)), "stack"),
                                               "SOCpredict_",
                                               depth_num_list[[i]], 
                                               ".tif"),
                             overwrite = TRUE)
    }
}


## ------------------------------------------------------------------------------------------------
testpt <- terra::buffer(spatSample(mas_stack_scale, 1, as.points = TRUE, na.rm = TRUE), 200)

test <- terra::crop(mas_stack_scale, testpt)


## ------------------------------------------------------------------------------------------------
test_pred <- SOC_pred_func(test, model = model, depths = 3, site = "MAS", path = SOC_pred_path)

#plot(rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/testSOCpredict_60.tif"))
#plot(test)


## ------------------------------------------------------------------------------------------------
#| eval: false
## SOC_pred_func(hoh_stack_scale, model = model, depths = 5, site = "HOH", path = SOC_pred_path)
## SOC_pred_func(mas_stack_scale, model = model, depths = 5, site = "MAS", path = SOC_pred_path)
## SOC_pred_func(col_stack_scale, model = model, depths = 5, site = "COL", path = SOC_pred_path)

## ------------------------------------------------------------------------------------------------
plot(10**rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/col__scaleSOCpredict_15.tif"))


## ------------------------------------------------------------------------------------------------

SOC_sum_units <- function(site, depth, path){
  stopifnot(is.character(depth))
  depth_list <- list(15, 30, 60, 100, 200)
  stopifnot(as.numeric(depth) %in% depth_list)
  
  depth_remove <- if(as.numeric(depth) < 200){
    depth_list[depth_list>as.numeric(depth)]
  } else {
    1000
  }
  
  filelist <- list.files(path, pattern = tolower(site), full.names = TRUE)
  filelist_depth <- filelist[!str_detect(filelist,pattern=as.character(depth_remove))]
  
  SOC <- rast(filelist_depth) |> 
      terra::app(fun=function(i) 10**(i)) |> sum()*100
  writeRaster(SOC, filename = paste0(getwd(), 
                                     "/", 
                                     path, 
                                     tolower(site), 
                                     "SOCsum",
                                     depth,
                                     ".tif"), overwrite = TRUE)
    
}


## ------------------------------------------------------------------------------------------------
#| eval: false
## hoh_SOC_sum200 <- SOC_sum_units(site = "hoh__", depth = "200", path = SOC_pred_path)
## mas_SOC_sum200 <- SOC_sum_units(site = "mas__", depth = "200", path = SOC_pred_path)
## col_SOC_sum200 <- SOC_sum_units(site = "col__", depth = "200", path = SOC_pred_path)
## 
## hoh_SOC_sum100 <- SOC_sum_units(site = "hoh__", depth = "100", path = SOC_pred_path)
## mas_SOC_sum100 <- SOC_sum_units(site = "mas__", depth = "100", path = SOC_pred_path)
## col_SOC_sum100 <- SOC_sum_units(site = "col__", depth = "100", path = SOC_pred_path)


## ------------------------------------------------------------------------------------------------
  test <- if((tools::toTitleCase(str_extract(deparse(substitute(col_stack$WIP)), "hoh|mas|col")) == "Hoh")){
    ("Hoh_WIP_Final_2855_Mask.tif")
    } else {
    ("False")
    }
test


## ------------------------------------------------------------------------------------------------
#| eval: false
## hoh_mndwi <- rast("SOIL CARBON/All_WA/data/Rasters/Hoh_spectral.tif", lyr = "MNDWI_median") |>
##   resample(y = hoh_stack_scale$WIP)
## mas_mndwi <- rast("SOIL CARBON/All_WA/data/Rasters/Mashel_spectral.tif", lyr = "MNDWI_median") |>
##   resample(y = mas_stack_scale$WIP)
## col_mndwi <- rast("SOIL CARBON/All_WA/data/Rasters/Colville_spectral.tif", lyr = "MNDWI_median") |>
##   resample(y = col_stack_scale$WIP)
## 
## SurfaceWater_Mask <- function(mndwi, lyr, path) {
## 
##   if(str_detect(deparse(substitute(lyr)), "WIP") == TRUE) {
##      wipname <-  if((tools::toTitleCase(str_extract(deparse(substitute(lyr)),
##                                                     "hoh|mas|col")) == "Hoh")){
##        "Hoh_WIP_Final_2855_Mask.tif"
##        } else if((tools::toTitleCase(str_extract(deparse(substitute(lyr)),
##                                                  "hoh|mas|col")) == "Mas")){
##      "Mashel_WIP_Final_2856_Mask.tif"
##        } else if((tools::toTitleCase(str_extract(deparse(substitute(lyr)),
##                                                 "hoh|mas|col")) == "Col")){
##       "Colvile_WIP_Final_2855_Mask.tif"
##        } else {
##        break
##        }
##     masked <- mask(lyr, mndwi > 0.0 | is.na(mndwi), maskvalues = TRUE,
##                            filename = paste0(paste0(getwd(),
##                                      "/",
##                                      "SOIL CARBON/All_WA/data/Rasters/",
##                                      wipname)),
##                            overwrite = TRUE)
##   } else {
##     masked <- mask(lyr, mndwi > 0.0 | is.na(mndwi), maskvalues = TRUE,
##                            filename = paste0(paste0(getwd(),
##                                      "/",
##                                      path,
##                                      str_remove_all(deparse(substitute(lyr)),
##                                                     pattern = "_"),
##                                      "mask",
##                                      ".tif")),
##                            overwrite = TRUE)
##   }
## }
## 
## hoh_SOC_sum100mask <- SurfaceWater_Mask(hoh_mndwi, hoh_SOC_sum100, path = SOC_pred_path)
## mas_SOC_sum100mask <-SurfaceWater_Mask(mas_mndwi, mas_SOC_sum100, path = SOC_pred_path)
## col_SOC_sum100mask <- SurfaceWater_Mask(col_mndwi, col_SOC_sum100, path = SOC_pred_path)
## 
## hoh_WIP_mask <- SurfaceWater_Mask(hoh_mndwi, hoh_stack$WIP, path = SOC_pred_path)
## mas_WIP_mask <-SurfaceWater_Mask(mas_mndwi, mas_stack$WIP, path = SOC_pred_path)
## col_WIP_mask <- SurfaceWater_Mask(col_mndwi, col_stack$WIP, path = SOC_pred_path)


## ------------------------------------------------------------------------------------------------
hoh_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/hohSOCsum100mask.tif")
mas_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/masSOCsum100mask.tif")
col_SOC_sum100mask <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/colSOCsum100mask.tif")

hoh_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Hoh_WIP_Final_2855_Mask.tif", lyrs = "WIP")
mas_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Mashel_WIP_Final_2856_Mask.tif", lyrs = "WIP")
col_wip <- rast("SOIL CARBON/All_WA/data/Rasters/Colville_WIP_Final_2855_Mask.tif", lyrs = "WIP")

plot(hoh_SOC_sum100mask)
plot(mas_SOC_sum100mask)
plot(col_SOC_sum100mask)


## ------------------------------------------------------------------------------------------------
C_map_simp <- function(C_map){
    name <- deparse(substitute(C_map))
    gt <- (C_map > -999)
    cell_size <- cellSize(gt, unit = "ha") |> mask(mask = gt)
    carbon_cell <- C_map*cellSize(gt, unit = "ha") # carbon in Mg per cell which is then added up 
    
    area_tot <- sum(values(cell_size), na.rm = T)
    C_mean <- mean(values(C_map), na.rm = T) #mean value of all values of Mg/ha cells
    TotalC_sum <- sum(values(carbon_cell), na.rm =T)
    
    
    return(data.frame("Name" = name, 
                      "Total_area" = area_tot, 
                      "AverageSOC_Mgha" = C_mean, 
                      "total_Carbon_Tg" = TotalC_sum/1e6,
                      stringsAsFactors = T))
}

C_map_wet_fractions <- function(C_map, WIP){
    name <- deparse(substitute(C_map))
    gt <- (C_map > -999)
    WIP_wetupl <- (WIP >= 0.50)
    WIP_mid <- (WIP >= 0.25 & WIP <= 0.75)
    
    C_map_wet <- mask(C_map, mask = WIP_wetupl, maskvalues = FALSE)
    C_map_upl <- mask(C_map, mask = WIP_wetupl, maskvalues = TRUE)
    C_map_mid <- mask(C_map, mask = WIP_mid, maskvalues = FALSE)
    
    cell_size_all <- cellSize(gt, unit = "ha") |> mask(mask = gt)
    cell_size_wet <- cellSize(gt, unit = "ha") |> mask(mask = C_map_wet)
    cell_size_upl <- cellSize(gt, unit = "ha") |> mask(mask = C_map_upl)
    cell_size_mid <- cellSize(gt, unit = "ha") |> mask(mask = C_map_mid)
    
    carbon_cell_all <- C_map*cell_size_all # carbon in Mg per cell which is then added up 
    carbon_cell_wet <- C_map_wet*cell_size_wet # carbon in Mg per cell which is then added up 
    carbon_cell_upl <- C_map_upl*cell_size_upl # carbon in Mg per cell which is then added up 
    carbon_cell_mid <- C_map_mid*cell_size_mid # carbon in Mg per cell which is then added up 
    
    area_tot <- sum(values(cell_size_all), na.rm = T)
    area_wet <- sum(values(cell_size_wet), na.rm = T)
    area_upl <- sum(values(cell_size_upl), na.rm = T)
    area_mid <- sum(values(cell_size_mid), na.rm = T)
    
    C_mean_all <- mean(values(C_map), na.rm = T) #mean value of all values of Mg/ha cells
    C_mean_wet <- mean(values(C_map_wet), na.rm = T)
    C_mean_upl <- mean(values(C_map_upl), na.rm = T)
    C_mean_mid <- mean(values(C_map_mid), na.rm = T)
    
    TotalC_sum <- sum(values(carbon_cell_all), na.rm =T)
    TotalC_sum_wet <- sum(values(carbon_cell_wet), na.rm =T)
    TotalC_sum_upl <- sum(values(carbon_cell_upl), na.rm =T)
    TotalC_sum_mid <- sum(values(carbon_cell_mid), na.rm =T)
    
    return(data.frame("Name" = c(name, paste0(name, "wet"), paste0(name, "upl"), paste0(name, "mid")), 
                      "Total_area" = c(area_tot, area_wet, area_upl, area_mid),
                      "AverageSOC_Mgha" = c(C_mean_all, C_mean_wet, C_mean_upl, C_mean_mid), 
                      "total_Carbon_Tg" = c(TotalC_sum/1e6, TotalC_sum_wet/1e6, 
                                            TotalC_sum_upl/1e6, TotalC_sum_mid/1e6),
                      stringsAsFactors = T))
}


## ------------------------------------------------------------------------------------------------
#| eval: false
## hoh_soc100_df <- C_map_wet_fractions(hoh_SOC_sum100mask, hoh_wip)
## mas_soc100_df <- C_map_wet_fractions(mas_SOC_sum100mask, mas_wip)
## col_soc100_df <- C_map_wet_fractions(col_SOC_sum100mask, col_wip)
## 
## all_df <- rbind(hoh_soc100_df, mas_soc100_df, col_soc100_df)
## 
## readr::write_csv(all_df, "SOIL CARBON/All_WA/data/dataframes/All_WA_MapSOC_wetuplmid_100mask.csv")


## ------------------------------------------------------------------------------------------------
#| eval: false
## SOC_df <- data.frame(Name = character(), Total_area = double(), AverageSOC_Mgha = double(), total_Carbon_Tg = double())
## 
## Map_SOC_func <- function(path, sumDepth){
##   filelist <-list.files(path, pattern = sumDepth, full.names = TRUE)
##   soc_list <- list()
##   for(i in 1:length(filelist)){
##     r <- rast(filelist[[i]])
##     soc_list[[i]] <-  C_map_simp(r)
##     #plot(r)
##   }
##   SOC_df <- rbind(soc_list)
##   return(SOC_df)
## }
## 
## testdf <- Map_SOC_func(SOC_pred_path, "sum100mask")
## 


## ------------------------------------------------------------------------------------------------
#| eval: false
## hoh_100soc <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/hohSOCsum100mask.tif")
## mas_100soc <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/masSOCsum100mask.tif")
## col_100soc <- rast("SOIL CARBON/All_WA/data/Rasters/SOC_Predictions/colSOCsum100mask.tif")
## 
## SOC_df <- rbind(C_map_simp(hoh_100soc),
## C_map_simp(mas_100soc),
## C_map_simp(col_100soc))
## 
## readr::write_csv(SOC_df, "SOIL CARBON/All_WA/data/dataframes/All_WA_MapSOC100mask.csv")


## ------------------------------------------------------------------------------------------------
#| eval: false
## SOC_df <- readr::read_csv("SOIL CARBON/All_WA/data/dataframes/All_WA_MapSOC_wetuplmid_100mask.csv")
## SOC_df |> mutate(Name = case_match( Name,
##              "hoh_100soc" ~ "Hoh",
##              "mas_100soc" ~ "Mashel",
##              "col_100soc" ~ "Colville"),
##              Total_area = signif(Total_area, digits = 5),
##              AverageSOC_Mgha = signif(AverageSOC_Mgha, digits = 3),
##              total_Carbon_Tg = signif(total_Carbon_Tg, digits = 2)) |>
##         rename("Study Area" = Name,
##                  "Total Area (ha)" = Total_area,
##                  "Mean SOC Stock (Mg ha^-1)" = AverageSOC_Mgha,
##                  "Total SOC Stock (Tg)" = total_Carbon_Tg)
## 
## 


## ------------------------------------------------------------------------------------------------
#| eval: false
## set.seed(11)
## hoh_stack_sample <- spatSample(hoh_stack, size = 1000, na.rm = TRUE)
## 
## kmod <- kmeans(hoh_stack_sample, centers = 5)
## summary(kmod)
## kmap <- terra::k_means(hoh_stack, centers = 5, maxcell = 1e5)
## 
## plot(kmap)

