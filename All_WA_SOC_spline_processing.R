library(terra)
library(aqp)
library(ithir)
library(mpspline2)
library(dplyr)

setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")

wa_dat <- read.csv("SOIL CARBON/ANALYSIS/R/WA_SOC_controls/All_WA_horizons_spec_geoagelith.csv")
str(wa_dat)

wa_dat_horC <- wa_dat |>
    dplyr::group_by(sample_ID) |> 
    mutate(top = case_when(is.na(depth_cm - lag(depth_cm)) ~ 0,
                           .default = lag(depth_cm)),
           bottom = depth_cm,
           center = abs(top - (top - bottom)/2)) |> 
    dplyr::select(sample_ID, top, center, bottom, carbon_perc, carbon_stock_g_cm2) 
str(wa_dat_horC)

#create a subset of the data with only the predictor variables
    # this will be joined to the new spline dataframes
    # should be 96 rows (number of pedons collected in WA)
wa_hor_dat_sub <- wa_dat |> dplyr::select(-c(X, depth_cm, BD_g_cm3, carbon_perc, carbon_stock_g_cm2)) |>
    group_by(sample_ID) |>
    filter(row_number() == 1) |> 
    dplyr::select(sample_ID, everything()) |>
    dplyr::arrange(site)
str(wa_hor_dat_sub)


wa_spl <- mpspline_tidy(obj = wa_dat_horC, var_name = "carbon_stock_g_cm2", d = c(0, 30, 60, 100)) 
wa_spl_dat_0_30_60_100 <- wa_spl$est_dcm |>
    left_join(y = wa_hor_dat_sub, by = join_by(sample_ID), relationship = "many-to-one") |>
    rename(SOC_stock_spline = SPLINED_VALUE,
           upper_depth = UD,
           lower_depth = LD)
str(wa_spl_dat_0_30_60_100)
write.csv(wa_spl_dat_0_30_60_100, file = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/wa_spl_dat_0_30_60_100.csv")

wa_spl_30_100 <- mpspline_tidy(obj = wa_dat_horC, var_name = "carbon_stock_g_cm2", d = c(0, 30, 100)) 
wa_spl_dat_30_100 <- wa_spl_30_100$est_dcm |> 
    left_join(y = wa_hor_dat_sub, by = join_by(sample_ID), relationship = "many-to-one") |>
    rename(SOC_stock_spline = SPLINED_VALUE,
           upper_depth = UD,
           lower_depth = LD)
str(wa_spl_dat_30_100)
write.csv(wa_spl_dat_30_100, file = "SOIL CARBON/ANALYSIS/R/WA_SOC_controls/wa_spl_dat_30_100.csv")
