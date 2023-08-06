library(lme4)
library(MASS)
library(lmerTest)
library(RLRsim)
library(terra)
library(sf)
library(tidyverse)
library(car)
library(ggplot2)
library(sjPlot)
library(sjstats)
library(ggeffects)
library(merTools)
library(glmnet)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(cowplot)


setwd('/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')


wa <- read.csv("SOIL CARBON/ANALYSIS/WA_SOC_pts_spec.csv")
wa_pts <- vect(wa, geom = c("x", "y"), crs = "EPSG:26910")
# 
# #wa <- read.csv("SOIL CARBON/ANALYSIS/WA_SOC.csv")
# hoh_pts <- wa_pts[wa_pts$site == "HOH",]
# mas_pts <- wa_pts[wa_pts$site == "MAS",]
# col_pts <- wa_pts[wa_pts$site == "COL",]
# # wa_geom <- geom(wa_pts)
# wa$x <- wa_geom[,'x']
# wa$y <- wa_geom[,'y']
# wa1 <- wa[,c(4:length(names(wa)))]
h_wip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/Hoh_WIP_Mask0_10_2022.tif")
m_wip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mashel_WIP_clip.tif")
c_wip <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/colville_NWI_WIP_clip_INV.tif")
hoh_poly <- vect("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/HOH_POLYGON_7_11_2022/HOH_POLYGON_711.shp")
mas_poly <- vect("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/MAS_poly/mashel_poly.shp")
col_poly <- vect("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/colville_nwi_poly.shp")

wa_geo <- vect("WA_Geo/WA_Geology_100K.gpkg")
wa_geo_df <- as.data.frame(wa_geo)
wa_geo_proj <- terra::project(wa_geo, "EPSG:26910")
wa_pts_geo <- terra::extract(wa_geo_proj, wa_pts)
wa_pts$age_lithol_extract <- wa_pts_geo$AGE_LITHOLOGY
wa_pts$x <- wa$x.x
wa_pts$y <- wa$y
write.csv(wa_pts, "SOIL CARBON/ANALYSIS/WA_SOC_pts_spec.csv")
#wa_geo_proj_lith <- wa_geo_proj["LITHOLOGY"]
#geo_age_lith <- crop(wa_geo_proj, col_poly)
# plot(hoh_geoage, type = "classes", "GEO_AGE_SUM")
#geo_lith_rsp <- rasterize(geo_lith, c_wip, field = "LITHOLOGY", filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/col_lithol.tif")
# # plot(hoh_geoagerast)
# # hoh_sitetif <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_site.tif")
# # writeRaster(hoh_sitetif, "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_sitecat.tif", overwrite = T)
# geo_extract <- terra::extract(wa_geo_proj, wa_pts)
# wa$LITHOL <- geo_extract$LITHOLOGY
# wa <- na.omit(wa)
#write.csv(wa, "SOIL CARBON/ANALYSIS/WA_SOC_geo.csv")
hoh <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/all_spec_GEE/hoh_Spectral_NDYI.tif")
mas <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/all_spec_GEE/mas_Spectral_NDYI.tif")
col <- rast("SOIL CARBON/SPATIAL LAYERS/GEE/all_spec_GEE/col_Spectral_NDYI.tif")

mas_polyrpj <- project(mas_poly, "EPSG:4326")
col_polyrpj <- project(col_poly, "EPSG:4326")

mascrop <- crop(mas, mas_polyrpj)
colcrop <- crop(col, col_polyrpj)
colmask <- mask(colcrop, col_polyrpj)

hohprj <- terra::project(hoh, "EPSG:26910")
masprj <- terra::project(mascrop, m_wip)
colprj <- terra::project(colmask, c_wip)


writeRaster(m_wip_mask, "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mashel_WIP_clip.tif", overwrite = T)
# hoh <- wa[wa$site == "HOH",]
# hoh_pts <- vect(hoh, geom = c("x", "y"), crs = "EPSG:26910")
# hoh_extract <- terra::extract(hoh20prj, hoh_pts, ID=F)
# wa[wa$site == "HOH", "EVI"] <- hoh_extract$EVI






masdtm <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mashel_dem_resamp_clip.tif")
library(whitebox)
library(randomcoloR)
col = randomcoloR::distinctColorPalette(k = 10)
wbt_geomorphons(
    dem = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/colville_dem_resamp.tif",
    output = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/col_geomorphons.tif",
    search = 125,
    threshold = 5,
    fdist = 100,
    skip = 25
)


# hoh_geomorphons <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_geomorphons.tif")
# plot(hoh_geomorphons, type = "classes", col = col)
# hoh_extract <- terra::extract(hoh_geomorphons, hoh_pts, bind = T)
# wa_dat[wa_dat$site == "HOH", "geomorph"] <- hoh_extract$hoh_geomorphons
# 
# mas_geomorphons <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/MAS/mashel_geomorphons.tif")
# plot(mas_geomorphons, type = "classes", col = col)
# plot(mas_pts, add = T)
# mas_extract <- terra::extract(mas_geomorphons, mas_pts, bind = T)
# wa_dat[wa_dat$site == "MAS", "geomorph"] <- mas_extract$mashel_geomorphons
# 
# col_geomorphons <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/COL/col_geomorphons.tif")
# plot(col_geomorphons, type = "classes", col = col)
# col_extract <- terra::extract(col_geomorphons, col_pts, bind = T)
# wa_dat[wa_dat$site == "COL", "geomorph"] <- col_extract$col_geomorphons
# 
# wa$geomorph <- wa_dat$geomorph

#subset for easier viewing and access
wa_dat <- wa |> dplyr::select( "site.x", "sample_name","CHN_30cm_Cstock_Mg_cm2",
                                            "CHN_90cm_Cstock_Mg_ha",  "CHN_1m_Cstock_Mg_ha",   "CHN_120cm_Cstock_Mg_ha",
                                            "WIP.x", "geo_extract.Label",  "hli.x" ) |>
    mutate_if(is.character, as.factor) |>
    stats::setNames(c("site", "name", "SOC30",
                      "SOC90",  "SOC1",   "SOC120",
                      "WIP", "GEO_Label",  "hli")) |>
    as.data.frame()
str(wa_dat)

hist(wa_dat$SOC1)


# Quick R^2 function
r.sq <- function(y,y.fitted){
    res <- y-y.fitted
    1-sum(res^2)/sum((y-mean(y))^2)
}



#### LMER #####
wa_mod1 <- lmer(log(SOC1) ~ WIP + (1|site) + (1|GEO_Label), data = wa_dat)
wa_mod2 <- lmer((SOC1) ~ WIP + (1|GEO_Label), data = wa_dat) 
wa_mod3 <-  lmer((SOC1) ~ WIP + (1|site), data = wa_dat)
wa_mod4 <-  lm((SOC1) ~ WIP, data = wa_dat)

exactRLRT(wa_mod2, wa_mod1, wa_mod3)
exactRLRT(wa_mod3, wa_mod1, wa_mod2)

wa_mod2F <- lmer((SOC1) ~ WIP + (1|GEO_Label), data = wa_dat, REML = F) 
wa_mod3F <-  lmer((SOC1) ~ WIP + (1|site), data = wa_dat, REML = F)

exactLRT(wa_mod2F, wa_mod4)
exactLRT(wa_mod3F, wa_mod4)

#### LMER Summary####

summary(wa_mod1)
plot(wa_mod1)
r.sq(log(wa_dat$SOC1), fitted(wa_mod1))
plot(log(wa_dat$SOC1), fitted(wa_mod1), ylim = c(3, 7), xlim = c(3,7),col = unique(wa_dat$site))
abline(0,1)
# Add a legend
legend("bottomright", pch=16,
       legend = levels(wa_dat$site),
       col = unique(wa_dat$site))

RMSE.merMod(wa_mod1, scale = T)*sd(wa_dat$SOC1)

hist(wa_dat$SOC1)

##### GLMER #####

#wa_gmod1.1 <- glmer( (SOC1) ~   WIP*site + hli+ (1|GEO_Label), family = Gamma(link = "log"), data = wa_dat);summary(wa_gmod1.1)
#wa_gmod1.2 <- glmer( (SOC1) ~   WIP*GEO_Label + hli+ (1|site), family = Gamma(link = "log"), data = wa_dat);summary(wa_gmod1.2)
wa_gmod1.3 <- glmer( (SOC1) ~   WIP + hli+ (1|GEO_Label:site), family = Gamma(link = "log"), data = wa_dat);summary(wa_gmod1.3)
wa_gmod1.4 <- glmer( (SOC1) ~   hli + (1|GEO_Label:site), family = Gamma(link = "log"), data = wa_dat);summary(wa_gmod1.4)
wa_gmod2 <- glmer((SOC1) ~ WIP + (1|GEO_Label:site), family = Gamma(link = "log"),data = wa_dat);summary(wa_gmod2)
wa_gmod3 <-  glmer((SOC1) ~ 1 + (1|GEO_Label:site), family = Gamma(link = "log"),data = wa_dat);summary(wa_gmod3)
wa_gmod4 <-  glm((SOC1) ~ WIP, family = Gamma(link = "log"), data = wa_dat);summary(wa_gmod4)

#exactRLRT(wa_gmod2, wa_gmod1.3, wa_gmod4) #doesn't work for glm?
#exactLRT(wa_gmod1.3, wa_gmod1.4)
gamma.dispersion(wa_gmod4)
anova(wa_gmod1.3, wa_gmod1.4)

dd <- lme4::fortify.merMod(wa_gmod1.3)
ggplot(dd, aes( sample =.resid)) + stat_qq() 

#### spatial autocorrelation ####
wa_pts <- st_as_sf(wa_dat, coords = c("x", "y"), crs = "EPSG:26910")
wa_pts <- na.omit(wa_pts)

spatialEco::correlogram(wa_pts, resid(wa_gmod1.3), dist = 85000)
plot(wa_pts['SOC1'])

#### spatial autocorrelation by study area #####
hoh_pts <- st_as_sf(wa_dat[wa_dat$site == "HOH",], coords = c("x", "y"), crs = "EPSG:26910")
mas_pts <- st_as_sf(wa_dat[wa_dat$site == "MAS",], coords = c("x", "y"), crs = "EPSG:26910")
col_pts <- st_as_sf(wa_dat[wa_dat$site == "COL",], coords = c("x", "y"), crs = "EPSG:26910")

spatialEco::correlogram(hoh_pts, hoh_pts$SOC1, dist = 1000)
spatialEco::correlogram(mas_pts, mas_pts$SOC1, dist = 2000)
spatialEco::correlogram(col_pts, col_pts$SOC1, dist = 800)

#### Gamma Summary ####
summary(wa_gmod1.3)
deviance(wa_gmod1.3) #deviance - measure of fit to the saturated model
            # may need to read Faraway for Gamma deviance

#Dharma package?
#install.packages("DHARMa")
library(DHARMa)
simulationOut <- simulateResiduals(wa_gmod1.3)
plot(simulationOut)

#Model Diagnostics -> full global model 
plot(wa_gmod1.3)
deviance(wa_gmod1.3)
plot(cooks.distance(wa_gmod1.3))
qqPlot(fitted(wa_gmod1.3))

qqnorm(unlist(ranef(wa_gmod1.3)))
qqline(unlist(ranef(wa_gmod1.3)))

pchisq(deviance(wa_gmod1.3),length(wa_dat$site)-1,lower.tail = F) #p-value deviance


#fitted vs. actual
r.sq((wa_dat$SOC1), fitted(wa_gmod1.3))
plot(wa_dat$SOC1, fitted(wa_gmod1.3), col = as.factor(wa_dat$GEO_Label) )
abline(0,1)
abline(lm((fitted(wa_gmod1.3)~ wa_dat$SOC1)), col = "red", lty = 3)
#lines(seq(0,700,7.4), predict(wa_gmod1, newdata = wa_dat, type = "response"), col = "black")
coef(wa_gmod1.3)
# Add a legend
legend("topleft", pch=16,
       legend = levels(wa_dat$GEO_Label),
       col = unique(wa_dat$GEO_Label))

rmse(wa_gmod1.3)#*sd(wa_dat$SOC1)
confint(wa_gmod1.3, method = "boot",nsim = 500)


#confint(wa_gmod1.3)


ggplot(wa_dat, aes(y = fitted(wa_gmod1.3), x = WIP, col = GEO_Label)) +
    geom_point(aes(shape = as.factor(WIP>=0.5) ),size = 6, stroke = 1, show.legend = F)+
    geom_smooth(aes(y = predict(wa_gmod1.3, type = "response")), se = F)


gm <- ggeffects::ggpredict(wa_gmod1, type = "re", terms = c("WIP","GEO_Label"))
plot(gm)




wa_pts <- vect("SOIL CARBON/All_WA/WA_SOC_GEE/WA_SOC_GEE_prj.gpkg")#vect(wa, geom=c("x", "y"), "EPSG:26910")
wa_df <- as.data.frame(wa_pts)
wa_dfn <- dplyr::select_if(wa_df, is.numeric)
wa_dfc <- cor_pmat(wa_dfn)
corrplot(wa_dfn, cutpts=c(0,0.25,0.5,0.75), color=TRUE)
ggcorrplot(wa_dfc)
plot(wa_dfn$CHN_1m_Cst, wa_dfn$NDYI)

gee_m <- glmer(CHN_1m_Cst ~ WIP:EVI + WIP:NDYI + hli + site +  (1|LITHOL:geo_extrac), family = Gamma(link = "log"), data = wa_df)
summary(gee_m)
max(fitted(gee_m))

{plot(wa_df$CHN_1m_Cst, fitted(gee_m),  ylim = c(0,700), xlim = c(0, 700))
abline(0,1)}

#wa_large <- merge(wa_df, wa, by = "y")
#(plot(wa_large$CHN_1m_Cst, wa_large$CHN_1m_Cstock_Mg_ha))
#write.csv(wa_large, "SOIL CARBON/ANALYSIS/WA_SOC_GEE.csv")
