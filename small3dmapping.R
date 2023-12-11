library(rayshader)
library(terra)
library(MetBrewer)
library(ggplot2)
library(rasterVis)
library(RColorBrewer)
library(magick)
library(rgl)
setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")


#### 2D carbon Plotting function ####
TwoD_plot <- function(layer, col){
    gplot(layer, maxpixels=200000) +geom_tile(aes(fill = value)) +
        col +
        theme(legend.position="bottom", legend.box = "horizontal",
              legend.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "white"))
}


allC_col <- scale_fill_gradientn("Soil Carbon Stock Mg/ha", 
                                 colors = met.brewer('Demuth', type = "continuous", direction = -1)[3:10], 
                                 na.value = "white", 
                                 limits = c(0, 1.5))
   
test_point <- vect(cbind(416141.650431608, 5295489.37534548), crs = "EPSG:26910")
test_circle <- terra::buffer(test_point, 2500)
DEM <- rast("SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/hoh_dtm_resamp.tif")
DEM_small <- crop(DEM, test_circle, mask = T, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/small_test/DEM_small.tif", overwrite = T)



Dmat <- DEM_small |> raster_to_matrix()
Dambmat <- ambient_shade(Dmat, maxsearch = 100, zscale = 5)
Draymat <- ray_shade(Dmat, zscale = 5, lambert = TRUE)

c30 <- rast("SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict30.tif")
c60 <- rast("SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict60.tif")
c100 <- rast("SOIL CARBON/ANALYSIS/R/WA_SOC_controls/mod2.1_predict100.tif")


c30mat <- c30 |> crop(test_circle, mask = T, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/small_test/c30_small.tif") |> raster_to_matrix()
c60mat <- c60 |> crop(test_circle, mask = T, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/small_test/c60_small.tif") |> raster_to_matrix()
c100mat <- c100 |> crop(test_circle, mask = T, filename = "SOIL CARBON/SPATIAL LAYERS/SPATIAL_LAYERS_7_11_22/HOH/small_test/c100_small.tif") |> raster_to_matrix()


Dmat_30 <- Dmat %>%
    height_shade(texture = (grDevices::gray.colors(256))) %>% 
    add_overlay(height_shade(c30mat, range = c(0,1.5), texture = MetBrewer::met.brewer("Demuth", type = "continuous", direction = -1)), 
                alphalayer = 1) %>%
    #add_overlay(height_shade(Dmat, texture = (grDevices::gray.colors(256)[20:205])), alphalayer = 0.5) %>%
    add_shadow(Draymat, max_darken = 0.3) %>%
    add_shadow(Dambmat, max_darken = 0.3) #%>%

Dmat_60 <- Dmat %>%
    height_shade(texture = (grDevices::gray.colors(256))) %>% 
    add_overlay(height_shade(c60mat, range = c(0,1.5), texture = MetBrewer::met.brewer("Demuth", type = "continuous", direction = -1)), 
                alphalayer = 1) %>%
    #add_overlay(height_shade(Dmat, texture = (grDevices::gray.colors(256)[20:205])), alphalayer = 0.5) %>%
    add_shadow(Draymat, max_darken = 0.3) %>%
    add_shadow(Dambmat, max_darken = 0.3)

Dmat_100 <- Dmat %>%
    height_shade(texture = (grDevices::gray.colors(256))) %>% 
    add_overlay(height_shade(c100mat, range = c(0,1.5), texture = MetBrewer::met.brewer("Demuth", type = "continuous", direction = -1)), 
                alphalayer = 1) %>%
    #add_overlay(height_shade(Dmat, texture = (grDevices::gray.colors(256)[20:205])), alphalayer = 0.5) %>%
    add_shadow(Draymat, max_darken = 0.3) %>%
    add_shadow(Dambmat, max_darken = 0.3)



plot_3d(Dmat_30, Dmat, zscale = 7, shadowdepth = -50, solid = F,  
        shadow = F, linewidth = 0.5, 
        zoom=.75, phi=30,theta=315,fov=5,
        baseshape = "circle", background = "black") 



render_floating_overlay(overlay = Dmat_60, 
                        heightmap = Dmat, 
                        altitude = -4000,
                        #baseshape = "circle",
                        zscale = 10, remove_na = F, 
                        clear_layers = F, 
                        horizontal_offset = c(-800, -800),
                        reorient = F)
render_floating_overlay(overlay = Dmat_100, 
                        heightmap = Dmat,
                        altitude = -8000, 
                        #baseshape = "circle",
                        zscale = 10, remove_na = F, 
                        clear_layers = F, 
                        horizontal_offset = c(-1600, -1600),
                        reorient = F)


filename_movie = ('/Users/Anthony/OneDrive - UW/University of Washington/Presentations/Conferences/AGU2023/ExploringControls/Figures and Tables/3dmap_gif.gif')
render_movie(filename = filename_movie)

########## ggplot ##########


