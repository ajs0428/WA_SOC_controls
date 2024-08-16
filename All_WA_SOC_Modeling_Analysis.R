## pre {
##   max-height: 300px;
##   overflow-y: auto;
## }
## 
## h1, h3, h4 {
##   text-align: center;
## }

## ----setup, include=FALSE------------------------------------------------------------------------------
library(formatR)
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.show = "hold", time_it = TRUE, dpi = 100)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 60), tidy = T, collapse = TRUE)
knitr::opts_knit$set(root.dir = '/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/')
library(rgl)
library(terra)
library(lme4)
library(MASS)
library(mgcv)
library(lmerTest)
library(MuMIn)
library(RLRsim)
library(terra)
library(spatialEco)
library(sf)
library(mapview)
library(car)
library(ggplot2)
library(forcats)
library(sjPlot)
library(sjstats)
library(DHARMa)
library(ggeffects)
library(merTools)
library(glmnet)
library(stats)
library(ggcorrplot)
library(RColorBrewer)
library(cowplot)
library(webshot)
library(kableExtra)
library(pdp)
library(vip)
library(formatR)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)

knitr::knit_hooks$set(webgl = hook_webgl)
rgl::setupKnitr(autoprint = TRUE)


## ----include=FALSE-------------------------------------------------------------------------------------

# Quick R^2 function
r.sq <- function(y,y.fitted){
    res <- y-y.fitted
    1-sum(res^2)/sum((y-mean(y))^2)
}


## ----echo=FALSE----------------------------------------------------------------------------------------
hoh_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/hoh_pts_2855.gpkg"))
mas_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/mas_pts_2856.gpkg"))
col_dat <- data.frame(vect("SOIL CARBON/All_WA/data/points/col_pts_2855.gpkg"))

wa_dat <- rbind(hoh_dat, mas_dat, col_dat) |> 
  mutate(
  GEO = as.factor(GEO),
  geomorphons = as.factor(geomorphons),
  site = as.factor(site),
  site = fct_reorder(site, SOC_stock_spline, .fun = "median")) |>
  dplyr::rename_with(~gsub("_median", "", .x, fixed = TRUE))
wa_dat


## ------------------------------------------------------------------------------------------------------
unique(hoh_dat$GEO)
unique(mas_dat$GEO)
unique(col_dat$GEO)


## ----echo=FALSE----------------------------------------------------------------------------------------
wa_dat_num <- wa_dat |> 
    dplyr::select(SOC_stock_spline,
               #lower_depth,
               DTM,
               HLI,
               # NDYI_median, # Normalized difference yellow index
               # NDVI_median, # Normalized difference vegetation index
               #MNDWI, #Modified Normalized Difference Water Index
               EVI, # Enhanced vegetation index
               WIP) |> as.matrix()

ggcorrplot(cor(wa_dat_num), method = "square", type = "full", lab = T, lab_size = 3)



## ----echo=FALSE, fig.show='hold', out.width="35%"------------------------------------------------------

columns_to_exclude <- c("SOC_stock_spline", "lower_depth") 

wa_dat_scale <- wa_dat |> 
    dplyr::select(sample_ID, lower_depth, SOC_stock_spline, 
                  site, DTM, geomorphons, GEO, WIP,
                  EVI, MNDWI, HLI) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
                  ~dplyr::case_when(TRUE ~ scale(.))),
        site = fct_reorder(site, SOC_stock_spline, .fun = "median")) 
(wa_dat_scale)
hist(wa_dat_scale$SOC_stock_spline)
hist(log10(wa_dat_scale$SOC_stock_spline))



## ------------------------------------------------------------------------------------------------------
#Full model with all parameters
mod1 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP+EVI+DTM + (GEO) + geomorphons + site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
# No geomorphonssum
mod2 <-lmer(log10(SOC_stock_spline) ~ 
                WIP+EVI+DTM+HLI + (GEO) + site +
                 lower_depth + (lower_depth|sample_ID),  
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#No spectral, geomorphons
mod3 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP+DTM+HLI + (GEO) + site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#No geology or geomorphons
mod4 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP+DTM+HLI + site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#No WIP, geomorphons
mod5 <- lmer(log10(SOC_stock_spline) ~ 
                  DTM+HLI + (GEO) + site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#HLI interaction within site
mod6 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP + DTM +(GEO) + HLI*site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#DTM interaction with site
mod7 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP + HLI + site*DTM + (GEO) + 
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#Just geology and WIP
mod8 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP + (GEO) + 
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#Just WIP
mod9 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP + 
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#NULL
mod10 <- lmer(log10(SOC_stock_spline) ~ 
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))




## ----tidy=FALSE, echo=FALSE----------------------------------------------------------------------------
anova_table_mod1 <- rbind(anova(mod1, mod2), 
                            anova(mod1, mod3), 
                            anova(mod1, mod4), 
                            anova(mod1, mod5),
                            anova(mod1, mod6), 
                            anova(mod1, mod7), 
                            anova(mod1, mod8), 
                            anova(mod1, mod9), 
                            anova(mod1, mod10)) 

anova_table_mod1 |>
    tibble::rownames_to_column(var = "models") |> 
    dplyr::mutate(Significant = case_when(`Pr(>Chisq)` < 0.05 & `Pr(>Chisq)`>  0.01~ "***",
                                `Pr(>Chisq)` < 0.01 & `Pr(>Chisq)` > 0.001~ "***",
                                `Pr(>Chisq)` < 0.001 ~ "***",
                                .default = ""),
           Models = case_when(models == "mod11" ~ "mod1",
                              models == "mod12" ~ "mod1",
                              models == "mod13" ~ "mod1",
                              models == "mod14" ~ "mod1",
                              models == "mod15" ~ "mod1",
                              models == "mod16" ~ "mod1",
                              models == "mod17" ~ "mod1",
                              models == "mod18" ~ "mod1",
                              .default = models)) |> 
    dplyr::select(-models) |> 
    dplyr::select(Models, everything()) |> 
    #dplyr::arrange((AIC)) |> 
    kbl() |>
    #row_spec(, background = "lightblue") |> 
    kable_classic_2("hover", full_width = F) 


## ------------------------------------------------------------------------------------------------------
#DTM interaction with site
#No HLI
mod7.1 <- lmer(log10(SOC_stock_spline) ~ 
                 WIP + HLI + DTM*site + (GEO) +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#No GEO, HLI
mod7.2 <- lmer(log10(SOC_stock_spline) ~ 
                WIP + DTM +site/HLI + (GEO) +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))
#No WIP
mod7.3 <- lmer(log10(SOC_stock_spline) ~ 
                 DTM*site +
                 lower_depth + (lower_depth|sample_ID), 
            data = wa_dat_scale, REML = F,
            control=lmerControl(optimizer="bobyqa"))


## ----tidy=FALSE, echo=FALSE----------------------------------------------------------------------------
anova_table_mod2 <- rbind(anova(mod7, mod7.1), #Not Significant
anova(mod7, mod7.2), 
anova(mod7, mod7.3), 
make.row.names = T) #Not Significant

anova_table_mod2 |>
    tibble::rownames_to_column(var = "models") |> 
    mutate(Significant = case_when(`Pr(>Chisq)` < 0.05 & `Pr(>Chisq)`>  0.01~ "***",
                                `Pr(>Chisq)` < 0.01 & `Pr(>Chisq)` > 0.001~ "***",
                                `Pr(>Chisq)` < 0.001 ~ "***",
                                .default = ""),
           Models = case_when(models == "mod71" ~ "mod7",
                              models == "mod72" ~ "mod7",
                              models == "mod73" ~ "mod7",
                              .default = models)) |>
    dplyr::select(-models) |> 
    dplyr::select(Models, everything()) |> 
    kbl() |>
    kable_classic_2("hover", full_width = F)



## ----include=FALSE-------------------------------------------------------------------------------------
round(r.sq(log10(wa_dat_scale$SOC_stock_spline), fitted(mod7)), 3)


## ----cache=TRUE, echo=FALSE, out.width="100%"----------------------------------------------------------
tab_model(mod7, bootstrap = T, seed = 11, dv.labels = c("Model 7"))


## ----anova---------------------------------------------------------------------------------------------
lmerTest::ranova(mod7)
lmerTest::difflsmeans(mod7, test.effs = "site", which = "site")
emmeans::emmeans(mod7, list(pairwise ~ site), adjust = "tukey")


## ----predplot, echo=FALSE, warning=FALSE, message=FALSE------------------------------------------------
library(RColorBrewer)
(mod7graph <- ggplot(wa_dat, 
                      aes(y = 10**(fitted(mod7)), x = (SOC_stock_spline))) +
    geom_jitter(color='black', 
                aes(fill = (wa_dat$WIP*100), 
                    shape = as.factor(site)),
               size = 4, stroke = 0.9, alpha = 0.8) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"),
                         name = "Wetland \nProbability %", n.breaks = 5, limits = c(0, 100)) +
     scale_shape_manual(name = "Study Area", 
                        values = c(21, 22, 23, 24)) +
    geom_smooth(aes(y = 10**(fitted(mod7)), x = (SOC_stock_spline)), 
                method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                linewidth = 0.9, linetype = 5, alpha = 0.3, se = T) +
    xlab(expression('Sampled SOC Stock (g cm'^-2*')')) + 
    ylab(expression('Predicted SOC Stock (g cm'^-2*')')) + 
    geom_abline(intercept = 0, slope = 1, linewidth = 0.9, linetype = "dashed") +
    annotate("text", 
             label = paste("R^{2} == ", 
                           round(r.sq((wa_dat$SOC_stock_spline),
                                      10**fitted(mod7)), 3)), 
             x = 0.5, y = 3.5, size = 4, parse = T) + 
    annotate("text", label = "Model Fit", 
             x = 1.5, y = 3.5, size = 4) +
    annotate("segment", color = "#fa3e3e",
             x = 1.3, y = 3.3, linewidth = 0.9, linetype = 5, 
             xend = 1.75) + 
    annotate("text", label = "1:1", 
             x = 2.5, y = 3.5, size = 4) +
    annotate("segment", color = "black",
             x = 2.3, y = 3.3, linewidth = 0.9, linetype = "dashed",
             xend = 2.75) +
    xlim(0, 4) +
    ylim(0, 4)  +
    theme(legend.position = 'right', 
          legend.key.size = unit(0.6, "cm"),
          legend.spacing.x = unit(1.2, "cm"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank(),
          text = element_text(size = 16)) +
    guides(guide_legend(byrow = TRUE))
)



## ----echo=FALSE, warning=FALSE, message=FALSE, include=FALSE-------------------------------------------
ggplot2::ggsave(plot = mod7graph, filename = paste0("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/SOIL CARBON/All_WA/writing/Figures/model7_nonlog_SOCprediction.png"),
                width = 9, height = 7.5, units = "in", dpi = 500)


## ----fig.show='hold', fig.align='center', echo=FALSE, message=FALSE------------------------------------
library(ggpubr)

partial_plot <- function(data, model, x_var, color_var, shape_var){
  if(is.factor(data[x_var][[1]]) == TRUE) {
  x_var <- sym(x_var)
  shape_var <- sym(shape_var)
  color_var <- sym(color_var) 
  fitted <- predict(model, newdata = data)
  
  ggplot(data, aes(y = 10**fitted, x = as.factor(!!x_var))) +
    geom_violin(show.legend = FALSE, scale = "width", linewidth = 0.9)+
    geom_jitter(aes(shape = !!shape_var, , color = !!color_var), width = 0.2, size = 2)+
    xlab("Study Area\n") + 
    ylab(expression('Fitted Model SOC Stock (g cm'^-2*')')) +
    theme(legend.position = 'right', 
          legend.key.size = unit(0.6, "cm"),
          legend.spacing.x = unit(1.2, "cm"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank(),
          text = element_text(size = 16))
  } else {
    if(x_var == "WIP"){
      label <- "Wetland Intrinsic Potential \n Scaled"
    } else if(x_var == "DTM") {
      label <- "Elevation Scaled\n"
    } else {
      label <- paste0("\n", x_var)
    }
  x_var <- sym(x_var)
  shape_var <- sym(shape_var)
  color_var <- sym(color_var) 
  fitted <- predict(model, newdata = data)
  
  ggplot(data, aes(y = 10**fitted, x = !!x_var)) +
    geom_point(aes(shape = !!shape_var, color = !!color_var), size = 2)+
    xlab(label) +
    ylab(expression('Fitted Model SOC Stock (g cm'^-2*')')) +
    geom_smooth(aes(y = 10**fitted, x = !!x_var), 
                method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                linewidth = 0.9, linetype = 5, alpha = 0.3, se = T) + 
    theme(legend.position = 'top', 
          legend.key.size = unit(0.8, "cm"),
          legend.spacing.x = unit(1.2, "cm"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank(),
          text = element_text(size = 16))
  }
}

mod7site <- partial_plot(data = wa_dat_scale, model = mod7, x_var = "site", shape_var = "site", color_var = "HLI")

mod7wip <- partial_plot(data = wa_dat_scale, model = mod7, x_var = "WIP", shape_var = "site", color_var = "DTM")

mod7dtm <- partial_plot(data = wa_dat_scale, model = mod7, x_var = "DTM", shape_var = "site", color_var = "HLI")

mod7hli <- partial_plot(data = wa_dat_scale, model = mod7, x_var = "HLI", shape_var = "site", color_var = "DTM")

ggarrange(mod7site, mod7wip+rremove("ylab"), mod7hli+rremove("ylab"), ncol = 3, nrow=1, common.legend = T)


## ----echo=FALSE, message=FALSE-------------------------------------------------------------------------
wa_dat_scale_wet <- wa_dat |> 
    filter((WIP >= 0.5 )) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
                  ~dplyr::case_when(TRUE ~ scale(.))))
    
wa_dat_scale_upl <- wa_dat |> 
    filter((WIP < 0.5 )) |> 
    dplyr::mutate(across(
        dplyr::where(is.numeric) & !all_of(columns_to_exclude),
                  ~dplyr::case_when(TRUE ~ scale(.))))
    


## ----echo=FALSE, message=FALSE-------------------------------------------------------------------------
mod7_wet <- update(mod7, data = wa_dat_scale_wet)

mod7wet_pred <- predict(mod7, newdata = wa_dat_scale_wet)
#summary(mod2.1_wet)
mod7_upl <- update(mod7, data = wa_dat_scale_upl)

mod7upl_pred <- predict(mod7, newdata = wa_dat_scale_upl)
#summary(mod2.1_upl)
tab_model(mod7_wet, mod7_upl, dv.labels = c("Wetlands (WIP > 0.50", "Uplands (WIP < 0.50)"))


## ----echo=FALSE, message=FALSE-------------------------------------------------------------------------
library(ppcor)

rbind(pcor.test(wa_dat_scale_wet$WIP, wa_dat_scale_wet$SOC_stock_spline, wa_dat_scale_wet$lower_depth),
pcor.test(wa_dat_scale_upl$WIP, wa_dat_scale_upl$SOC_stock_spline, wa_dat_scale_upl$lower_depth))

rbind(pcor.test(wa_dat_scale_wet$DTM, wa_dat_scale_wet$SOC_stock_spline, wa_dat_scale_wet$lower_depth),
pcor.test(wa_dat_scale_upl$DTM, wa_dat_scale_upl$SOC_stock_spline, wa_dat_scale_upl$lower_depth))

rbind(pcor.test(wa_dat_scale_wet$HLI, wa_dat_scale_wet$SOC_stock_spline, wa_dat_scale_wet$lower_depth),
pcor.test(wa_dat_scale_upl$HLI, wa_dat_scale_upl$SOC_stock_spline, wa_dat_scale_upl$lower_depth))



## ----echo=FALSE, message=FALSE-------------------------------------------------------------------------
library(ggpubr)

DTM_wet <- partial_plot(data = wa_dat_scale_wet, 
             model = mod7, 
             x_var = "DTM", 
             color_var = "site", 
             shape_var = "GEO")

wip_wet <- partial_plot(data = wa_dat_scale_wet, 
             model = mod7, 
             x_var = "WIP", 
             color_var = "site", 
             shape_var = "GEO")

DTM_upl <- partial_plot(data = wa_dat_scale_upl, 
             model = mod7, 
             x_var = "DTM", 
             color_var = "site", 
             shape_var = "GEO")

wip_upl <- partial_plot(data = wa_dat_scale_upl, 
             model = mod7, 
             x_var = "WIP", 
             color_var = "site", 
             shape_var = "GEO")

wetupl_graphs <- ggarrange(DTM_wet+rremove("xlab"), wip_wet+rremove("xlab")+rremove("ylab"), 
                           DTM_upl, wip_upl+rremove("ylab"), ncol = 2, nrow=2, common.legend = T)
annotate_figure(
    annotate_figure(wetupl_graphs, 
                    left = grid::textGrob("Upland", 
                                      rot = 90, vjust = 0.5, hjust = 3, 
                                      gp = grid::gpar(cex = 1.3))), 
    left = grid::textGrob("Wetland", 
                    rot = 90, vjust = 2.1, hjust = -1.5, gp = grid::gpar(cex = 1.3)))

ggplot2::ggsave(plot = wetupl_graphs, filename = paste0("/Users/Anthony/OneDrive - UW/University of Washington/Presentations/Conferences/AGU2023/ExploringControls/Figures and Tables/mod2.1_wet_upl.jpeg"),
                width = 9, height = 7.5, units = "in", dpi = 500)



## ----fig.show='hold', message=FALSE, echo=FALSE--------------------------------------------------------
library(interactions)
DTMsite <- interactions::interact_plot(mod7, pred = DTM, modx = site, plot.points = T, plot.shape = T,
                            x.label = "DTM", 
                            y.label = "Fitted Model SOC Stock",
                            legend.main = "Site",
                            interval = TRUE, data = wa_dat_scale)
DTMsitewet <- interactions::interact_plot(mod7_wet, pred = DTM, modx = site, plot.points = T, plot.shape = T,
                            x.label = "DTM", 
                            y.label = "Fitted Model SOC Stock",
                            legend.main = "Site",
                            interval = TRUE, data = wa_dat_scale_wet)
DTMsiteupl <- interactions::interact_plot(mod7_upl, pred = DTM, modx = site, plot.points = T, plot.shape = T,
                            x.label = "DTM", 
                            y.label = "Fitted Model SOC Stock",
                            legend.main = "Site",
                            interval = TRUE, data = wa_dat_scale_upl)

DTMsite
ggarrange(DTMsite, DTMsiteupl, DTMsitewet)


## ------------------------------------------------------------------------------------------------------
hoh_dat_scale <- wa_dat_scale |> filter(site == "HOH")
mas_dat_scale <- wa_dat_scale |> filter(site == "MAS")
col_dat_scale <- wa_dat_scale |> filter(site == "COL")

hoh_datpred <- predict(mod7, newdata = hoh_dat_scale)
mas_datpred <- predict(mod7, newdata = mas_dat_scale)
col_datpred <- predict(mod7, newdata = col_dat_scale)


## ------------------------------------------------------------------------------------------------------
partial_plot(data = hoh_dat_scale, 
             model = mod7, 
             x_var = "HLI", 
             color_var = "lower_depth", 
             shape_var = "GEO")
partial_plot(data = mas_dat_scale, 
             model = mod7, 
             x_var = "HLI", 
             color_var = "lower_depth", 
             shape_var = "GEO")

partial_plot(data = col_dat_scale, 
             model = mod7, 
             x_var = "HLI", 
             color_var = "lower_depth", 
             shape_var = "GEO")

