library(terra)
library(sf)
library(dplyr)
library(tidyterra)
library(corrplot)
library(lme4)
library(lmerTest)
library(emmeans)
library(r2glmm)
library(ggcorrplot)
library(merTools)
library(MuMIn)
library(mgcv)
library(mpspline2)

setwd("/Users/Anthony/OneDrive - UW/University of Washington/Data and Modeling/")
# Quick R^2 function
r.sq <- function(y,y.fitted){
    res <- y-y.fitted
    1-sum(res^2)/sum((y-mean(y))^2)
}






wa_hor_dat <- read.csv("SOIL CARBON/ANALYSIS/All_WA_horizons_lithsum2.csv")
str(wa_hor_dat)

wa_dat <- wa_hor_dat |> dplyr::select("sample_ID", "lat","lon","depth_cm", "thickness_cm",
                                      "rock_perc","BD_g_cm3","field_texture","field_texture_binned",
                                      "redox","carbon_perc","carbon_stock_g_cm2",
                                      "ANDWI_median", 
                                      "DSI_median",
                                      "DSWI1_median",           
                                      "EMBI_median",
                                      "EVI_median",
                                      "LITHOL","LSWI_median",          
                                      "MNDWI_median",         
                                      "NDVI_median",          
                                      "NDYI_median",          
                                      "SAVI_median",          
                                      "SCI_median",           
                                      "WIP",                 
                                      "site", "lithsum",
                                      "temp", "precip") |>
    #mutate_at("geomorph", as.factor) |>
    mutate_if(is.character, as.factor) |>
    #mutate_if(is.numeric, scale) |>
    stats::setNames(c("sample_ID", "lat","lon","depth", "thickness",
                      "rock_perc","BD","field_texture","field_texture_binned",
                      "redox","SOC_perc","SOC_stock",
                      "ANDWI_median", 
                      "DSI_median",
                      "DSWI1_median",           
                      "EMBI_median",
                      "EVI_median",
                      "LITHOL","LSWI_median",          
                      "MNDWI_median",         
                      "NDVI_median",          
                      "NDYI_median",          
                      "SAVI_median",          
                      "SCI_median",           
                      "WIP",                 
                      "site", "lithsum",
                      "temp", "precip")) |> 
    mutate_at(c("temp", "precip", "depth"), scale) |>
    # mutate(TCC_fact = case_when(tree_canopy_cover>= 70 ~ "HFOR", 
    #                             tree_canopy_cover >= 50 & tree_canopy_cover < 70 ~ "MFOR",
    #                             tree_canopy_cover >= 25 & tree_canopy_cover < 50 ~ "LFOR", 
    #                             tree_canopy_cover < 25 ~ "NON",)) |> 
    as.data.frame()

str(wa_dat)



###### Parameter selection #######
    # Look for high corr values with SOC stock
wa_dat_num <- wa_dat |> select_if(is.numeric) |> as.matrix()

ggcorrplot(cor(wa_dat_num), method = "square", type = "full", lab = T, lab_size = 2)



###### Model parameters #######
hist(wa_dat$SOC_perc)
hist(wa_dat$SOC_stock)

gmod <- lmer(log(SOC_stock) ~  WIP +(lithsum)*temp*precip + (1+(depth)|sample_ID), data = wa_dat, na.action = 'na.fail')
x <- dredge(gmod, fixed = c("WIP", "lithsum", "temp", "precip"), beta = "sd" )
print(x, abbrev.names = F)
dmod <- get.models(x, 1)[[1]]


qqnorm(y = resid(dmod))
qqline(y = resid(dmod))

summary(dmod)
summary(gmod)
r.sq(log(wa_dat$SOC_stock), (fitted(dmod)))
plot(log(wa_dat$SOC_stock), (fitted(dmod)))
abline(0,1)

emmeans(dmod, "lithsum", type = "response", na.rm = T)
confint(dmod, oldNames = F)
merTools::RMSE.merMod(dmod, scale = T)*sd((wa_dat$SOC_stock))
r.squaredGLMM(dmod)
r2glmm::r2beta(dmod, partial = T, method = 'nsj')

gmod2 <- lmer(log(SOC_stock) ~  WIP  +(lithsum) +  (site) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod3 <- lmer(log(SOC_stock) ~  NDYI_median  +(lithsum) +  (site) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod4 <- lmer(log(SOC_stock) ~  NDYI_median + WIP   +  (site) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod5 <- lmer(log(SOC_stock) ~  NDYI_median + WIP  +(lithsum) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod6 <- lmer(log(SOC_stock) ~  (lithsum) +  (site) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod7 <- lmer(log(SOC_stock) ~  NDYI_median + WIP   + (1+scale(depth)|sample_ID), data = wa_dat)
gmod8 <- lmer(log(SOC_stock) ~  (site) + (1+scale(depth)|sample_ID), data = wa_dat)
gmod9 <- lmer(log(SOC_stock) ~  NDYI_median:WIP  +(lithsum) +  (site) + (1+scale(depth)|sample_ID), data = wa_dat)

anova(gmod, gmod2)
anova(gmod, gmod3)
anova(gmod, gmod4)
anova(gmod, gmod5)
anova(gmod, gmod6)
anova(gmod, gmod7)
anova(gmod, gmod8)
anova(gmod, gmod9)


############################################################################
library(RColorBrewer)
ggplot(wa_dat, aes(y = (fitted(dmod)), x = log(SOC_stock))) +
    geom_point(color='black', shape=21, aes(fill = (WIP*100)),
               # #shape = as.factor(WIP >= 0.5), 
               size = 3, stroke = 0.5, alpha = 0.8) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), 
                         name = "WIP %", n.breaks = 5, limits = c(0, 100)) +
    geom_smooth(aes(y = (fitted(dmod)), x = log(SOC_stock)), 
                method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                alpha = 0.2, size = 1.2, linetype = 5, alpha = 0.3, se = T) +
    xlab(expression('Sampled SOC % (MgC ha'^-1*')')) +  
    ylab(expression('Predicted 1m SOC % (MgC ha'^-1*')')) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed") +
    #xlim(-4, 5) +
    #ylim(-4, 5)  +
    guides(guide_legend(byrow = TRUE)) +
    theme(legend.position = 'right', 
          legend.key.size = unit(1.7, "cm"),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank(),
          text = element_text(size = 18))


ggplot(wa_dat, aes(y = (fitted(dmod)), x = log(SOC_stock))) +
    geom_jitter(color='black', shape=21, aes(fill = site),
               # #shape = as.factor(WIP >= 0.5), 
               size = 3, stroke = 0.5, alpha = 0.8) +
    # scale_fill_gradientn(colours = brewer.pal(9, "YlGnBu"), 
    #                      name = "WIP %", n.breaks = 5, limits = c(0, 100)) +
    geom_smooth(aes(y = (fitted(dmod)), x = log(SOC_stock)), 
                method = "lm", color = "#fa3e3e", fill = "#fa3e3e", 
                alpha = 0.2, size = 1.2, linetype = 5, alpha = 0.3, se = T) +
    xlab(expression('Sampled SOC % (MgC ha'^-1*')')) +  
    ylab(expression('Predicted 1m SOC % (MgC ha'^-1*')')) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, linetype = "dashed") +
    #xlim(-4, 5) +
    #ylim(-4, 5)  +
    guides(guide_legend(byrow = TRUE)) +
    theme(legend.position = 'right', 
          legend.key.size = unit(1, "cm"),
          legend.key = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey80"),
          axis.ticks = element_blank(),
          text = element_text(size = 10))

wa_dat_sum <- wa_dat |> dplyr::summarise(stock_sum = sum(SOC_stock)*100, 
                                         #stock_check = mean(SOC120),
                                         WIP = mean(WIP),
                                         n = n(), 
                                         .by = sample_ID)


########## Wetlands only defined by WIP #########

wa_wet <- wa_dat |> filter(WIP >= 0.50)

dmod_wet <- update(dmod, data = wa_wet)
summary(dmod_wet)
emmeans(dmod_wet, "lithsum", type = "response")
confint(dmod_wet, oldNames = F)
merTools::RMSE.merMod(dmod_wet, scale = T)*sd((wa_dat$SOC_stock))
r.squaredGLMM(dmod_wet)
r2glmm::r2beta(dmod_wet, partial = T, method = 'nsj')

########## Uplands only defined by WIP #########

wa_upl <- wa_dat |> filter(WIP < 0.50)

dmod_upl <- update(dmod, data = wa_upl)
summary(dmod_upl)
emmeans(dmod_upl, "lithsum", type = "response")
confint(dmod_upl, oldNames = F)
merTools::RMSE.merMod(dmod_upl, scale = T)*sd((wa_dat$SOC_stock))
r.squaredGLMM(dmod_upl)
r2glmm::r2beta(dmod_upl, partial = T, method = 'nsj')


########### GAM ###############

gam_mod <- gamm(log(SOC_stock) ~ s(WIP) + site + lithsum + s(depth, sample_ID, bs = "re"), data = wa_dat)
summary(gam_mod$gam)
mgcv::qq.gam(gam_mod$gam)

plot(log(wa_dat$SOC_stock), fitted(gam_mod$gam))
abline(0,1)
r.sq(log(wa_dat$SOC_stock), fitted(gam_mod$gam))
gam.vcomp(gam_mod$gam, conf.lev = .95)
plot.gam(gam_mod$gam)


########## Random Forest ############
library(randomForest)
set.seed(11)

# Validation Set 
train.index <- as.vector(sample(c(1:481), 337, replace=F))
train <- wa_dat[train.index, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]
test <- wa_dat[-train.index, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]

train.indexw <- as.vector(sample(c(1:180), 126, replace=F))
trainw <- wa_wet[train.indexw, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]
testw <- wa_wet[-train.indexw, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]

train.indexu <- as.vector(sample(c(1:301), 211, replace=F))
trainu <- wa_dat[train.indexu, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]
testu <- wa_dat[-train.indexu, c("SOC_stock", "WIP", "lithsum", "temp", "precip", "depth")]


rf_model <- randomForest(log(SOC_stock) ~ .,
                         ntree = 1000, importance = TRUE, data = train)
rf_modelw <- randomForest(log(SOC_stock) ~ .,
                         ntree = 1000, importance = TRUE, data = trainw)
rf_modelu <- randomForest(log(SOC_stock) ~ .,
                         ntree = 1000, importance = TRUE, data = trainu)


rf_model
rf_modelw
rf_modelu

plot(rf_model)
plot(rf_modelw)
plot(rf_modelu)

varImpPlot(rf_model)
varImpPlot(rf_modelw)
varImpPlot(rf_modelu)

rf.pred <- predict(rf_model, newdata = test)
mean((rf.pred - log(test$SOC_stock))^2)
rf.full <- predict(rf_model, newdata = wa_dat)

plot(log(wa_dat$SOC_stock), rf.full)
r.sq(log(wa_dat$SOC_stock), rf.full)

############# Boosted Regression Trees ##############
library(gbm)

boost <- gbm(log(SOC_stock) ~ ., data = train,
                                distribution = "gaussian", 
                                n.trees = 500, interaction.depth=4,
                                shrinkage = 0.1) # learning rate
boostw <- gbm(log(SOC_stock) ~ ., data = trainw,
             distribution = "gaussian", 
             n.trees = 500, interaction.depth=4,
             shrinkage = 0.1) # learning rate
boostu <- gbm(log(SOC_stock) ~ ., data = trainu,
             distribution = "gaussian", 
             n.trees = 500, interaction.depth=4,
             shrinkage = 0.1) # learning rate

summary(boost)
summary(boostw)
summary(boostu)

boost.pred <- predict(boost, newdata=test)
mean((boost.pred - log(test$SOC_stock))^2)
boost.full <- predict(boost, newdata = wa_dat)

plot(log(wa_dat$SOC_stock), boost.full)
r.sq(log(wa_dat$SOC_stock), boost.full)
abline(0,1)


############# Digital soil Mapping model #############
wa_spl_dat <- wa_hor_dat |> dplyr::select("sample_ID", "depth_cm", "lat","lon","thickness_cm",
                                      "rock_perc","BD_g_cm3","field_texture","field_texture_binned",
                                      "redox","carbon_perc","carbon_stock_g_cm2",
                                      "ANDWI_median", 
                                      "DSI_median",
                                      "DSWI1_median",           
                                      "EMBI_median",
                                      "EVI_median",
                                      "LITHOL","LSWI_median",          
                                      "MNDWI_median",         
                                      "NDVI_median",          
                                      "NDYI_median",          
                                      "SAVI_median",          
                                      "SCI_median",           
                                      "WIP",                 
                                      "site", "lithsum",
                                      "temp", "precip") |>
    #mutate_at("geomorph", as.factor) |>
    #mutate_if(is.character, as.factor) |>
    #mutate_if(is.numeric, scale) |>
    stats::setNames(c("sample_ID", "depth","lat","lon", "thickness",
                      "rock_perc","BD","field_texture","field_texture_binned",
                      "redox","SOC_perc","SOC_stock",
                      "ANDWI_median", 
                      "DSI_median",
                      "DSWI1_median",           
                      "EMBI_median",
                      "EVI_median",
                      "LITHOL","LSWI_median",          
                      "MNDWI_median",         
                      "NDVI_median",          
                      "NDYI_median",          
                      "SAVI_median",          
                      "SCI_median",           
                      "WIP",                 
                      "site", "lithsum",
                      "temp", "precip")) |> 
    #mutate_at(c("temp", "precip", "depth"), scale) |>
    dplyr::group_by(sample_ID) |> 
    mutate(top = case_when(is.na(depth - lag(depth)) ~ 0,
                           .default = lag(depth)),
           bottom = depth,
           center = abs(top - (top - bottom)/2)) |> 
    dplyr::select(sample_ID, top, bottom, center, everything()) |> 
    as.data.frame()
str(wa_spl_dat)

wa_spl <- mpspline_tidy(obj = wa_spl_dat, var_name = "SOC_stock")
t <- head(wa_spl$est_dcm, 6)
str(wa_spl)

test_spl <- mpspline_tidy(obj = wa_spl_dat, var_name = "SOC_stock", d = c(0, 100, 200))
head(test_spl$est_dcm)
