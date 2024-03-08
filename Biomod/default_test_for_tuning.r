#---------------------------------------------------------
# Load packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(terra)
library(mgcv)
#install.packages("rasterVis")
library(rasterVis)

#---------------------------------------------------------
# Set working directory and import all occ data for species
setwd("C:\\Users\\chris\\Desktop\\Meg full git repo")
data <- read.csv("input_data\\occ_data\\20240219_Sampled by raster.csv")
head(data)
print(nrow(data))

#---------------------------------------------------------
# Load current climate data (global)
ch_bio1_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio1_clipped.tif")
ch_bio4_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio4_clipped.tif")
ch_bio10_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio10_clipped.tif")
ch_bio11_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio11_clipped.tif")
ch_bio15_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio15_clipped.tif")
ch_bio16_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio16_clipped.tif")
ch_bio17_raster <- raster("input_data\\climate_data\\CURRENT_DATA\\bio17_clipped.tif")

plot(ch_bio1_raster)
#---------------------------------------------------------
# Create stack of climatic data
ch_current_clim <- stack(
                        ch_bio1_raster,
                        ch_bio4_raster,
                        ch_bio10_raster,
                        ch_bio11_raster,
                        ch_bio15_raster,
                        ch_bio16_raster,
                        ch_bio17_raster
                        )

names(ch_current_clim) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(ch_current_clim)

#---------------------------------------------------------
#Format input data, and select pseudo-absences if wanted, for usage in biomod2
ch_SpheCoop_data <- BIOMOD_FormatingData(
                        resp.var = rep(1, nrow(data)),
                        expl.var = ch_current_clim,
                        resp.xy = data[,c('decimalLongitude','decimalLatitude')],
                        resp.name = 'Sphaeropteris.cooperi',
                        PA.nb.rep = 3,
                        PA.nb.absences = 1281,
                        PA.strategy = 'random'
                        )
ch_SpheCoop_data
plot(ch_SpheCoop_data)
#---------------------------------------------------------
# Building the individual models
ch_SpheCoop_models <- BIOMOD_Modeling(
                        bm.format = ch_SpheCoop_data,
                        models = c("GLM","GBM","GAM","ANN","CTA","FDA","RF",'MARS', 'MAXENT'),
                        nb.rep = 4,
                        CV.nb.rep = 4,
                        data.split.perc = 80,
                        metric.eval = c('TSS', 'ROC', 'KAPPA'),
                        var.import = 3,
                        do.full.models = F,
                        scale.models = TRUE,
                        modeling.id = "SpheCoop-models",
                        seed.val = 42
                        )

# Getting model scores for the individual models trained
indiv_SpheCoop_scores <- get_evaluations(ch_SpheCoop_models)
indiv_SpheCoop_scores
write.csv(indiv_SpheCoop_scores, file = "DEFAULT-Model-evaluations.csv", row.names = FALSE)

#---------------------------------------------------------
# Varaiable importances
importances<-(SpheCoop_models_var_import <- get_variables_importance(ch_SpheCoop_models))
write.csv(importances, file = "DEFAULT-Model-importances.csv", row.names = FALSE)
head(importances)

#---------------------------------------------------------
# Visualizing the predictive accuracy of the models
# By PA
bm_PlotEvalMean(
    ch_SpheCoop_models, 
    metric.eval = c('TSS','ROC'),
    group.by ='PA',
    do.plot = TRUE,
    xlim = c(0,2),
    ylim = c(0.9,1.05)
    )

# By Run
bm_PlotEvalMean(
    ch_SpheCoop_models, 
    metric.eval = c('TSS','ROC'),
    group.by ='run',
    do.plot = TRUE,
    xlim = c(0,2),
    ylim = c(0.9,1.05)
    )

# By Algorithm
bm_PlotEvalMean(
    ch_SpheCoop_models, 
    metric.eval = c('TSS','ROC'),
    group.by ='algo',
    do.plot = TRUE,
    xlim = c(0,2),
    ylim = c(0.75,1.05)
    )

#---------------------------------------------------------
# Calculate mean variable importance by algorithm
library(dplyr)
library(tidyr)
mean_importance <- importances %>%
  group_by(algo, expl.var) %>%
  summarize(mean_var_imp = mean(var.imp, na.rm = TRUE)) %>%
  spread(key = algo, value = mean_var_imp)
print(mean_importance)

#---------------------------------------------------------
# Load models and plot response curve for all models
SpheCoop_glm_subset <- get_built_models(ch_SpheCoop_models, algo = 'GLM')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_glm_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )

SpheCoop_ann_subset <- get_built_models(ch_SpheCoop_models, algo = 'ANN')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_ann_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )

SpheCoop_gbm_subset <- get_built_models(ch_SpheCoop_models, algo = 'GBM')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_gbm_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )

SpheCoop_gam_subset <- get_built_models(ch_SpheCoop_models, algo = 'GAM')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_gam_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )
SpheCoop_cta_subset <- get_built_models(ch_SpheCoop_models, algo = 'CTA')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_cta_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )
SpheCoop_fda_subset <- get_built_models(ch_SpheCoop_models, algo = 'FDA')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_fda_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )
SpheCoop_mars_subset <- get_built_models(ch_SpheCoop_models, algo = 'MARS')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_mars_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )
SpheCoop_maxent_subset <- get_built_models(ch_SpheCoop_models, algo = 'MAXENT')

response_curves <- bm_PlotResponseCurves(
    bm.out = ch_SpheCoop_models,
    models.chosen = SpheCoop_maxent_subset,
    new.env = get_formal_data(ch_SpheCoop_models,'expl.var'),
    show.variables = get_formal_data(ch_SpheCoop_models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var = 'mean',
    legend = FALSE,
    display_title = FALSE,
    do.plot = TRUE
    )

#---------------------------------------------------------
# Train ensemble model
ch_SpheCoop_ensemble <- BIOMOD_EnsembleModeling(
                        bm.mod = ch_SpheCoop_models,
                        em.by = 'all',
                        em.algo = c('EMmean','EMmedian','EMca','EMwmean'),
                        metric.select.thresh = c(0.8),
                        metric.select = c('TSS'),
                        var.import = 0,
                        seed.val = 42
                        )

ensemble_score_csv <- (ensemble_scores <- get_evaluations(ch_SpheCoop_ensemble))
print(ensemble_score_csv)
# Write the dataframe to a CSV file
write.csv(ensemble_score_csv, file = "DEFAULT-Model-ensemble-evaluations.csv", row.names = FALSE)
