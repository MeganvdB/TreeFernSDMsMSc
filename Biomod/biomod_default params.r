#---------------------------------------------------------
# Load packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(terra)
library(mgcv)
install.packages("rasterVis")
library(rasterVis)

#---------------------------------------------------------
# Set working directory and import all occ data for species
setwd("C:\\Users\\chris\\Desktop\\Meg full git repo")
data <- read.csv("input_data\\occ_data\\glob_occ_data_with_bio.csv")
head(data)
print(nrow(data))

#---------------------------------------------------------
# Load current climate data (global)
ch_bio1_raster <- raster("input_data\\climate_data\\current\\bio1_clipped.tif")
ch_bio4_raster <- raster("input_data\\climate_data\\current\\bio4_clipped.tif")
ch_bio10_raster <- raster("input_data\\climate_data\\current\\bio10_clipped.tif")
ch_bio11_raster <- raster("input_data\\climate_data\\current\\bio11_clipped.tif")
ch_bio15_raster <- raster("input_data\\climate_data\\current\\bio15_clipped.tif")
ch_bio16_raster <- raster("input_data\\climate_data\\current\\bio16_clipped.tif")
ch_bio17_raster <- raster("input_data\\climate_data\\current\\bio17_clipped.tif")

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
                        PA.nb.absences = 2964,
                        PA.strategy = 'random'
                        )
ch_SpheCoop_data
plot(ch_SpheCoop_data)

#---------------------------------------------------------
# Default hyperparameters
ch_SpheCoop_model_options <- BIOMOD_ModelingOptions(
                        MAXENT = list(
                            path_to_maxent.jar = 'C:\\Users\\chris\\Desktop\\MegMSc\\maxent',
                            visual = FALSE
                        )
                    )
                        
#---------------------------------------------------------
# Building the individual models
ch_SpheCoop_models <- BIOMOD_Modeling(
                        bm.format = ch_SpheCoop_data,
                        models = c("GLM","GBM","GAM","ANN","CTA","FDA",'MARS', 'MAXENT'),
                        bm.options = ch_SpheCoop_model_options,
                        nb.rep = 4,
                        CV.nb.rep = 4,
                        data.split.perc = 80,
                        metric.eval = c('TSS', 'ROC', 'KAPPA'),
                        var.import = 3,
                        do.full.models = F,
                        scale.models = TRUE,
                        modeling.id = "GFDL-ESM4-model_run",
                        seed.val = 42
                        )

# Getting model scores for the individual models trained
indiv_SpheCoop_scores <- get_evaluations(ch_SpheCoop_models)
indiv_SpheCoop_scores
write.csv(indiv_SpheCoop_scores, file = "default-Model-evaluations.csv", row.names = FALSE)

#---------------------------------------------------------
# Varaiable importances
importances<-(SpheCoop_models_var_import <- get_variables_importance(ch_SpheCoop_models))
write.csv(importances, file = "default-Model-importances.csv", row.names = FALSE)
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
    ylim = c(0.9,1.05)
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
write.csv(ensemble_score_csv, file = "default-Model-ensemble-evaluations.csv", row.names = FALSE)

#---------------------------------------------------------
# Loading raster data
# Load Current GR data
bio1_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio1_GR.tif")
bio4_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio4_GR.tif")
bio10_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio10_GR.tif")
bio11_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio11_GR.tif")
bio15_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio15_GR.tif")
bio16_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio16_GR.tif")
bio17_GR_current <- raster("input_data\\climate_data\\current\\garden_route\\bio17_GR.tif")

# Load future data ssp126
bio1_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio1_GR.tif")
bio4_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio4_GR.tif")
bio10_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio10_GR.tif")
bio11_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio11_GR.tif")
bio15_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio15_GR.tif")
bio16_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio16_GR.tif")
bio17_GR_2040_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp126\\garden_route\\bio17_GR.tif")

bio1_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio1_GR.tif")
bio4_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio4_GR.tif")
bio10_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio10_GR.tif")
bio11_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio11_GR.tif")
bio15_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio15_GR.tif")
bio16_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio16_GR.tif")
bio17_GR_2100_126 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp126\\garden_route\\bio17_GR.tif")

# Load future data ssp585
bio1_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio1_GR.tif")
bio4_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio4_GR.tif")
bio10_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio10_GR.tif")
bio11_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio11_GR.tif")
bio15_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio15_GR.tif")
bio16_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio16_GR.tif")
bio17_GR_2040_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2040-2070\\ssp585\\garden_route\\bio17_GR.tif")

bio1_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio1_GR.tif")
bio4_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio4_GR.tif")
bio10_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio10_GR.tif")
bio11_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio11_GR.tif")
bio15_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio15_GR.tif")
bio16_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio16_GR.tif")
bio17_GR_2100_585 <- raster("input_data\\climate_data\\GFDL-ESM4\\2071-2100\\ssp585\\garden_route\\bio17_GR.tif")

# Create stacks of climatic data
GR_clim_current <- stack(
                        bio1_GR_current,
                        bio4_GR_current,
                        bio10_GR_current,
                        bio11_GR_current,
                        bio15_GR_current,
                        bio16_GR_current,
                        bio17_GR_current
                        )
GR_clim_2040_126 <- stack(
                        bio1_GR_2040_126,
                        bio4_GR_2040_126,
                        bio10_GR_2040_126,
                        bio11_GR_2040_126,
                        bio15_GR_2040_126,
                        bio16_GR_2040_126,
                        bio17_GR_2040_126
                        )
GR_clim_2100_126 <- stack(
                        bio1_GR_2100_126,
                        bio4_GR_2100_126,
                        bio10_GR_2100_126,
                        bio11_GR_2100_126,
                        bio15_GR_2100_126,
                        bio16_GR_2100_126,
                        bio17_GR_2100_126
                        )
GR_clim_2040_585 <- stack(
                        bio1_GR_2040_585,
                        bio4_GR_2040_585,
                        bio10_GR_2040_585,
                        bio11_GR_2040_585,
                        bio15_GR_2040_585,
                        bio16_GR_2040_585,
                        bio17_GR_2040_585
                        )
GR_clim_2100_585 <- stack(
                        bio1_GR_2100_585,
                        bio4_GR_2100_585,
                        bio10_GR_2100_585,
                        bio11_GR_2100_585,
                        bio15_GR_2100_585,
                        bio16_GR_2100_585,
                        bio17_GR_2100_585
                        )

# Name layers of stacks
names(GR_clim_current) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_126) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(GR_clim_current)
plot(GR_clim_2040_126)
plot(GR_clim_2100_126)
plot(GR_clim_2040_585)
plot(GR_clim_2100_585)

#---------------------------------------------------------
# Future projections
# Do projection of current data
ch_SpheCoop_proj_current <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_current,
                        proj.name = "current",
                        binary.meth = "TSS",
                        do.stack = FALSE
)

ch_SpheCoop_proj_ensamble_current <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_current,
                        projection.output = ch_SpheCoop_proj_current,
                        proj.name = "ens_current",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)
# Do projection of 2040 ssp126
SpheCoop_proj_2040_126 <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126,
                        proj.name = "proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

SpheCoop_proj_ensamble_2040_126 <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126,
                        projection.output = SpheCoop_proj_2040_126,
                        proj.name = "proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
SpheCoop_proj_2100_126 <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126,
                        proj.name = "proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

SpheCoop_proj_ensamble_2100_126 <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126,
                        projection.output = SpheCoop_proj_2040_126,
                        proj.name = "proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
SpheCoop_proj_2040_585 <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585,
                        proj.name = "proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

SpheCoop_proj_ensamble_2040_585 <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585,
                        projection.output = SpheCoop_proj_2040_585,
                        proj.name = "proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
SpheCoop_proj_2100_585 <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585,
                        proj.name = "proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

SpheCoop_proj_ensamble_2100_585 <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585,
                        projection.output = SpheCoop_proj_2040_585,
                        proj.name = "proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
stk_SpheCoop_ef_current <- get_predictions(ch_SpheCoop_proj_ensamble_current)
plot(stk_SpheCoop_ef_current)
stk_SpheCoop_ef_2040_126 <- get_predictions(SpheCoop_proj_ensamble_2040_126)
plot(stk_SpheCoop_ef_2040_126)
stk_SpheCoop_ef_2100_126 <- get_predictions(SpheCoop_proj_ensamble_2100_126)
plot(stk_SpheCoop_ef_2100_126)
stk_SpheCoop_ef_2040_585 <- get_predictions(SpheCoop_proj_ensamble_2040_585)
plot(stk_SpheCoop_ef_2040_585)
stk_SpheCoop_ef_2100_585 <- get_predictions(SpheCoop_proj_ensamble_2100_585)
plot(stk_SpheCoop_ef_2100_585)

# Set it to only use EMca and EMwmean
subset(stk_SpheCoop_ef_current ,grep('EMca|EMwmean', names(stk_SpheCoop_ef_current )))
subset(stk_SpheCoop_ef_2040_126 ,grep('EMca|EMwmean', names(stk_SpheCoop_ef_2040_126 )))
subset(stk_SpheCoop_ef_2100_126 ,grep('EMca|EMwmean', names(stk_SpheCoop_ef_2100_126 )))
subset(stk_SpheCoop_ef_2040_585 ,grep('EMca|EMwmean', names(stk_SpheCoop_ef_2040_585 )))
subset(stk_SpheCoop_ef_2100_585 ,grep('EMca|EMwmean', names(stk_SpheCoop_ef_2100_585 )))

#---------------------------------------------------------
# Species range change setup
# Garden route
Sphe_coop_bin_proj_GR_current <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_126 <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_126 <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_585 <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_585 <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

#---------------------------------------------------------
# Garden route ssp126
SRC_current_2040_126 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_126
)
SRC_current_2040_126$Compt.By.Models

SRC_current_2100_126 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_126
)
SRC_current_2100_126$Compt.By.Models

SRC_2040_2100_126 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_126,
    Sphe_coop_bin_proj_GR_2100_126
)
SRC_2040_2100_126$Compt.By.Models

# Garden route ssp585
SRC_current_2040_585 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_585
)
SRC_current_2040_585$Compt.By.Models

SRC_current_2100_585 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_585
)
SRC_current_2100_585$Compt.By.Models

SRC_2040_2100_585 <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_585,
    Sphe_coop_bin_proj_GR_2100_585
)
SRC_2040_2100_585$Compt.By.Models

#---------------------------------------------------------
# Plot species range change

# Load libraries to change levelplot theme
install.packages("viridis")
library(viridis)

# Garden route ssp126
raster_2040_126 <- SRC_current_2040_126$Diff.By.Pixel
raster_2100_126 <- SRC_current_2100_126$Diff.By.Pixel
raster_future_126 <- SRC_2040_2100_126$Diff.By.Pixel

plot(raster_2040_126)
plot(raster_2100_126)
plot(raster_future_126)
layers_126<-c(raster_2040_126,raster_2100_126, raster_future_126)
summary(layers_126)
SpheCoop_src_map_126 <- stack(layers_126)
summary(SpheCoop_src_map_126)
names(SpheCoop_src_map_126) <- c("current_2040_126_ca", "current_2040_126_wm","current_2100_126_ca","current_2100_126_wm","2040_2100_126_ca","2040_2100_126_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(SpheCoop_src_map_126, main = 'SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


# Garden route 585
raster_2040_585 <- SRC_current_2040_585$Diff.By.Pixel
raster_2100_585 <- SRC_current_2100_585$Diff.By.Pixel
raster_future_585 <- SRC_2040_2100_585$Diff.By.Pixel

plot(raster_2040_585)
plot(raster_2100_585)
plot(raster_future_585)
layers_585<-c(raster_2040_585,raster_2100_585, raster_future_585)
summary(layers_585)
SpheCoop_src_map_585 <- stack(layers_585)
summary(SpheCoop_src_map_585)
names(SpheCoop_src_map_585) <- c("current_2040_585_ca", "current_2040_585_wm","current_2100_585_ca","current_2100_585_wm","2040_2100_585_ca","2040_2100_585_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(SpheCoop_src_map_585, main = 'SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )

#---------------------------------------------------------