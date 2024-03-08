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
setwd("C:\\Users\\megan\\neDrive - Stellenbosch University\\Chapter 2 - Tree ferns and climate change\\0015054-240202131308920\\20240304_Biomod")
data <- read.csv("20240219_Sampled by raster.csv")
head(data)
print(nrow(data))

#---------------------------------------------------------
# Load current climate data (global)
ch_bio1_raster <- raster("climate_data\\CURRENT_DATA\\bio1_clipped.tif")
ch_bio4_raster <- raster("climate_data\\CURRENT_DATA\\bio4_clipped.tif")
ch_bio10_raster <- raster("climate_data\\CURRENT_DATA\\bio10_clipped.tif")
ch_bio11_raster <- raster("climate_data\\CURRENT_DATA\\bio11_clipped.tif")
ch_bio15_raster <- raster("climate_data\\CURRENT_DATA\\bio15_clipped.tif")
ch_bio16_raster <- raster("climate_data\\CURRENT_DATA\\bio16_clipped.tif")
ch_bio17_raster <- raster("climate_data\\CURRENT_DATA\\bio17_clipped.tif")

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
# Caret tuning
BIOMOD_Tuning(
  ch_SpheCoop_data,
  bm.options = BIOMOD_ModelingOptions(),
  models = c("GLM", "GBM", "GAM", "CTA", "ANN", "SRE", "FDA", "MARS", "RF", "MAXENT"),
  metric.eval = "ROC",
  ctrl.train = NULL,
  ctrl.train.tuneLength = 30,
  ctrl.ANN = NULL,
  ctrl.CTA = NULL,
  ctrl.FDA = NULL,
  ctrl.GAM = NULL,
  ctrl.GBM = NULL,
  ctrl.GLM = NULL,
  ctrl.MARS = NULL,
  ctrl.RF = NULL,
  ANN.method = "avNNet",
  ANN.decay.tune = c(0.001, 0.01, 0.05, 0.1),
  ANN.size.tune = c(2, 4, 6, 8),
  ANN.maxit = 500,
  ANN.MaxNWts = 10 * (ncol(bm.format@data.env.var) + 1) + 10 + 1,
  MARS.method = "earth",
  GAM.method = "gam",
  GLM.method = "glmStepAIC",
  GLM.type = c("simple", "quadratic", "polynomial", "s_smoother"),
  GLM.interaction = c(0, 1),
  ME.cvmethod = NULL,
  ME.overlap = FALSE,
  ME.kfolds = 10,
  ME.n.bg = 10000,
  ME.env = NULL,
  ME.metric = "ROC",
  ME.clamp = TRUE,
  ME.parallel = FALSE,
  ME.numCores = NULL,
  RF.method = "rf",
  weights = NULL
)


#---------------------------------------------------------
# Caret tuned hyperparameters
ch_SpheCoop_model_options <- BIOMOD_ModelingOptions(
                        GLM = list( type = 'polynomial',
                                    interaction.level = 1,
                                    myFormula = Sphaeropteris.cooperi ~ AnnMeanTemp + I(AnnMeanTemp^2) + I(AnnMeanTemp^3) + TempSeason + I(TempSeason^2) + I(TempSeason^3) + MeanTempWarQ + I(MeanTempWarQ^2) + MeanTempColdQ + I(MeanTempColdQ^3) + PrecipSeason + I(PrecipSeason^2) + I(PrecipSeason^3) + PrecipWetQ + I(PrecipWetQ^2) + I(PrecipWetQ^3) + PrecipDryQ + I(PrecipDryQ^2) + AnnMeanTemp:TempSeason + AnnMeanTemp:MeanTempWarQ + AnnMeanTemp:MeanTempColdQ + AnnMeanTemp:PrecipSeason + AnnMeanTemp:PrecipDryQ + TempSeason:MeanTempWarQ + TempSeason:MeanTempColdQ + TempSeason:PrecipSeason +
                                                TempSeason:PrecipWetQ + TempSeason:PrecipDryQ + MeanTempWarQ:MeanTempColdQ + MeanTempWarQ:PrecipSeason + MeanTempWarQ:PrecipWetQ + MeanTempWarQ:PrecipDryQ + MeanTempColdQ:PrecipSeason + MeanTempColdQ:PrecipWetQ + MeanTempColdQ:PrecipDryQ + PrecipSeason:PrecipWetQ + PrecipSeason:PrecipDryQ + PrecipWetQ:PrecipDryQ,
                                    test = 'none',
                                    family = binomial(link = 'logit'),
                                    mustart = 0.5,
                                    control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
                        GBM = list( distribution = 'bernoulli',
                                    n.trees = 10000,
                                    interaction.depth = 8,
                                    n.minobsinnode = 5,
                                    shrinkage = 0.01,
                                    bag.fraction = 0.5,
                                    train.fraction = 1,
                                    cv.folds = 3,
                                    keep.data = FALSE,
                                    verbose = FALSE,
                                    perf.method = 'cv',
                                    n.cores = 1),
                        GAM = list( algo = 'GAM_mgcv',
                                    type = 's_smoother',
                                    k = -1,
                                    interaction.level = 0,
                                    myFormula = NULL,
                                    family = binomial(link = 'logit'),
                                    method = 'GCV.Cp',
                                    optimizer = c('outer','newton'),
                                    select = FALSE,
                                    knots = NULL,
                                    paraPen = NULL,
                                    control = list(nthreads = 1, ncv.threads = 1, irls.reg = 0
                                                , epsilon = 1e-07, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                                , rank.tol = 1.49011611938477e-08
                                                , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                                , optim = list(factr=1e+07)
                                                , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                                , idLinksBases = TRUE, scalePenalty = TRUE, efs.lspmax = 15, efs.tol = 0.1
                                                , keepData = FALSE, edge.correct = FALSE) ),
                        FDA = list( method = 'mars',
                                    add_args = list(degree = 2, nprune = 19)),
                        MARS = list( type = 'simple',
                                        interaction.level = 0,
                                        myFormula = NULL,
                                        nk = NULL,
                                        penalty = 2,
                                        thresh = 0.001,
                                        nprune = 15,
                                        pmethod = 'backward'),
                        MAXENT = list( path_to_maxent.jar = 'C:\\Users\\chris\\Desktop\\MegMSc\\maxent',
                                        memory_allocated = 512,
                                        initial_heap_size = NULL,
                                        max_heap_size = NULL,
                                        background_data_dir = 'default',
                                        maximumbackground = 'default',
                                        maximumiterations = 200,
                                        visible = FALSE,
                                        linear = TRUE,
                                        quadratic = FALSE,
                                        product = FALSE,
                                        threshold = FALSE,
                                        hinge = FALSE,
                                        lq2lqptthreshold = 80,
                                        l2lqthreshold = 10,
                                        hingethreshold = 15,
                                        beta_threshold = -1,
                                        beta_categorical = -1,
                                        beta_lqp = -1,
                                        beta_hinge = -1,
                                        betamultiplier = 2,
                                        defaultprevalence = 0.5),
                        RF = list( do.classif = TRUE,
                                    ntree = 500,
                                    mtry = 3,
                                    sampsize = NULL,
                                    nodesize = 5,
                                    maxnodes = NULL),
                        CTA = list( method = 'class',
                                    parms = 'default',
                                    cost = NULL,
                                    control = list(xval = 5, minbucket = 5, minsplit = 5
                                                    , cp = list(cp=0), maxdepth = list(maxdepth=8)) ),
                        ANN = list( NbCV = 5,
                                    size = 6,
                                    decay = 0.001,
                                    rang = 0.1,
                                    maxit = 500)
                    )
                        
#---------------------------------------------------------
# Building the individual models
ch_SpheCoop_models <- BIOMOD_Modeling(
                        bm.format = ch_SpheCoop_data,
                        models = c("GLM","GBM","GAM","ANN","CTA","FDA","RF",'MARS', 'MAXENT'),
                        bm.options = ch_SpheCoop_model_options,
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
write.csv(indiv_SpheCoop_scores, file = "TUNED-Model-evaluations.csv", row.names = FALSE)

#---------------------------------------------------------
# Varaiable importances
importances<-(SpheCoop_models_var_import <- get_variables_importance(ch_SpheCoop_models))
write.csv(importances, file = "TUNED-Model-importances.csv", row.names = FALSE)
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
write.csv(ensemble_score_csv, file = "TUNED-Model-ensemble-evaluations.csv", row.names = FALSE)

#=========================================================
# CURRENT DATA
#=========================================================
# Load Current GR data
bio1_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio1.tif")
bio4_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio4.tif")
bio10_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio10.tif")
bio11_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio11.tif")
bio15_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio15.tif")
bio16_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio16.tif")
bio17_GR_current <- raster("input_data\\climate_data\\CURRENT_DATA\\garden_route\\bio17.tif")

# Load Current REU data
bio1_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio1.tif")
bio4_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio4.tif")
bio10_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio10.tif")
bio11_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio11.tif")
bio15_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio15.tif")
bio16_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio16.tif")
bio17_REU_current <- raster("input_data\\climate_data\\CURRENT_DATA\\reunion\\bio17.tif")

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

REU_clim_current <- stack(
                        bio1_REU_current,
                        bio4_REU_current,
                        bio10_REU_current,
                        bio11_REU_current,
                        bio15_REU_current,
                        bio16_REU_current,
                        bio17_REU_current
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
names(REU_clim_current) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
# Do projection of current data
GR_SpheCoop_proj_current <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_current,
                        proj.name = "GR_current",
                        binary.meth = "TSS",
                        do.stack = FALSE
)
GR_SpheCoop_proj_ensamble_current <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_current,
                        projection.output = GR_SpheCoop_proj_current,
                        proj.name = "GR_ens_current",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

REU_SpheCoop_proj_current <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_current,
                        proj.name = "REU_current",
                        binary.meth = "TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_current <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_current,
                        projection.output = REU_SpheCoop_proj_current,
                        proj.name = "REU_ens_current",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)
#=========================================================
# GFDL-ESM4
#=========================================================
# Garden Route
#---------------------------------------------------------
# Load future data ssp126
bio1_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio1_gfdlesm4_ssp126_2041_2070_GR.tif")
bio4_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio4_gfdlesm4_ssp126_2041_2070_GR.tif")
bio10_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio10_gfdlesm4_ssp126_2041_2070_GR.tif")
bio11_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio11_gfdlesm4_ssp126_2041_2070_GR.tif")
bio15_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio15_gfdlesm4_ssp126_2041_2070_GR.tif")
bio16_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio16_gfdlesm4_ssp126_2041_2070_GR.tif")
bio17_GR_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2041-2070\\bio17_gfdlesm4_ssp126_2041_2070_GR.tif")

bio1_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio1_gfdlesm4_ssp126_2071_2100_GR.tif")
bio4_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio4_gfdlesm4_ssp126_2071_2100_GR.tif")
bio10_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio10_gfdlesm4_ssp126_2071_2100_GR.tif")
bio11_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio11_gfdlesm4_ssp126_2071_2100_GR.tif")
bio15_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio15_gfdlesm4_ssp126_2071_2100_GR.tif")
bio16_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio16_gfdlesm4_ssp126_2071_2100_GR.tif")
bio17_GR_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp126\\2071-2100\\bio17_gfdlesm4_ssp126_2071_2100_GR.tif")

# Load future data ssp585
bio1_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio1_gfdlesm4_ssp585_2040_2071_GR.tif")
bio4_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio4_gfdlesm4_ssp585_2040_2071_GR.tif")
bio10_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio10_gfdlesm4_ssp585_2040_2071_GR.tif")
bio11_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio11_gfdlesm4_ssp585_2040_2071_GR.tif")
bio15_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio15_gfdlesm4_ssp585_2040_2071_GR.tif")
bio16_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio16_gfdlesm4_ssp585_2040_2071_GR.tif")
bio17_GR_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2041-2070\\bio17_gfdlesm4_ssp585_2040_2071_GR.tif")

bio1_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio1_gfdlesm4_ssp585_2071_2100_GR.tif")
bio4_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio4_gfdlesm4_ssp585_2071_2100_GR.tif")
bio10_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio10_gfdlesm4_ssp585_2071_2100_GR.tif")
bio11_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio11_gfdlesm4_ssp585_2071_2100_GR.tif")
bio15_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio15_gfdlesm4_ssp585_2071_2100_GR.tif")
bio16_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio16_gfdlesm4_ssp585_2071_2100_GR.tif")
bio17_GR_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Garden route\\ssp585\\2071-2100\\bio17_gfdlesm4_ssp585_2071_2100_GR.tif")

# Create raster stacks
GR_clim_2040_126_GFDL <- stack(
                        bio1_GR_2040_126_GFDL,
                        bio4_GR_2040_126_GFDL,
                        bio10_GR_2040_126_GFDL,
                        bio11_GR_2040_126_GFDL,
                        bio15_GR_2040_126_GFDL,
                        bio16_GR_2040_126_GFDL,
                        bio17_GR_2040_126_GFDL
                        )
GR_clim_2100_126_GFDL <- stack(
                        bio1_GR_2100_126_GFDL,
                        bio4_GR_2100_126_GFDL,
                        bio10_GR_2100_126_GFDL,
                        bio11_GR_2100_126_GFDL,
                        bio15_GR_2100_126_GFDL,
                        bio16_GR_2100_126_GFDL,
                        bio17_GR_2100_126_GFDL
                        )
GR_clim_2040_585_GFDL <- stack(
                        bio1_GR_2040_585_GFDL,
                        bio4_GR_2040_585_GFDL,
                        bio10_GR_2040_585_GFDL,
                        bio11_GR_2040_585_GFDL,
                        bio15_GR_2040_585_GFDL,
                        bio16_GR_2040_585_GFDL,
                        bio17_GR_2040_585_GFDL
                        )
GR_clim_2100_585_GFDL <- stack(
                        bio1_GR_2100_585_GFDL,
                        bio4_GR_2100_585_GFDL,
                        bio10_GR_2100_585_GFDL,
                        bio11_GR_2100_585_GFDL,
                        bio15_GR_2100_585_GFDL,
                        bio16_GR_2100_585_GFDL,
                        bio17_GR_2100_585_GFDL
                        )

# Name layers of stacks
names(GR_clim_2040_126_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(GR_clim_2040_126_GFDL)
plot(GR_clim_2100_126_GFDL)
plot(GR_clim_2040_585_GFDL)
plot(GR_clim_2100_585_GFDL)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
GR_SpheCoop_proj_2040_126_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126_GFDL,
                        proj.name = "GR_GFDL_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_126_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126_GFDL,
                        projection.output = GR_SpheCoop_proj_2040_126_GFDL,
                        proj.name = "GR_GFDL_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
GR_SpheCoop_proj_2100_126_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126_GFDL,
                        proj.name = "GR_GFDL_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_126_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126_GFDL,
                        projection.output = GR_SpheCoop_proj_2040_126_GFDL,
                        proj.name = "GR_GFDL_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
GR_SpheCoop_proj_2040_585_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585_GFDL,
                        proj.name = "GR_GFDL_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_585_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585_GFDL,
                        projection.output = GR_SpheCoop_proj_2040_585_GFDL,
                        proj.name = "GR_GFDL_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
GR_SpheCoop_proj_2100_585_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585_GFDL,
                        proj.name = "GR_GFDL_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_585_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585_GFDL,
                        projection.output = GR_SpheCoop_proj_2040_585_GFDL,
                        proj.name = "GR_GFDL_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
GR_stk_SpheCoop_ef_current_GFDL <- get_predictions(GR_SpheCoop_proj_ensamble_current)
plot(GR_stk_SpheCoop_ef_current_GFDL)
GR_stk_SpheCoop_ef_2040_126_GFDL <- get_predictions(GR_SpheCoop_proj_ensamble_2040_126_GFDL)
plot(GR_stk_SpheCoop_ef_2040_126_GFDL)
GR_stk_SpheCoop_ef_2100_126_GFDL <- get_predictions(GR_SpheCoop_proj_ensamble_2100_126_GFDL)
plot(GR_stk_SpheCoop_ef_2100_126_GFDL)
GR_stk_SpheCoop_ef_2040_585_GFDL <- get_predictions(GR_SpheCoop_proj_ensamble_2040_585_GFDL)
plot(GR_stk_SpheCoop_ef_2040_585_GFDL)
GR_stk_SpheCoop_ef_2100_585_GFDL <- get_predictions(GR_SpheCoop_proj_ensamble_2100_585_GFDL)
plot(GR_stk_SpheCoop_ef_2100_585_GFDL)

# Set it to only use EMca and EMwmean
subset(GR_stk_SpheCoop_ef_current_GFDL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_current_GFDL )))
subset(GR_stk_SpheCoop_ef_2040_126_GFDL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_126_GFDL )))
subset(GR_stk_SpheCoop_ef_2100_126_GFDL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_126_GFDL )))
subset(GR_stk_SpheCoop_ef_2040_585_GFDL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_585_GFDL )))
subset(GR_stk_SpheCoop_ef_2100_585_GFDL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_585_GFDL )))

#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
# Loading raster data
# Load future data ssp126
bio1_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio1_gfdlesm4_ssp126_2041_2070_REU.tif")
bio4_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio4_gfdlesm4_ssp126_2041_2070_REU.tif")
bio10_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio10_gfdlesm4_ssp126_2041_2070_REU.tif")
bio11_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio11_gfdlesm4_ssp126_2041_2070_REU.tif")
bio15_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio15_gfdlesm4_ssp126_2041_2070_REU.tif")
bio16_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio16_gfdlesm4_ssp126_2041_2070_REU.tif")
bio17_REU_2040_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2041-2070\\bio17_gfdlesm4_ssp126_2041_2070_REU.tif")

bio1_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio1_gfdlesm4_ssp126_2071_2100_REU.tif")
bio4_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio4_gfdlesm4_ssp126_2071_2100_REU.tif")
bio10_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio10_gfdlesm4_ssp126_2071_2100_REU.tif")
bio11_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio11_gfdlesm4_ssp126_2071_2100_REU.tif")
bio15_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio15_gfdlesm4_ssp126_2071_2100_REU.tif")
bio16_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio16_gfdlesm4_ssp126_2071_2100_REU.tif")
bio17_REU_2100_126_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp126\\2071-2100\\bio17_gfdlesm4_ssp126_2071_2100_REU.tif")

# Load future data ssp585
bio1_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio1_gfdlesm4_ssp585_2040_2071_REU.tif")
bio4_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio4_gfdlesm4_ssp585_2040_2071_REU.tif")
bio10_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio10_gfdlesm4_ssp585_2040_2071_REU.tif")
bio11_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio11_gfdlesm4_ssp585_2040_2071_REU.tif")
bio15_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio15_gfdlesm4_ssp585_2040_2071_REU.tif")
bio16_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio16_gfdlesm4_ssp585_2040_2071_REU.tif")
bio17_REU_2040_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2041-2070\\bio17_gfdlesm4_ssp585_2040_2071_REU.tif")

bio1_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio1_gfdlesm4_ssp585_2071_2100_REU.tif")
bio4_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio4_gfdlesm4_ssp585_2071_2100_REU.tif")
bio10_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio10_gfdlesm4_ssp585_2071_2100_REU.tif")
bio11_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio11_gfdlesm4_ssp585_2071_2100_REU.tif")
bio15_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio15_gfdlesm4_ssp585_2071_2100_REU.tif")
bio16_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio16_gfdlesm4_ssp585_2071_2100_REU.tif")
bio17_REU_2100_585_GFDL <- raster("input_data\\climate_data\\GFDL-ESM4\\Reunion\\ssp585\\2071-2100\\bio17_gfdlesm4_ssp585_2071_2100_REU.tif")

# Create stacks of climatic data
REU_clim_2040_126_GFDL <- stack(
                        bio1_REU_2040_126_GFDL,
                        bio4_REU_2040_126_GFDL,
                        bio10_REU_2040_126_GFDL,
                        bio11_REU_2040_126_GFDL,
                        bio15_REU_2040_126_GFDL,
                        bio16_REU_2040_126_GFDL,
                        bio17_REU_2040_126_GFDL
                        )
REU_clim_2100_126_GFDL <- stack(
                        bio1_REU_2100_126_GFDL,
                        bio4_REU_2100_126_GFDL,
                        bio10_REU_2100_126_GFDL,
                        bio11_REU_2100_126_GFDL,
                        bio15_REU_2100_126_GFDL,
                        bio16_REU_2100_126_GFDL,
                        bio17_REU_2100_126_GFDL
                        )
REU_clim_2040_585_GFDL <- stack(
                        bio1_REU_2040_585_GFDL,
                        bio4_REU_2040_585_GFDL,
                        bio10_REU_2040_585_GFDL,
                        bio11_REU_2040_585_GFDL,
                        bio15_REU_2040_585_GFDL,
                        bio16_REU_2040_585_GFDL,
                        bio17_REU_2040_585_GFDL
                        )
REU_clim_2100_585_GFDL <- stack(
                        bio1_REU_2100_585_GFDL,
                        bio4_REU_2100_585_GFDL,
                        bio10_REU_2100_585_GFDL,
                        bio11_REU_2100_585_GFDL,
                        bio15_REU_2100_585_GFDL,
                        bio16_REU_2100_585_GFDL,
                        bio17_REU_2100_585_GFDL
                        )

# Name layers of stacks
names(REU_clim_2040_126_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_126_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2040_585_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_585_GFDL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(REU_clim_current)
plot(REU_clim_2040_126_GFDL)
plot(REU_clim_2100_126_GFDL)
plot(REU_clim_2040_585_GFDL)
plot(REU_clim_2100_585_GFDL)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
REU_SpheCoop_proj_2040_126_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_126_GFDL,
                        proj.name = "REU_GFDL_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_126_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_126_GFDL,
                        projection.output = REU_SpheCoop_proj_2040_126_GFDL,
                        proj.name = "REU_GFDL_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
REU_SpheCoop_proj_2100_126_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_126_GFDL,
                        proj.name = "REU_GFDL_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_126_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_126_GFDL,
                        projection.output = REU_SpheCoop_proj_2100_126_GFDL,
                        proj.name = "REU_GFDL_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
REU_SpheCoop_proj_2040_585_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_585_GFDL,
                        proj.name = "REU_GFDL_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_585_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_585_GFDL,
                        projection.output = REU_SpheCoop_proj_2040_585_GFDL,
                        proj.name = "REU_GFDL_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
REU_SpheCoop_proj_2100_585_GFDL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_585_GFDL,
                        proj.name = "REU_GFDL_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_585_GFDL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_585_GFDL,
                        projection.output = REU_SpheCoop_proj_2100_585_GFDL,
                        proj.name = "REU_GFDL_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
REU_stk_SpheCoop_ef_current_GFDL <- get_predictions(REU_SpheCoop_proj_ensamble_current)
plot(REU_stk_SpheCoop_ef_current_GFDL)
REU_stk_SpheCoop_ef_2040_126_GFDL <- get_predictions(REU_SpheCoop_proj_ensamble_2040_126_GFDL)
plot(REU_stk_SpheCoop_ef_2040_126_GFDL)
REU_stk_SpheCoop_ef_2100_126_GFDL <- get_predictions(REU_SpheCoop_proj_ensamble_2100_126_GFDL)
plot(REU_stk_SpheCoop_ef_2100_126_GFDL)
REU_stk_SpheCoop_ef_2040_585_GFDL <- get_predictions(REU_SpheCoop_proj_ensamble_2040_585_GFDL)
plot(REU_stk_SpheCoop_ef_2040_585_GFDL)
REU_stk_SpheCoop_ef_2100_585_GFDL <- get_predictions(REU_SpheCoop_proj_ensamble_2100_585_GFDL)
plot(REU_stk_SpheCoop_ef_2100_585_GFDL)

# Set it to only use EMca and EMwmean
subset(REU_stk_SpheCoop_ef_current_GFDL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_current_GFDL )))
subset(REU_stk_SpheCoop_ef_2040_126_GFDL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_126_GFDL )))
subset(REU_stk_SpheCoop_ef_2100_126_GFDL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_126_GFDL )))
subset(REU_stk_SpheCoop_ef_2040_585_GFDL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_585_GFDL )))
subset(REU_stk_SpheCoop_ef_2100_585_GFDL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_585_GFDL )))

#=========================================================
# IPSL-CM6A
#=========================================================
# Garden Route
#---------------------------------------------------------
# Load future data ssp126
bio1_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio1_ipslcm6a_ssp126_2041_2070_GR.tif")
bio4_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio4_ipslcm6a_ssp126_2041_2070_GR.tif")
bio10_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio10_ipslcm6a_ssp126_2041_2070_GR.tif")
bio11_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio11_ipslcm6a_ssp126_2041_2070_GR.tif")
bio15_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio15_ipslcm6a_ssp126_2041_2070_GR.tif")
bio16_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio16_ipslcm6a_ssp126_2041_2070_GR.tif")
bio17_GR_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2041-2070\\bio17_ipslcm6a_ssp126_2041_2070_GR.tif")

bio1_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio1_ipslcm6a_ssp126_2071_2100_GR.tif")
bio4_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio4_ipslcm6a_ssp126_2071_2100_GR.tif")
bio10_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio10_ipslcm6a_ssp126_2071_2100_GR.tif")
bio11_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio11_ipslcm6a_ssp126_2071_2100_GR.tif")
bio15_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio15_ipslcm6a_ssp126_2071_2100_GR.tif")
bio16_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio16_ipslcm6a_ssp126_2071_2100_GR.tif")
bio17_GR_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp126\\2071-2100\\bio17_ipslcm6a_ssp126_2071_2100_GR.tif")

# Load future data ssp585
bio1_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio1_ipslcm6a_ssp585_2041_2070_GR.tif")
bio4_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio4_ipslcm6a_ssp585_2041_2070_GR.tif")
bio10_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio10_ipslcm6a_ssp585_2041_2070_GR.tif")
bio11_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio11_ipslcm6a_ssp585_2041_2070_GR.tif")
bio15_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio15_ipslcm6a_ssp585_2041_2070_GR.tif")
bio16_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio16_ipslcm6a_ssp585_2041_2070_GR.tif")
bio17_GR_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2041-2070\\bio17_ipslcm6a_ssp585_2041_2070_GR.tif")

bio1_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio1_ipslcm6a_ssp585_2071_2100_GR.tif")
bio4_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio4_ipslcm6a_ssp585_2071_2100_GR.tif")
bio10_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio10_ipslcm6a_ssp585_2071_2100_GR.tif")
bio11_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio11_ipslcm6a_ssp585_2071_2100_GR.tif")
bio15_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio15_ipslcm6a_ssp585_2071_2100_GR.tif")
bio16_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio16_ipslcm6a_ssp585_2071_2100_GR.tif")
bio17_GR_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Garden route\\ssp585\\2071-2100\\bio17_ipslcm6a_ssp585_2071_2100_GR.tif")

# Create raster stacks
GR_clim_2040_126_IPSL <- stack(
                        bio1_GR_2040_126_IPSL,
                        bio4_GR_2040_126_IPSL,
                        bio10_GR_2040_126_IPSL,
                        bio11_GR_2040_126_IPSL,
                        bio15_GR_2040_126_IPSL,
                        bio16_GR_2040_126_IPSL,
                        bio17_GR_2040_126_IPSL
                        )
GR_clim_2100_126_IPSL <- stack(
                        bio1_GR_2100_126_IPSL,
                        bio4_GR_2100_126_IPSL,
                        bio10_GR_2100_126_IPSL,
                        bio11_GR_2100_126_IPSL,
                        bio15_GR_2100_126_IPSL,
                        bio16_GR_2100_126_IPSL,
                        bio17_GR_2100_126_IPSL
                        )
GR_clim_2040_585_IPSL <- stack(
                        bio1_GR_2040_585_IPSL,
                        bio4_GR_2040_585_IPSL,
                        bio10_GR_2040_585_IPSL,
                        bio11_GR_2040_585_IPSL,
                        bio15_GR_2040_585_IPSL,
                        bio16_GR_2040_585_IPSL,
                        bio17_GR_2040_585_IPSL
                        )
GR_clim_2100_585_IPSL <- stack(
                        bio1_GR_2100_585_IPSL,
                        bio4_GR_2100_585_IPSL,
                        bio10_GR_2100_585_IPSL,
                        bio11_GR_2100_585_IPSL,
                        bio15_GR_2100_585_IPSL,
                        bio16_GR_2100_585_IPSL,
                        bio17_GR_2100_585_IPSL
                        )

# Name layers of stacks
names(GR_clim_2040_126_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(GR_clim_2040_126_IPSL)
plot(GR_clim_2100_126_IPSL)
plot(GR_clim_2040_585_IPSL)
plot(GR_clim_2100_585_IPSL)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
GR_SpheCoop_proj_2040_126_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126_IPSL,
                        proj.name = "GR_IPSL_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_126_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126_IPSL,
                        projection.output = GR_SpheCoop_proj_2040_126_IPSL,
                        proj.name = "GR_IPSL_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
GR_SpheCoop_proj_2100_126_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126_IPSL,
                        proj.name = "GR_IPSL_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_126_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126_IPSL,
                        projection.output = GR_SpheCoop_proj_2040_126_IPSL,
                        proj.name = "GR_IPSL_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
GR_SpheCoop_proj_2040_585_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585_IPSL,
                        proj.name = "GR_IPSL_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_585_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585_IPSL,
                        projection.output = GR_SpheCoop_proj_2040_585_IPSL,
                        proj.name = "GR_IPSL_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
GR_SpheCoop_proj_2100_585_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585_IPSL,
                        proj.name = "GR_IPSL_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_585_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585_IPSL,
                        projection.output = GR_SpheCoop_proj_2040_585_IPSL,
                        proj.name = "GR_IPSL_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
GR_stk_SpheCoop_ef_current_IPSL <- get_predictions(GR_SpheCoop_proj_ensamble_current)
plot(GR_stk_SpheCoop_ef_current_IPSL)
GR_stk_SpheCoop_ef_2040_126_IPSL <- get_predictions(GR_SpheCoop_proj_ensamble_2040_126_IPSL)
plot(GR_stk_SpheCoop_ef_2040_126_IPSL)
GR_stk_SpheCoop_ef_2100_126_IPSL <- get_predictions(GR_SpheCoop_proj_ensamble_2100_126_IPSL)
plot(GR_stk_SpheCoop_ef_2100_126_IPSL)
GR_stk_SpheCoop_ef_2040_585_IPSL <- get_predictions(GR_SpheCoop_proj_ensamble_2040_585_IPSL)
plot(GR_stk_SpheCoop_ef_2040_585_IPSL)
GR_stk_SpheCoop_ef_2100_585_IPSL <- get_predictions(GR_SpheCoop_proj_ensamble_2100_585_IPSL)
plot(GR_stk_SpheCoop_ef_2100_585_IPSL)

# Set it to only use EMca and EMwmean
subset(GR_stk_SpheCoop_ef_current_IPSL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_current_IPSL )))
subset(GR_stk_SpheCoop_ef_2040_126_IPSL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_126_IPSL )))
subset(GR_stk_SpheCoop_ef_2100_126_IPSL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_126_IPSL )))
subset(GR_stk_SpheCoop_ef_2040_585_IPSL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_585_IPSL )))
subset(GR_stk_SpheCoop_ef_2100_585_IPSL ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_585_IPSL )))


#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
# Loading raster data
# Load future data ssp126
bio1_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio1_ipslcm6a_ssp126_2041_2070_REU.tif")
bio4_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio4_ipslcm6a_ssp126_2041_2070_REU.tif")
bio10_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio10_ipslcm6a_ssp126_2041_2070_REU.tif")
bio11_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio11_ipslcm6a_ssp126_2041_2070_REU.tif")
bio15_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio15_ipslcm6a_ssp126_2041_2070_REU.tif")
bio16_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio16_ipslcm6a_ssp126_2041_2070_REU.tif")
bio17_REU_2040_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2041-2070\\bio17_ipslcm6a_ssp126_2041_2070_REU.tif")

bio1_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio1_ipslcm6a_ssp126_2071_2100_REU.tif")
bio4_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio4_ipslcm6a_ssp126_2071_2100_REU.tif")
bio10_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio10_ipslcm6a_ssp126_2071_2100_REU.tif")
bio11_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio11_ipslcm6a_ssp126_2071_2100_REU.tif")
bio15_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio15_ipslcm6a_ssp126_2071_2100_REU.tif")
bio16_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio16_ipslcm6a_ssp126_2071_2100_REU.tif")
bio17_REU_2100_126_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp126\\2071-2100\\bio17_ipslcm6a_ssp126_2071_2100_REU.tif")

# Load future data ssp585
bio1_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio1_ipslcm6a_ssp585_2041_2070_REU.tif")
bio4_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio4_ipslcm6a_ssp585_2041_2070_REU.tif")
bio10_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio10_ipslcm6a_ssp585_2041_2070_REU.tif")
bio11_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio11_ipslcm6a_ssp585_2041_2070_REU.tif")
bio15_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio15_ipslcm6a_ssp585_2041_2070_REU.tif")
bio16_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio16_ipslcm6a_ssp585_2041_2070_REU.tif")
bio17_REU_2040_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2041-2070\\bio17_ipslcm6a_ssp585_2041_2070_REU.tif")

bio1_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio1_ipslcm6a_ssp585_2071_2100_REU.tif")
bio4_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio4_ipslcm6a_ssp585_2071_2100_REU.tif")
bio10_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio10_ipslcm6a_ssp585_2071_2100_REU.tif")
bio11_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio11_ipslcm6a_ssp585_2071_2100_REU.tif")
bio15_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio15_ipslcm6a_ssp585_2071_2100_REU.tif")
bio16_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio16_ipslcm6a_ssp585_2071_2100_REU.tif")
bio17_REU_2100_585_IPSL <- raster("input_data\\climate_data\\IPSL-CM6A\\Reunion\\ssp585\\2071-2100\\bio17_ipslcm6a_ssp585_2071_2100_REU.tif")

# Create stacks of climatic data
REU_clim_2040_126_IPSL <- stack(
                        bio1_REU_2040_126_IPSL,
                        bio4_REU_2040_126_IPSL,
                        bio10_REU_2040_126_IPSL,
                        bio11_REU_2040_126_IPSL,
                        bio15_REU_2040_126_IPSL,
                        bio16_REU_2040_126_IPSL,
                        bio17_REU_2040_126_IPSL
                        )
REU_clim_2100_126_IPSL <- stack(
                        bio1_REU_2100_126_IPSL,
                        bio4_REU_2100_126_IPSL,
                        bio10_REU_2100_126_IPSL,
                        bio11_REU_2100_126_IPSL,
                        bio15_REU_2100_126_IPSL,
                        bio16_REU_2100_126_IPSL,
                        bio17_REU_2100_126_IPSL
                        )
REU_clim_2040_585_IPSL <- stack(
                        bio1_REU_2040_585_IPSL,
                        bio4_REU_2040_585_IPSL,
                        bio10_REU_2040_585_IPSL,
                        bio11_REU_2040_585_IPSL,
                        bio15_REU_2040_585_IPSL,
                        bio16_REU_2040_585_IPSL,
                        bio17_REU_2040_585_IPSL
                        )
REU_clim_2100_585_IPSL <- stack(
                        bio1_REU_2100_585_IPSL,
                        bio4_REU_2100_585_IPSL,
                        bio10_REU_2100_585_IPSL,
                        bio11_REU_2100_585_IPSL,
                        bio15_REU_2100_585_IPSL,
                        bio16_REU_2100_585_IPSL,
                        bio17_REU_2100_585_IPSL
                        )

# Name layers of stacks
names(REU_clim_2040_126_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_126_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2040_585_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_585_IPSL) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(REU_clim_current)
plot(REU_clim_2040_126_IPSL)
plot(REU_clim_2100_126_IPSL)
plot(REU_clim_2040_585_IPSL)
plot(REU_clim_2100_585_IPSL)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
REU_SpheCoop_proj_2040_126_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_126_IPSL,
                        proj.name = "REU_IPSL_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_126_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_126_IPSL,
                        projection.output = REU_SpheCoop_proj_2040_126_IPSL,
                        proj.name = "REU_IPSL_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
REU_SpheCoop_proj_2100_126_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_126_IPSL,
                        proj.name = "REU_IPSL_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_126_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_126_IPSL,
                        projection.output = REU_SpheCoop_proj_2100_126_IPSL,
                        proj.name = "REU_IPSL_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
REU_SpheCoop_proj_2040_585_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_585_IPSL,
                        proj.name = "REU_IPSL_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_585_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_585_IPSL,
                        projection.output = REU_SpheCoop_proj_2040_585_IPSL,
                        proj.name = "REU_IPSL_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
REU_SpheCoop_proj_2100_585_IPSL <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_585_IPSL,
                        proj.name = "REU_IPSL_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_585_IPSL <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_585_IPSL,
                        projection.output = REU_SpheCoop_proj_2100_585_IPSL,
                        proj.name = "REU_IPSL_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
REU_stk_SpheCoop_ef_current_IPSL <- get_predictions(REU_SpheCoop_proj_ensamble_current)
plot(REU_stk_SpheCoop_ef_current_IPSL)
REU_stk_SpheCoop_ef_2040_126_IPSL <- get_predictions(REU_SpheCoop_proj_ensamble_2040_126_IPSL)
plot(REU_stk_SpheCoop_ef_2040_126_IPSL)
REU_stk_SpheCoop_ef_2100_126_IPSL <- get_predictions(REU_SpheCoop_proj_ensamble_2100_126_IPSL)
plot(REU_stk_SpheCoop_ef_2100_126_IPSL)
REU_stk_SpheCoop_ef_2040_585_IPSL <- get_predictions(REU_SpheCoop_proj_ensamble_2040_585_IPSL)
plot(REU_stk_SpheCoop_ef_2040_585_IPSL)
REU_stk_SpheCoop_ef_2100_585_IPSL <- get_predictions(REU_SpheCoop_proj_ensamble_2100_585_IPSL)
plot(REU_stk_SpheCoop_ef_2100_585_IPSL)

# Set it to only use EMca and EMwmean
subset(REU_stk_SpheCoop_ef_current_IPSL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_current_IPSL )))
subset(REU_stk_SpheCoop_ef_2040_126_IPSL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_126_IPSL )))
subset(REU_stk_SpheCoop_ef_2100_126_IPSL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_126_IPSL )))
subset(REU_stk_SpheCoop_ef_2040_585_IPSL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_585_IPSL )))
subset(REU_stk_SpheCoop_ef_2100_585_IPSL ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_585_IPSL )))

#=========================================================
# MRI
#=========================================================
# Garden Route
#---------------------------------------------------------
# Load future data ssp126
bio1_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio1_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio4_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio4_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio10_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio10_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio11_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio11_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio15_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio15_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio16_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio16_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio17_GR_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2041-2070\\bio17_2041-2070_mri-esm2-0_ssp126_V.2.1_GR.tif")

bio1_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio1_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio4_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio4_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio10_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio10_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio11_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio11_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio15_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio15_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio16_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio16_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")
bio17_GR_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp126\\2071-2100\\bio17_2071-2100_mri-esm2-0_ssp126_V.2.1_GR.tif")

# Load future data ssp585
bio1_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio1_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio4_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio4_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio10_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio10_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio11_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio11_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio15_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio15_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio16_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio16_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio17_GR_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2041-2070\\bio17_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif")

bio1_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio1_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio4_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio4_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio10_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio10_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio11_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio11_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio15_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio15_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio16_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio16_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")
bio17_GR_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Garden route\\ssp585\\2071-2100\\bio17_2071-2100_mri-esm2-0_ssp585_V.2.1_GR.tif")

# Create raster stacks
GR_clim_2040_126_MRI <- stack(
                        bio1_GR_2040_126_MRI,
                        bio4_GR_2040_126_MRI,
                        bio10_GR_2040_126_MRI,
                        bio11_GR_2040_126_MRI,
                        bio15_GR_2040_126_MRI,
                        bio16_GR_2040_126_MRI,
                        bio17_GR_2040_126_MRI
                        )
GR_clim_2100_126_MRI <- stack(
                        bio1_GR_2100_126_MRI,
                        bio4_GR_2100_126_MRI,
                        bio10_GR_2100_126_MRI,
                        bio11_GR_2100_126_MRI,
                        bio15_GR_2100_126_MRI,
                        bio16_GR_2100_126_MRI,
                        bio17_GR_2100_126_MRI
                        )
GR_clim_2040_585_MRI <- stack(
                        bio1_GR_2040_585_MRI,
                        bio4_GR_2040_585_MRI,
                        bio10_GR_2040_585_MRI,
                        bio11_GR_2040_585_MRI,
                        bio15_GR_2040_585_MRI,
                        bio16_GR_2040_585_MRI,
                        bio17_GR_2040_585_MRI
                        )
GR_clim_2100_585_MRI <- stack(
                        bio1_GR_2100_585_MRI,
                        bio4_GR_2100_585_MRI,
                        bio10_GR_2100_585_MRI,
                        bio11_GR_2100_585_MRI,
                        bio15_GR_2100_585_MRI,
                        bio16_GR_2100_585_MRI,
                        bio17_GR_2100_585_MRI
                        )

# Name layers of stacks
names(GR_clim_2040_126_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(GR_clim_2040_126_MRI)
plot(GR_clim_2100_126_MRI)
plot(GR_clim_2040_585_MRI)
plot(GR_clim_2100_585_MRI)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
GR_SpheCoop_proj_2040_126_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126_MRI,
                        proj.name = "GR_MRI_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_126_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126_MRI,
                        projection.output = GR_SpheCoop_proj_2040_126_MRI,
                        proj.name = "GR_MRI_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
GR_SpheCoop_proj_2100_126_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126_MRI,
                        proj.name = "GR_MRI_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_126_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126_MRI,
                        projection.output = GR_SpheCoop_proj_2040_126_MRI,
                        proj.name = "GR_MRI_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
GR_SpheCoop_proj_2040_585_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585_MRI,
                        proj.name = "GR_MRI_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_585_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585_MRI,
                        projection.output = GR_SpheCoop_proj_2040_585_MRI,
                        proj.name = "GR_MRI_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
GR_SpheCoop_proj_2100_585_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585_MRI,
                        proj.name = "GR_MRI_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_585_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585_MRI,
                        projection.output = GR_SpheCoop_proj_2040_585_MRI,
                        proj.name = "GR_MRI_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
GR_stk_SpheCoop_ef_current_MRI <- get_predictions(GR_SpheCoop_proj_ensamble_current)
plot(GR_stk_SpheCoop_ef_current_MRI)
GR_stk_SpheCoop_ef_2040_126_MRI <- get_predictions(GR_SpheCoop_proj_ensamble_2040_126_MRI)
plot(GR_stk_SpheCoop_ef_2040_126_MRI)
GR_stk_SpheCoop_ef_2100_126_MRI <- get_predictions(GR_SpheCoop_proj_ensamble_2100_126_MRI)
plot(GR_stk_SpheCoop_ef_2100_126_MRI)
GR_stk_SpheCoop_ef_2040_585_MRI <- get_predictions(GR_SpheCoop_proj_ensamble_2040_585_MRI)
plot(GR_stk_SpheCoop_ef_2040_585_MRI)
GR_stk_SpheCoop_ef_2100_585_MRI <- get_predictions(GR_SpheCoop_proj_ensamble_2100_585_MRI)
plot(GR_stk_SpheCoop_ef_2100_585_MRI)

# Set it to only use EMca and EMwmean
subset(GR_stk_SpheCoop_ef_current_MRI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_current_MRI )))
subset(GR_stk_SpheCoop_ef_2040_126_MRI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_126_MRI )))
subset(GR_stk_SpheCoop_ef_2100_126_MRI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_126_MRI )))
subset(GR_stk_SpheCoop_ef_2040_585_MRI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_585_MRI )))
subset(GR_stk_SpheCoop_ef_2100_585_MRI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_585_MRI )))


#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
# Loading raster data
# Load future data ssp126
bio1_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio4_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio4_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio10_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio10_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio11_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio11_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio15_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio15_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio16_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio16_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")
bio17_REU_2040_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio17_2041-2070_mri-esm2-0_ssp126_V.2.1_R.tif")

bio1_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio1_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio4_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio4_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio10_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio10_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio11_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio11_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio15_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio15_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio16_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio16_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")
bio17_REU_2100_126_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio17_2071-2100_mri-esm2-0_ssp126_V.2.1_R.tif")

# Load future data ssp585
bio1_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio4_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio4_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio10_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio10_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio11_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio11_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio15_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio15_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio16_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio16_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")
bio17_REU_2040_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio17_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif")

bio1_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio1_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio4_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio4_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio10_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio10_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio11_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio11_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio15_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio15_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio16_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio16_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")
bio17_REU_2100_585_MRI <- raster("input_data\\climate_data\\MRI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio17_2071-2100_mri-esm2-0_ssp585_V.2.1_R.tif")

# Create stacks of climatic data
REU_clim_2040_126_MRI <- stack(
                        bio1_REU_2040_126_MRI,
                        bio4_REU_2040_126_MRI,
                        bio10_REU_2040_126_MRI,
                        bio11_REU_2040_126_MRI,
                        bio15_REU_2040_126_MRI,
                        bio16_REU_2040_126_MRI,
                        bio17_REU_2040_126_MRI
                        )
REU_clim_2100_126_MRI <- stack(
                        bio1_REU_2100_126_MRI,
                        bio4_REU_2100_126_MRI,
                        bio10_REU_2100_126_MRI,
                        bio11_REU_2100_126_MRI,
                        bio15_REU_2100_126_MRI,
                        bio16_REU_2100_126_MRI,
                        bio17_REU_2100_126_MRI
                        )
REU_clim_2040_585_MRI <- stack(
                        bio1_REU_2040_585_MRI,
                        bio4_REU_2040_585_MRI,
                        bio10_REU_2040_585_MRI,
                        bio11_REU_2040_585_MRI,
                        bio15_REU_2040_585_MRI,
                        bio16_REU_2040_585_MRI,
                        bio17_REU_2040_585_MRI
                        )
REU_clim_2100_585_MRI <- stack(
                        bio1_REU_2100_585_MRI,
                        bio4_REU_2100_585_MRI,
                        bio10_REU_2100_585_MRI,
                        bio11_REU_2100_585_MRI,
                        bio15_REU_2100_585_MRI,
                        bio16_REU_2100_585_MRI,
                        bio17_REU_2100_585_MRI
                        )

# Name layers of stacks
names(REU_clim_2040_126_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_126_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2040_585_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_585_MRI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(REU_clim_current)
plot(REU_clim_2040_126_MRI)
plot(REU_clim_2100_126_MRI)
plot(REU_clim_2040_585_MRI)
plot(REU_clim_2100_585_MRI)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
REU_SpheCoop_proj_2040_126_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_126_MRI,
                        proj.name = "REU_MRI_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_126_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_126_MRI,
                        projection.output = REU_SpheCoop_proj_2040_126_MRI,
                        proj.name = "REU_MRI_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
REU_SpheCoop_proj_2100_126_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_126_MRI,
                        proj.name = "REU_MRI_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_126_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_126_MRI,
                        projection.output = REU_SpheCoop_proj_2100_126_MRI,
                        proj.name = "REU_MRI_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
REU_SpheCoop_proj_2040_585_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_585_MRI,
                        proj.name = "REU_MRI_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_585_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_585_MRI,
                        projection.output = REU_SpheCoop_proj_2040_585_MRI,
                        proj.name = "REU_MRI_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
REU_SpheCoop_proj_2100_585_MRI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_585_MRI,
                        proj.name = "REU_MRI_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_585_MRI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_585_MRI,
                        projection.output = REU_SpheCoop_proj_2100_585_MRI,
                        proj.name = "REU_MRI_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
REU_stk_SpheCoop_ef_current_MRI <- get_predictions(REU_SpheCoop_proj_ensamble_current)
plot(REU_stk_SpheCoop_ef_current_MRI)
REU_stk_SpheCoop_ef_2040_126_MRI <- get_predictions(REU_SpheCoop_proj_ensamble_2040_126_MRI)
plot(REU_stk_SpheCoop_ef_2040_126_MRI)
REU_stk_SpheCoop_ef_2100_126_MRI <- get_predictions(REU_SpheCoop_proj_ensamble_2100_126_MRI)
plot(REU_stk_SpheCoop_ef_2100_126_MRI)
REU_stk_SpheCoop_ef_2040_585_MRI <- get_predictions(REU_SpheCoop_proj_ensamble_2040_585_MRI)
plot(REU_stk_SpheCoop_ef_2040_585_MRI)
REU_stk_SpheCoop_ef_2100_585_MRI <- get_predictions(REU_SpheCoop_proj_ensamble_2100_585_MRI)
plot(REU_stk_SpheCoop_ef_2100_585_MRI)

# Set it to only use EMca and EMwmean
subset(REU_stk_SpheCoop_ef_current_MRI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_current_MRI )))
subset(REU_stk_SpheCoop_ef_2040_126_MRI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_126_MRI )))
subset(REU_stk_SpheCoop_ef_2100_126_MRI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_126_MRI )))
subset(REU_stk_SpheCoop_ef_2040_585_MRI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_585_MRI )))
subset(REU_stk_SpheCoop_ef_2100_585_MRI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_585_MRI )))

#=========================================================
# MPI-CM6A
#=========================================================
# Garden Route
#---------------------------------------------------------
# Load future data ssp126
bio1_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio1_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio4_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio4_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio10_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio10_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio11_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio11_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio15_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio15_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio16_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio16_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio17_GR_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2041-2070\\CHELSA_bio17_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")

bio1_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio1_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio4_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio4_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio10_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio11_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio11_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio15_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio16_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio16_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")
bio17_GR_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp126\\2071-2100\\CHELSA_bio17_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_GR.tif")

# Load future data ssp585
bio1_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio1_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio4_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio4_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio10_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio10_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio11_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio11_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio15_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio15_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio16_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio16_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio17_GR_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2041-2070\\CHELSA_bio17_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")

bio1_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio1_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio4_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio4_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio10_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio11_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio11_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio15_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio16_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio16_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")
bio17_GR_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Garden route\\ssp585\\2071-2100\\CHELSA_bio17_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_GR.tif")

# Create raster stacks
GR_clim_2040_126_MPI <- stack(
                        bio1_GR_2040_126_MPI,
                        bio4_GR_2040_126_MPI,
                        bio10_GR_2040_126_MPI,
                        bio11_GR_2040_126_MPI,
                        bio15_GR_2040_126_MPI,
                        bio16_GR_2040_126_MPI,
                        bio17_GR_2040_126_MPI
                        )
GR_clim_2100_126_MPI <- stack(
                        bio1_GR_2100_126_MPI,
                        bio4_GR_2100_126_MPI,
                        bio10_GR_2100_126_MPI,
                        bio11_GR_2100_126_MPI,
                        bio15_GR_2100_126_MPI,
                        bio16_GR_2100_126_MPI,
                        bio17_GR_2100_126_MPI
                        )
GR_clim_2040_585_MPI <- stack(
                        bio1_GR_2040_585_MPI,
                        bio4_GR_2040_585_MPI,
                        bio10_GR_2040_585_MPI,
                        bio11_GR_2040_585_MPI,
                        bio15_GR_2040_585_MPI,
                        bio16_GR_2040_585_MPI,
                        bio17_GR_2040_585_MPI
                        )
GR_clim_2100_585_MPI <- stack(
                        bio1_GR_2100_585_MPI,
                        bio4_GR_2100_585_MPI,
                        bio10_GR_2100_585_MPI,
                        bio11_GR_2100_585_MPI,
                        bio15_GR_2100_585_MPI,
                        bio16_GR_2100_585_MPI,
                        bio17_GR_2100_585_MPI
                        )

# Name layers of stacks
names(GR_clim_2040_126_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(GR_clim_2040_126_MPI)
plot(GR_clim_2100_126_MPI)
plot(GR_clim_2040_585_MPI)
plot(GR_clim_2100_585_MPI)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
GR_SpheCoop_proj_2040_126_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126_MPI,
                        proj.name = "GR_MPI_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_126_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126_MPI,
                        projection.output = GR_SpheCoop_proj_2040_126_MPI,
                        proj.name = "GR_MPI_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
GR_SpheCoop_proj_2100_126_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126_MPI,
                        proj.name = "GR_MPI_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_126_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126_MPI,
                        projection.output = GR_SpheCoop_proj_2040_126_MPI,
                        proj.name = "GR_MPI_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
GR_SpheCoop_proj_2040_585_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585_MPI,
                        proj.name = "GR_MPI_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_585_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585_MPI,
                        projection.output = GR_SpheCoop_proj_2040_585_MPI,
                        proj.name = "GR_MPI_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
GR_SpheCoop_proj_2100_585_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585_MPI,
                        proj.name = "GR_MPI_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_585_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585_MPI,
                        projection.output = GR_SpheCoop_proj_2040_585_MPI,
                        proj.name = "GR_MPI_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
GR_stk_SpheCoop_ef_current_MPI <- get_predictions(GR_SpheCoop_proj_ensamble_current)
plot(GR_stk_SpheCoop_ef_current_MPI)
GR_stk_SpheCoop_ef_2040_126_MPI <- get_predictions(GR_SpheCoop_proj_ensamble_2040_126_MPI)
plot(GR_stk_SpheCoop_ef_2040_126_MPI)
GR_stk_SpheCoop_ef_2100_126_MPI <- get_predictions(GR_SpheCoop_proj_ensamble_2100_126_MPI)
plot(GR_stk_SpheCoop_ef_2100_126_MPI)
GR_stk_SpheCoop_ef_2040_585_MPI <- get_predictions(GR_SpheCoop_proj_ensamble_2040_585_MPI)
plot(GR_stk_SpheCoop_ef_2040_585_MPI)
GR_stk_SpheCoop_ef_2100_585_MPI <- get_predictions(GR_SpheCoop_proj_ensamble_2100_585_MPI)
plot(GR_stk_SpheCoop_ef_2100_585_MPI)

# Set it to only use EMca and EMwmean
subset(GR_stk_SpheCoop_ef_current_MPI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_current_MPI )))
subset(GR_stk_SpheCoop_ef_2040_126_MPI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_126_MPI )))
subset(GR_stk_SpheCoop_ef_2100_126_MPI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_126_MPI )))
subset(GR_stk_SpheCoop_ef_2040_585_MPI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_585_MPI )))
subset(GR_stk_SpheCoop_ef_2100_585_MPI ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_585_MPI )))


#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
# Loading raster data
# Load future data ssp126
bio1_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio1_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio4_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio4_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio10_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio10_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio11_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio11_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio15_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio15_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio16_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio16_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio17_REU_2040_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2041-2070\\CHELSA_bio17_2041-2070_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")

bio1_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio1_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio4_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio4_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio10_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio11_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio11_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio15_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio16_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio16_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")
bio17_REU_2100_126_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp126\\2071-2100\\CHELSA_bio17_2071-2100_mpi-esm1-2-hr_ssp126_V.2.1_R.tif")

# Load future data ssp585
bio1_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio1_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio4_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio4_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio10_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio10_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio11_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio11_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio15_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio15_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio16_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio16_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio17_REU_2040_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2041-2070\\CHELSA_bio17_2041-2070_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")

bio1_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio1_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio4_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio4_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio10_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio10_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio11_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio11_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio15_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio15_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio16_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio16_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")
bio17_REU_2100_585_MPI <- raster("input_data\\climate_data\\MPI-CM6A\\Reunion\\ssp585\\2071-2100\\CHELSA_bio17_2071-2100_mpi-esm1-2-hr_ssp585_V.2.1_R.tif")

# Create stacks of climatic data
REU_clim_2040_126_MPI <- stack(
                        bio1_REU_2040_126_MPI,
                        bio4_REU_2040_126_MPI,
                        bio10_REU_2040_126_MPI,
                        bio11_REU_2040_126_MPI,
                        bio15_REU_2040_126_MPI,
                        bio16_REU_2040_126_MPI,
                        bio17_REU_2040_126_MPI
                        )
REU_clim_2100_126_MPI <- stack(
                        bio1_REU_2100_126_MPI,
                        bio4_REU_2100_126_MPI,
                        bio10_REU_2100_126_MPI,
                        bio11_REU_2100_126_MPI,
                        bio15_REU_2100_126_MPI,
                        bio16_REU_2100_126_MPI,
                        bio17_REU_2100_126_MPI
                        )
REU_clim_2040_585_MPI <- stack(
                        bio1_REU_2040_585_MPI,
                        bio4_REU_2040_585_MPI,
                        bio10_REU_2040_585_MPI,
                        bio11_REU_2040_585_MPI,
                        bio15_REU_2040_585_MPI,
                        bio16_REU_2040_585_MPI,
                        bio17_REU_2040_585_MPI
                        )
REU_clim_2100_585_MPI <- stack(
                        bio1_REU_2100_585_MPI,
                        bio4_REU_2100_585_MPI,
                        bio10_REU_2100_585_MPI,
                        bio11_REU_2100_585_MPI,
                        bio15_REU_2100_585_MPI,
                        bio16_REU_2100_585_MPI,
                        bio17_REU_2100_585_MPI
                        )

# Name layers of stacks
names(REU_clim_2040_126_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_126_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2040_585_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_585_MPI) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(REU_clim_current)
plot(REU_clim_2040_126_MPI)
plot(REU_clim_2100_126_MPI)
plot(REU_clim_2040_585_MPI)
plot(REU_clim_2100_585_MPI)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
REU_SpheCoop_proj_2040_126_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_126_MPI,
                        proj.name = "REU_MPI_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_126_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_126_MPI,
                        projection.output = REU_SpheCoop_proj_2040_126_MPI,
                        proj.name = "REU_MPI_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
REU_SpheCoop_proj_2100_126_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_126_MPI,
                        proj.name = "REU_MPI_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_126_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_126_MPI,
                        projection.output = REU_SpheCoop_proj_2100_126_MPI,
                        proj.name = "REU_MPI_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
REU_SpheCoop_proj_2040_585_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_585_MPI,
                        proj.name = "REU_MPI_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_585_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_585_MPI,
                        projection.output = REU_SpheCoop_proj_2040_585_MPI,
                        proj.name = "REU_MPI_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
REU_SpheCoop_proj_2100_585_MPI <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_585_MPI,
                        proj.name = "REU_MPI_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_585_MPI <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_585_MPI,
                        projection.output = REU_SpheCoop_proj_2100_585_MPI,
                        proj.name = "REU_MPI_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.clamping.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
REU_stk_SpheCoop_ef_current_MPI <- get_predictions(REU_SpheCoop_proj_ensamble_current)
plot(REU_stk_SpheCoop_ef_current_MPI)
REU_stk_SpheCoop_ef_2040_126_MPI <- get_predictions(REU_SpheCoop_proj_ensamble_2040_126_MPI)
plot(REU_stk_SpheCoop_ef_2040_126_MPI)
REU_stk_SpheCoop_ef_2100_126_MPI <- get_predictions(REU_SpheCoop_proj_ensamble_2100_126_MPI)
plot(REU_stk_SpheCoop_ef_2100_126_MPI)
REU_stk_SpheCoop_ef_2040_585_MPI <- get_predictions(REU_SpheCoop_proj_ensamble_2040_585_MPI)
plot(REU_stk_SpheCoop_ef_2040_585_MPI)
REU_stk_SpheCoop_ef_2100_585_MPI <- get_predictions(REU_SpheCoop_proj_ensamble_2100_585_MPI)
plot(REU_stk_SpheCoop_ef_2100_585_MPI)

# Set it to only use EMca and EMwmean
subset(REU_stk_SpheCoop_ef_current_MPI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_current_MPI )))
subset(REU_stk_SpheCoop_ef_2040_126_MPI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_126_MPI )))
subset(REU_stk_SpheCoop_ef_2100_126_MPI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_126_MPI )))
subset(REU_stk_SpheCoop_ef_2040_585_MPI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_585_MPI )))
subset(REU_stk_SpheCoop_ef_2100_585_MPI ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_585_MPI )))

#=========================================================
# UKESM
#=========================================================
# Garden Route
#---------------------------------------------------------
# Load future data ssp126
bio1_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio1_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio4_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio4_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio10_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio10_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio11_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio11_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio15_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio15_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio16_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio16_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio17_GR_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2041-2070\\CHELSA_bio17_2041-2070_ukesm1-0-ll_ssp126_V.2.1_GR.tif")

bio1_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio1_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio4_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio4_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio10_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio11_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio11_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio15_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio16_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio16_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")
bio17_GR_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp126\\2071-2100\\CHELSA_bio17_2071-2100_ukesm1-0-ll_ssp126_V.2.1_GR.tif")

# Load future data ssp585
bio1_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio1_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio4_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio4_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio10_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio10_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio11_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio11_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio15_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio15_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio16_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio16_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio17_GR_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2041-2070\\CHELSA_bio17_2041-2070_ukesm1-0-ll_ssp585_V.2.1_GR.tif")

bio1_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio1_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio4_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio4_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio10_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio11_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio11_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio15_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio16_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio16_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")
bio17_GR_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Garden route\\ssp585\\2071-2100\\CHELSA_bio17_2071-2100_ukesm1-0-ll_ssp585_V.2.1_GR.tif")

# Create raster stacks
GR_clim_2040_126_UKESM <- stack(
                        bio1_GR_2040_126_UKESM,
                        bio4_GR_2040_126_UKESM,
                        bio10_GR_2040_126_UKESM,
                        bio11_GR_2040_126_UKESM,
                        bio15_GR_2040_126_UKESM,
                        bio16_GR_2040_126_UKESM,
                        bio17_GR_2040_126_UKESM
                        )
GR_clim_2100_126_UKESM <- stack(
                        bio1_GR_2100_126_UKESM,
                        bio4_GR_2100_126_UKESM,
                        bio10_GR_2100_126_UKESM,
                        bio11_GR_2100_126_UKESM,
                        bio15_GR_2100_126_UKESM,
                        bio16_GR_2100_126_UKESM,
                        bio17_GR_2100_126_UKESM
                        )
GR_clim_2040_585_UKESM <- stack(
                        bio1_GR_2040_585_UKESM,
                        bio4_GR_2040_585_UKESM,
                        bio10_GR_2040_585_UKESM,
                        bio11_GR_2040_585_UKESM,
                        bio15_GR_2040_585_UKESM,
                        bio16_GR_2040_585_UKESM,
                        bio17_GR_2040_585_UKESM
                        )
GR_clim_2100_585_UKESM <- stack(
                        bio1_GR_2100_585_UKESM,
                        bio4_GR_2100_585_UKESM,
                        bio10_GR_2100_585_UKESM,
                        bio11_GR_2100_585_UKESM,
                        bio15_GR_2100_585_UKESM,
                        bio16_GR_2100_585_UKESM,
                        bio17_GR_2100_585_UKESM
                        )

# Name layers of stacks
names(GR_clim_2040_126_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_126_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2040_585_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(GR_clim_2100_585_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )

plot(GR_clim_2040_126_UKESM)
plot(GR_clim_2100_126_UKESM)
plot(GR_clim_2040_585_UKESM)
plot(GR_clim_2100_585_UKESM)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
GR_SpheCoop_proj_2040_126_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_126_UKESM,
                        proj.name = "GR_UKESM_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_126_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_126_UKESM,
                        projection.output = GR_SpheCoop_proj_2040_126_UKESM,
                        proj.name = "GR_UKESM_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
GR_SpheCoop_proj_2100_126_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_126_UKESM,
                        proj.name = "GR_UKESM_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_126_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_126_UKESM,
                        projection.output = GR_SpheCoop_proj_2040_126_UKESM,
                        proj.name = "GR_UKESM_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
GR_SpheCoop_proj_2040_585_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2040_585_UKESM,
                        proj.name = "GR_UKESM_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2040_585_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2040_585_UKESM,
                        projection.output = GR_SpheCoop_proj_2040_585_UKESM,
                        proj.name = "GR_UKESM_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
GR_SpheCoop_proj_2100_585_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = GR_clim_2100_585_UKESM,
                        proj.name = "GR_UKESM_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

GR_SpheCoop_proj_ensamble_2100_585_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = GR_clim_2100_585_UKESM,
                        projection.output = GR_SpheCoop_proj_2040_585_UKESM,
                        proj.name = "GR_UKESM_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
GR_stk_SpheCoop_ef_current_UKESM <- get_predictions(GR_SpheCoop_proj_ensamble_current)
plot(GR_stk_SpheCoop_ef_current_UKESM)
GR_stk_SpheCoop_ef_2040_126_UKESM <- get_predictions(GR_SpheCoop_proj_ensamble_2040_126_UKESM)
plot(GR_stk_SpheCoop_ef_2040_126_UKESM)
GR_stk_SpheCoop_ef_2100_126_UKESM <- get_predictions(GR_SpheCoop_proj_ensamble_2100_126_UKESM)
plot(GR_stk_SpheCoop_ef_2100_126_UKESM)
GR_stk_SpheCoop_ef_2040_585_UKESM <- get_predictions(GR_SpheCoop_proj_ensamble_2040_585_UKESM)
plot(GR_stk_SpheCoop_ef_2040_585_UKESM)
GR_stk_SpheCoop_ef_2100_585_UKESM <- get_predictions(GR_SpheCoop_proj_ensamble_2100_585_UKESM)
plot(GR_stk_SpheCoop_ef_2100_585_UKESM)

# Set it to only use EMca and EMwmean
subset(GR_stk_SpheCoop_ef_current_UKESM ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_current_UKESM )))
subset(GR_stk_SpheCoop_ef_2040_126_UKESM ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_126_UKESM )))
subset(GR_stk_SpheCoop_ef_2100_126_UKESM ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_126_UKESM )))
subset(GR_stk_SpheCoop_ef_2040_585_UKESM ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2040_585_UKESM )))
subset(GR_stk_SpheCoop_ef_2100_585_UKESM ,grep('EMca|EMwmean', names(GR_stk_SpheCoop_ef_2100_585_UKESM )))


#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
# Loading raster data
# Load future data ssp126
bio1_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio1_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio4_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio4_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio10_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio10_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio11_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio11_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio15_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio15_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio16_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio16_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio17_REU_2040_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2041-2070\\CHELSA_bio17_2041-2070_ukesm1-0-ll_ssp126_V.2.1_R.tif")

bio1_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio1_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio4_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio4_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio10_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio11_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio11_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio15_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio16_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio16_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")
bio17_REU_2100_126_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp126\\2071-2100\\CHELSA_bio17_2071-2100_ukesm1-0-ll_ssp126_V.2.1_R.tif")

# Load future data ssp585
bio1_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio1_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio4_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio4_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio10_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio10_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio11_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio11_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio15_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio15_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio16_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio16_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio17_REU_2040_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2041-2070\\CHELSA_bio17_2041-2070_ukesm1-0-ll_ssp585_V.2.1_R.tif")

bio1_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio1_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio4_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio4_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio10_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio10_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio11_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio11_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio15_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio15_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio16_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio16_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")
bio17_REU_2100_585_UKESM <- raster("input_data\\climate_data\\UKESM\\Reunion\\ssp585\\2071-2100\\CHELSA_bio17_2071-2100_ukesm1-0-ll_ssp585_V.2.1_R.tif")

# Create stacks of climatic data
REU_clim_2040_126_UKESM <- stack(
                        bio1_REU_2040_126_UKESM,
                        bio4_REU_2040_126_UKESM,
                        bio10_REU_2040_126_UKESM,
                        bio11_REU_2040_126_UKESM,
                        bio15_REU_2040_126_UKESM,
                        bio16_REU_2040_126_UKESM,
                        bio17_REU_2040_126_UKESM
                        )
REU_clim_2100_126_UKESM <- stack(
                        bio1_REU_2100_126_UKESM,
                        bio4_REU_2100_126_UKESM,
                        bio10_REU_2100_126_UKESM,
                        bio11_REU_2100_126_UKESM,
                        bio15_REU_2100_126_UKESM,
                        bio16_REU_2100_126_UKESM,
                        bio17_REU_2100_126_UKESM
                        )
REU_clim_2040_585_UKESM <- stack(
                        bio1_REU_2040_585_UKESM,
                        bio4_REU_2040_585_UKESM,
                        bio10_REU_2040_585_UKESM,
                        bio11_REU_2040_585_UKESM,
                        bio15_REU_2040_585_UKESM,
                        bio16_REU_2040_585_UKESM,
                        bio17_REU_2040_585_UKESM
                        )
REU_clim_2100_585_UKESM <- stack(
                        bio1_REU_2100_585_UKESM,
                        bio4_REU_2100_585_UKESM,
                        bio10_REU_2100_585_UKESM,
                        bio11_REU_2100_585_UKESM,
                        bio15_REU_2100_585_UKESM,
                        bio16_REU_2100_585_UKESM,
                        bio17_REU_2100_585_UKESM
                        )

# Name layers of stacks
names(REU_clim_2040_126_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_126_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2040_585_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
names(REU_clim_2100_585_UKESM) <- c(
                        "AnnMeanTemp",
                        "TempSeason",
                        "MeanTempWarQ",
                        "MeanTempColdQ",
                        "PrecipSeason",
                        "PrecipWetQ",
                        "PrecipDryQ"
                        )
plot(REU_clim_current)
plot(REU_clim_2040_126_UKESM)
plot(REU_clim_2100_126_UKESM)
plot(REU_clim_2040_585_UKESM)
plot(REU_clim_2100_585_UKESM)

#---------------------------------------------------------
# Future projections
# Do projection of 2040 ssp126
REU_SpheCoop_proj_2040_126_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_126_UKESM,
                        proj.name = "REU_UKESM_proj_2040_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_126_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_126_UKESM,
                        projection.output = REU_SpheCoop_proj_2040_126_UKESM,
                        proj.name = "REU_UKESM_proj_ens_2040_126",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp126
REU_SpheCoop_proj_2100_126_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_126_UKESM,
                        proj.name = "REU_UKESM_proj_2100_126",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_126_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_126_UKESM,
                        projection.output = REU_SpheCoop_proj_2100_126_UKESM,
                        proj.name = "REU_UKESM_proj_ens_2100_126",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2040 ssp585
REU_SpheCoop_proj_2040_585_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2040_585_UKESM,
                        proj.name = "REU_UKESM_proj_2040_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2040_585_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2040_585_UKESM,
                        projection.output = REU_SpheCoop_proj_2040_585_UKESM,
                        proj.name = "REU_UKESM_proj_ens_2040_585",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

# Do projection of 2100 ssp585
REU_SpheCoop_proj_2100_585_UKESM <- BIOMOD_Projection (
                        bm.mod = ch_SpheCoop_models,
                        new.env = REU_clim_2100_585_UKESM,
                        proj.name = "REU_UKESM_proj_2100_585",
                        metric.binary="TSS",
                        do.stack = FALSE
)

REU_SpheCoop_proj_ensamble_2100_585_UKESM <- BIOMOD_EnsembleForecasting(
                        bm.em = ch_SpheCoop_ensemble,
                        new.env = REU_clim_2100_585_UKESM,
                        projection.output = REU_SpheCoop_proj_2100_585_UKESM,
                        proj.name = "REU_UKESM_proj_ens_2100_585",
                        metric.binary="TSS",
                        build.claUKESMng.mask=TRUE,
                        do.stack=FALSE
)

#---------------------------------------------------------
# Get the ensemble models projection stack
REU_stk_SpheCoop_ef_current_UKESM <- get_predictions(REU_SpheCoop_proj_ensamble_current)
plot(REU_stk_SpheCoop_ef_current_UKESM)
REU_stk_SpheCoop_ef_2040_126_UKESM <- get_predictions(REU_SpheCoop_proj_ensamble_2040_126_UKESM)
plot(REU_stk_SpheCoop_ef_2040_126_UKESM)
REU_stk_SpheCoop_ef_2100_126_UKESM <- get_predictions(REU_SpheCoop_proj_ensamble_2100_126_UKESM)
plot(REU_stk_SpheCoop_ef_2100_126_UKESM)
REU_stk_SpheCoop_ef_2040_585_UKESM <- get_predictions(REU_SpheCoop_proj_ensamble_2040_585_UKESM)
plot(REU_stk_SpheCoop_ef_2040_585_UKESM)
REU_stk_SpheCoop_ef_2100_585_UKESM <- get_predictions(REU_SpheCoop_proj_ensamble_2100_585_UKESM)
plot(REU_stk_SpheCoop_ef_2100_585_UKESM)

# Set it to only use EMca and EMwmean
subset(REU_stk_SpheCoop_ef_current_UKESM ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_current_UKESM )))
subset(REU_stk_SpheCoop_ef_2040_126_UKESM ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_126_UKESM )))
subset(REU_stk_SpheCoop_ef_2100_126_UKESM ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_126_UKESM )))
subset(REU_stk_SpheCoop_ef_2040_585_UKESM ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2040_585_UKESM )))
subset(REU_stk_SpheCoop_ef_2100_585_UKESM ,grep('EMca|EMwmean', names(REU_stk_SpheCoop_ef_2100_585_UKESM )))
