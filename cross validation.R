#---------------------------------------------------------
### Load packages ###
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

# --------------------------------------------------------------
# 1. Formatting data
# Format Data with pseudo-absences 
myBiomodData <- BIOMOD_FormatingData(resp.var = rep(1, nrow(data)), 
                                     expl.var = ch_current_clim, 
                                     resp.xy = data[,c('decimalLongitude','decimalLatitude')], 
                                     resp.name = "Sphaeropteris.cooperi",
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 1281,
                                     PA.strategy = 'random',
                                     PA.dist.min = 0, 
                                     PA.dist.max = NULL,)
class(myBiomodData)
typeof(myBiomodData)
str(myBiomodData)
slotNames(myBiomodData)
plot(myBiomodData)
# --------------------------------------------------------------
# 2. Defining Models Options 
# 2.1 Caret-tuned hyper-parameters
BIOMOD_Tuning(
  myBiomodData,
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

# 2.2 Specify tuned parameters
myBiomodOption <- BIOMOD_ModelingOptions(
                        GLM = list( type = 'polynomial',
                                    interaction.level = 1,
                                    myFormula = Sphaeropteris.cooperi ~
                                                AnnMeanTemp +
                                                I(AnnMeanTemp^2) +
                                                I(AnnMeanTemp^3) +
                                                TempSeason +
                                                I(TempSeason^2) +
                                                I(TempSeason^3) +
                                                MeanTempWarQ +
                                                I(MeanTempWarQ^2) +
                                                MeanTempColdQ +
                                                I(MeanTempColdQ^3) +
                                                PrecipSeason +
                                                I(PrecipSeason^2) +
                                                I(PrecipSeason^3) +
                                                PrecipWetQ +
                                                I(PrecipWetQ^2) +
                                                I(PrecipWetQ^3) +
                                                PrecipDryQ +
                                                I(PrecipDryQ^2) +
                                                AnnMeanTemp:TempSeason +
                                                AnnMeanTemp:MeanTempWarQ +
                                                AnnMeanTemp:MeanTempColdQ +
                                                AnnMeanTemp:PrecipSeason +
                                                AnnMeanTemp:PrecipDryQ +
                                                TempSeason:MeanTempWarQ + 
                                                TempSeason:MeanTempColdQ + 
                                                TempSeason:PrecipSeason + 
                                                TempSeason:PrecipWetQ + 
                                                TempSeason:PrecipDryQ + 
                                                MeanTempWarQ:MeanTempColdQ + 
                                                MeanTempWarQ:PrecipSeason + 
                                                MeanTempWarQ:PrecipWetQ + 
                                                MeanTempWarQ:PrecipDryQ + 
                                                MeanTempColdQ:PrecipSeason + 
                                                MeanTempColdQ:PrecipWetQ + 
                                                MeanTempColdQ:PrecipDryQ + 
                                                PrecipSeason:PrecipWetQ + 
                                                PrecipSeason:PrecipDryQ + 
                                                PrecipWetQ:PrecipDryQ,
                                    test = 'none',
                                    family = binomial(link = 'logit'),
                                    mustart = 0.5,
                                    control = glm.control(epsilon = 1e-08,
                                                          maxit = 50,
                                                          trace = FALSE)),
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
                                    family = binomial(link = "logit"),
                                    method = 'GCV.Cp',
                                    optimizer = c("outer", "newton"),
                                    select = FALSE,
                                    knots = NULL,
                                    paraPen = NULL,
                                    control = list(nthreads = 1,
                                            ncv.threads = 1,
                                            irls.reg = 0,
                                            epsilon = 1e-07,
                                            maxit = 200,
                                            trace = FALSE,
                                            mgcv.tol = 1e-07,
                                            mgcv.half = 15,
                                            rank.tol = 1.49011611938477e-08,
                                            nlm = list(ndigit = 7,
                                                       gradtol = 1e-06,
                                                       stepmax = 2,
                                                       steptol = 1e-04,
                                                       iterlim = 200,
                                                       check.analyticals = 0),
                                            optim = list(factr = 1e+07),
                                            newton = list(conv.tol = 1e-06,
                                                          maxNstep = 5,
                                                          maxSstep = 2,
                                                          maxHalf = 30,
                                                          use.svd = 0),
                                            idLinksBases = TRUE,
                                            scalePenalty = TRUE,
                                            efs.lspmax = 15,
                                            efs.tol = 0.1,
                                            keepData = FALSE,
                                            edge.correct = FALSE) 
                                    ),
                        FDA = list(method = 'mars',
                                    add_args = list(degree = 2, 
                                                    nprune = 19)),
                        MARS = list(type = 'simple',
                                    interaction.level = 0,
                                    myFormula = NULL,
                                    nk = NULL,
                                    penalty = 2,
                                    thresh = 0.001,
                                    nprune = 15,
                                    pmethod = "backward"),
                        MAXENT = list(path_to_maxent.jar = "C:\\Users\\chris\\Desktop\\MegMSc\\maxent",
                                        memory_allocated = 512,
                                        initial_heap_size = NULL,
                                        max_heap_size = NULL,
                                        background_data_dir = "default",
                                        maximumbackground = "default",
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
                        RF = list(do.classif = TRUE,
                                    ntree = 500,
                                    mtry = 3,
                                    sampsize = NULL,
                                    nodesize = 5,
                                    maxnodes = NULL),
                        CTA = list(method = "class",
                                    parms = "default",
                                    cost = NULL,
                                    control = list(xval = 5, 
                                                  minbucket = 5, 
                                                  minsplit = 5, 
                                                  cp = list(cp=0), 
                                                  maxdepth = list(maxdepth=8))),
                        ANN = list(NbCV = 5,
                                    size = 6,
                                    decay = 0.001,
                                    rang = 0.1,
                                    maxit = 500),
                        SRE = list(quantile(0.25))
                    )


# --------------------------------------------------------------
# 3. Cross-validation
# 3.1 block selection
cv.b <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = "block")
head(cv.b)
apply(cv.b, 2, table)

# 3.2 k-fold selection
cv.k <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = "kfold",
                           nb.rep = 2,
                           k = 3)
head(cv.k)
apply(cv.k, 2, table)

# 3.3 random selection
cv.r <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = "random",
                           nb.rep = 3,
                           k = 0.8)
head(cv.r)
apply(cv.r, 2, table)

# 3.4 stratified selection (geographic)
cv.s <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = "strat",
                           k = 2,
                           balance = "presences",
                           strat = "x")
head(cv.s)
apply(cv.s, 2, table)

# 3.5 stratified selection (environmental)
cv.e <- bm_CrossValidation(bm.format = myBiomodData,
                           strategy = "env",
                           k = 2,
                           balance = "presences")
head(cv.e)
apply(cv.e, 2, table)


# 4. Modelling
# allModels <- c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXENT", "RF")
allModels<- c("GLM")
# 4.1 Model single models using block validation
myBiomodModelOutCV.S <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = "cv.s",
                                    models = allModels,
                                    CV.strategy = "strat",
                                    CV.perc = 0.80,
                                    CV.strat = 'x',
                                    CV.k = 2,
                                    CV.balance = "presences",
                                    bm.options = myBiomodOption,
                                    nb.rep = 3,
                                    data.split.table = myBiomodCV,
                                    metric.eval = c("TSS", "ROC"),
                                    var.import = 4,
                                    do.full.models = FALSE,
                                    seed.val = 42
                                    )
 
#To get/verify the cross-validation table and visualize it on your BIOMOD.formated.data or BIOMOD.formated.data.PA object.
#TRUE for calibration and FALSE for validation
myCalibLinesCV.S <- get_calib_lines(myBiomodModelOutCV.S)
class(myBiomodModelOutCV.E)
class(myCalibLinesCV.E)
typeof(myBiomodModelOutCV.E)
typeof(myCalibLinesCV.E)
str(myBiomodModelOutCV.E)
str(myCalibLinesCV.E)
plot(myBiomodModelOutCV.E)
# List available slots
slotNames(myBiomodModelOutCV.E)

# Access specific slots
myBiomodModelOutCV.E@models.evaluation
myBiomodModelOutCV.E@models.prediction
myBiomodModelOutCV.E@variables.importance
# And so on...

# Once you identify the data you want to plot, you can use standard R plotting functions to create your plots.

plot(myBiomodModelOutCV.B, calib.lines = myCalibLinesCV.B)

# Extract evaluation data from BIOMOD.models.out
evaluation_data <- get_evaluations(myBiomodModelOutCV.E)

# Extract calibration lines
calibration_lines <- myCalibLinesCV.E

# Example plot (replace with your actual plotting code)
# For example, let's plot sensitivity against specificity
plot(evaluation_data$sensitivity, evaluation_data$specificity, 
     xlab = "Sensitivity", ylab = "Specificity",
     main = "Sensitivity vs Specificity")

# Add calibration lines to the plot
abline(v = which(calibration_lines), col = "red", lty = 2)

# Get evaluation scores & variables importance
myEvalCV.B <- get_evaluations(myBiomodModelOutCV.B)
myEvalCV.B$CV.strategy <- "block" 
myEvalCV.B$CV.strategy[grepl("13", myEvalCV.B$full.name)] <- "Full"
myEvalCV.B$CV.strategy[grepl("11|12", myEvalCV.B$full.name)] <- "Block"
head(myEvalCV.B)

boxplot(myEvalCV.B$calibration ~ interaction(myEvalCV.B$algo, myEvalCV.B$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))

boxplot(myEvalCV.B$validation ~ interaction(myEvalCV.B$algo, myEvalCV.B$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))

# 4.2 Model single models using k-fold validation
myBiomodModelOutCV.K <- BIOMOD_Modeling(bm.format = myBiomodData,
                                        modeling.id = "cv.k",
                                        models = c("RF"),
                                        CV.strategy = "kfold",
                                        k = 2,
                                        bm.options = myBiomodOption,
                                        metric.eval = c("TSS", "ROC"),
                                        var.import = 1,
                                        seed.val = 42
                                    )
 
#To get/verify the cross-validation table and visualize it on your BIOMOD.formated.data or BIOMOD.formated.data.PA object.
#TRUE for calibration and FALSE for validation
myCalibLinesCV.K <- get_calib_lines(myBiomodModelOutCV.K)
plot(myBiomodModelOutCV.K, calib.lines = myCalibLinesCV.K)

# Get evaluation scores & variables importance
myEvalCV.K <- get_evaluations(myBiomodModelOutCV.K) 
myEvalCV.K$CV.strategy <- "k-fold" 
myEvalCV.K$CV.strategy[grepl("13", myEvalCV.K$full.name)] <- "Full"
myEvalCV.K$CV.strategy[grepl("11|12", myEvalCV.K$full.name)] <- "k-fold"
head(myEvalCV.K)

boxplot(myEvalCV.K$calibration ~ interaction(myEvalCV.K$algo, myEvalCV.K$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))

boxplot(myEvalCV.K$validation ~ interaction(myEvalCV.K$algo, myEvalCV.K$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))


# 4.3 Model single models using random selection validation
myBiomodModelOutCV.R <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = "cv.r",
                                    models = allModels,
                                    CV.strategy = "random",
                                    bm.options = myBiomodOption,
                                    nb.rep = 3,
                                    metric.eval = c("TSS", "ROC"),
                                    var.import = 4,
                                    do.full.models = FALSE,
                                    seed.val = 42
                                    )
 
#To get/verify the cross-validation table and visualize it on your BIOMOD.formated.data or BIOMOD.formated.data.PA object.
#TRUE for calibration and FALSE for validation
myCalibLinesCV.R <- get_calib_lines(myBiomodModelOutCV.R)
plot(myBiomodModelOutCV.R, calib.lines = myCalibLinesCV.R)

# Get evaluation scores & variables importance
myEvalCV.R <- get_evaluations(myBiomodModelOutCV.R) 
myEvalCV.R$CV.strategy <- "random" 
myEvalCV.R$CV.strategy[grepl("13", myEvalCV.R$full.name)] <- "Full"
myEvalCV.R$CV.strategy[grepl("11|12", myEvalCV.R$full.name)] <- "random"
head(myEvalCV.R)

boxplot(myEvalCV.R$calibration ~ interaction(myEvalCV.R$algo, myEvalCV.R$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))

boxplot(myEvalCV.R$validation ~ interaction(myEvalCV.R$algo, myEvalCV.R$CV.strategy),
        xlab = "",
        ylab = "ROC AUC",
        col = rep(c("brown", "cadetblue"), 3))

