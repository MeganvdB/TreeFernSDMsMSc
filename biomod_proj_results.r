library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(terra)
library(mgcv)
#install.packages("rasterVis")
library(rasterVis)

#=========================================================
# GFDL-ESM4
#=========================================================
# Garden Route
#---------------------------------------------------------
#---------------------------------------------------------
# Species range change setup
# Garden route
Sphe_coop_bin_proj_GR_current <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_126_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_126_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_585_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_585_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_GFDL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

#---------------------------------------------------------
# Garden route ssp126
GR_SRC_current_2040_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_126_GFDL
)
GR_SRC_current_2040_126_GFDL$Compt.By.Models

GR_SRC_current_2100_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_126_GFDL
)
GR_SRC_current_2100_126_GFDL$Compt.By.Models

GR_SRC_2040_2100_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_126_GFDL,
    Sphe_coop_bin_proj_GR_2100_126_GFDL
)
GR_SRC_2040_2100_126_GFDL$Compt.By.Models

# Garden route ssp585
GR_SRC_current_2040_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_585_GFDL
)
GR_SRC_current_2040_585_GFDL$Compt.By.Models

GR_SRC_current_2100_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_585_GFDL
)
GR_SRC_current_2100_585_GFDL$Compt.By.Models

GR_SRC_2040_2100_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_585_GFDL,
    Sphe_coop_bin_proj_GR_2100_585_GFDL
)
GR_SRC_2040_2100_585_GFDL$Compt.By.Models

#---------------------------------------------------------
# Plot species range change

# Load libraries to change levelplot theme
#install.packages("viridis")
library(viridis)

# Garden route ssp126
GR_raster_2040_126_GFDL <- GR_SRC_current_2040_126_GFDL$Diff.By.Pixel
GR_raster_2100_126_GFDL <- GR_SRC_current_2100_126_GFDL$Diff.By.Pixel
GR_raster_future_126_GFDL <- GR_SRC_2040_2100_126_GFDL$Diff.By.Pixel

plot(GR_raster_2040_126_GFDL)
plot(GR_raster_2100_126_GFDL)
plot(GR_raster_future_126_GFDL)
GR_layers_126_GFDL<-c(GR_raster_2040_126_GFDL,GR_raster_2100_126_GFDL, GR_raster_future_126_GFDL)
summary(GR_layers_126_GFDL)
GR_SpheCoop_src_map_126_GFDL <- stack(GR_layers_126_GFDL)
summary(GR_SpheCoop_src_map_126_GFDL)
names(GR_SpheCoop_src_map_126_GFDL) <- c("GR_GFDL_current_2040_126_ca", "GR_GFDL_current_2040_126_wm","GR_GFDL_current_2100_126_ca","GR_GFDL_current_2100_126_wm","GR_GFDL_2040_2100_126_ca","GR_GFDL_2040_2100_126_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(GR_SpheCoop_src_map_126_GFDL, main = 'SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


# Garden route 585
GR_raster_2040_585_GFDL <- GR_SRC_current_2040_585_GFDL$Diff.By.Pixel
GR_raster_2100_585_GFDL <- GR_SRC_current_2100_585_GFDL$Diff.By.Pixel
GR_raster_future_585_GFDL <- GR_SRC_2040_2100_585_GFDL$Diff.By.Pixel

plot(GR_raster_2040_585_GFDL)
plot(GR_raster_2100_585_GFDL)
plot(GR_raster_future_585_GFDL)
GR_layers_585_GFDL<-c(GR_raster_2040_585_GFDL,GR_raster_2100_585_GFDL, GR_raster_future_585_GFDL)
summary(GR_layers_585_GFDL)
GR_SpheCoop_src_map_585_GFDL <- stack(GR_layers_585_GFDL)
summary(GR_SpheCoop_src_map_585_GFDL)
names(GR_SpheCoop_src_map_585_GFDL) <- c("GR_GFDL_current_2040_585_ca", "GR_GFDL_current_2040_585_wm","GR_GFDL_current_2100_585_ca","GR_GFDL_current_2100_585_wm","GR_GFDL_2040_2100_585_ca","GR_GFDL_2040_2100_585_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(GR_SpheCoop_src_map_585_GFDL, main = 'GR_GFDL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )

#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
#---------------------------------------------------------
# Species range change setup
# Garden route
Sphe_coop_bin_proj_REU_current_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_REU_2040_126_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_REU_2100_126_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_REU_2040_585_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_REU_2100_585_GFDL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_GFDL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

#---------------------------------------------------------
# Reunion route ssp126
REU_SRC_current_2040_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_GFDL,
    Sphe_coop_bin_proj_REU_2040_126_GFDL
)
REU_SRC_current_2040_126_GFDL$Compt.By.Models

REU_SRC_current_2100_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_GFDL,
    Sphe_coop_bin_proj_REU_2100_126_GFDL
)
REU_SRC_current_2100_126_GFDL$Compt.By.Models

REU_SRC_2040_2100_126_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_2040_126_GFDL,
    Sphe_coop_bin_proj_REU_2100_126_GFDL
)
REU_SRC_2040_2100_126_GFDL$Compt.By.Models

# Reunion route ssp585
REU_SRC_current_2040_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_GFDL,
    Sphe_coop_bin_proj_REU_2040_585_GFDL
)
REU_SRC_current_2040_585_GFDL$Compt.By.Models

REU_SRC_current_2100_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_GFDL,
    Sphe_coop_bin_proj_REU_2100_585_GFDL
)
REU_SRC_current_2100_585_GFDL$Compt.By.Models

REU_SRC_2040_2100_585_GFDL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_2040_585_GFDL,
    Sphe_coop_bin_proj_REU_2100_585_GFDL
)
REU_SRC_2040_2100_585_GFDL$Compt.By.Models

#---------------------------------------------------------
# Plot species range change

# Load libraries to change levelplot theme
install.packages("viridis")
library(viridis)

# Garden route ssp126
REU_raster_2040_126_GFDL <- REU_SRC_current_2040_126_GFDL$Diff.By.Pixel
REU_raster_2100_126_GFDL <- REU_SRC_current_2100_126_GFDL$Diff.By.Pixel
REU_raster_future_126_GFDL <- REU_SRC_2040_2100_126_GFDL$Diff.By.Pixel

plot(REU_raster_2040_126_GFDL)
plot(REU_raster_2100_126_GFDL)
plot(REU_raster_future_126_GFDL)
REU_layers_126_GFDL<-c(REU_raster_2040_126_GFDL,REU_raster_2100_126_GFDL, REU_raster_future_126_GFDL)
summary(REU_layers_126_GFDL)
REU_SpheCoop_src_map_126_GFDL <- stack(REU_layers_126_GFDL)
summary(REU_SpheCoop_src_map_126_GFDL)
names(REU_SpheCoop_src_map_126_GFDL) <- c("REU_GFDL_current_2040_126_ca", "REU_GFDL_current_2040_126_wm","REU_GFDL_current_2100_126_ca","REU_GFDL_current_2100_126_wm","REU_GFDL_2040_2100_126_ca","REU_GFDL_2040_2100_126_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(REU_SpheCoop_src_map_126_GFDL, main = 'REU_GFDL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


# Garden route 585
REU_raster_2040_585_GFDL <- REU_SRC_current_2040_585_GFDL$Diff.By.Pixel
REU_raster_2100_585_GFDL <- REU_SRC_current_2100_585_GFDL$Diff.By.Pixel
REU_raster_future_585_GFDL <- REU_SRC_2040_2100_585_GFDL$Diff.By.Pixel

plot(REU_raster_2040_585_GFDL)
plot(REU_raster_2100_585_GFDL)
plot(REU_raster_future_585_GFDL)
REU_layers_585_GFDL<-c(REU_raster_2040_585_GFDL,REU_raster_2100_585_GFDL, REU_raster_future_585_GFDL)
summary(REU_layers_585_GFDL)
REU_SpheCoop_src_map_585_GFDL <- stack(REU_layers_585_GFDL)
summary(REU_SpheCoop_src_map_585_GFDL)
names(REU_SpheCoop_src_map_585_GFDL) <- c("REU_GFDL_current_2040_585_ca", "REU_GFDL_current_2040_585_wm","REU_GFDL_current_2100_585_ca","REU_GFDL_current_2100_585_wm","REU_GFDL_2040_2100_585_ca","REU_GFDL_2040_2100_585_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(REU_SpheCoop_src_map_585_GFDL, main = 'REU_GFDL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )

#=========================================================
# IPSL-CM6A
#=========================================================
# Garden Route
#---------------------------------------------------------
#---------------------------------------------------------
# Species range change setup
# Garden route
Sphe_coop_bin_proj_GR_current <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_126_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_126_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_GR_2040_585_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_GR_2100_585_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_GR_IPSL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

#---------------------------------------------------------
# Garden route ssp126
GR_SRC_current_2040_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_126_IPSL
)
GR_SRC_current_2040_126_IPSL$Compt.By.Models

GR_SRC_current_2100_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_126_IPSL
)
GR_SRC_current_2100_126_IPSL$Compt.By.Models

GR_SRC_2040_2100_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_126_IPSL,
    Sphe_coop_bin_proj_GR_2100_126_IPSL
)
GR_SRC_2040_2100_126_IPSL$Compt.By.Models

# Garden route ssp585
GR_SRC_current_2040_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2040_585_IPSL
)
GR_SRC_current_2040_585_IPSL$Compt.By.Models

GR_SRC_current_2100_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_current,
    Sphe_coop_bin_proj_GR_2100_585_IPSL
)
GR_SRC_current_2100_585_IPSL$Compt.By.Models

GR_SRC_2040_2100_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_GR_2040_585_IPSL,
    Sphe_coop_bin_proj_GR_2100_585_IPSL
)
GR_SRC_2040_2100_585_IPSL$Compt.By.Models

#---------------------------------------------------------
# Plot species range change

# Load libraries to change levelplot theme
library(viridis)

# Garden route ssp126
GR_raster_2040_126_IPSL <- GR_SRC_current_2040_126_IPSL$Diff.By.Pixel
GR_raster_2100_126_IPSL <- GR_SRC_current_2100_126_IPSL$Diff.By.Pixel
GR_raster_future_126_IPSL <- GR_SRC_2040_2100_126_IPSL$Diff.By.Pixel

plot(GR_raster_2040_126_IPSL)
plot(GR_raster_2100_126_IPSL)
plot(GR_raster_future_126_IPSL)
GR_layers_126_IPSL<-c(GR_raster_2040_126_IPSL,GR_raster_2100_126_IPSL, GR_raster_future_126_IPSL)
summary(GR_layers_126_IPSL)
GR_SpheCoop_src_map_126_IPSL <- stack(GR_layers_126_IPSL)
summary(GR_SpheCoop_src_map_126_IPSL)
names(GR_SpheCoop_src_map_126_IPSL) <- c("GR_IPSL_current_2040_126_ca", "GR_IPSL_current_2040_126_wm","GR_IPSL_current_2100_126_ca","GR_IPSL_current_2100_126_wm","GR_IPSL_2040_2100_126_ca","GR_IPSL_2040_2100_126_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(GR_SpheCoop_src_map_126_IPSL, main = 'SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


# Garden route 585
GR_raster_2040_585_IPSL <- GR_SRC_current_2040_585_IPSL$Diff.By.Pixel
GR_raster_2100_585_IPSL <- GR_SRC_current_2100_585_IPSL$Diff.By.Pixel
GR_raster_future_585_IPSL <- GR_SRC_2040_2100_585_IPSL$Diff.By.Pixel

plot(GR_raster_2040_585_IPSL)
plot(GR_raster_2100_585_IPSL)
plot(GR_raster_future_585_IPSL)
GR_layers_585_IPSL<-c(GR_raster_2040_585_IPSL,GR_raster_2100_585_IPSL, GR_raster_future_585_IPSL)
summary(GR_layers_585_IPSL)
GR_SpheCoop_src_map_585_IPSL <- stack(GR_layers_585_IPSL)
summary(GR_SpheCoop_src_map_585_IPSL)
names(GR_SpheCoop_src_map_585_IPSL) <- c("GR_IPSL_current_2040_585_ca", "GR_IPSL_current_2040_585_wm","GR_IPSL_current_2100_585_ca","GR_IPSL_current_2100_585_wm","GR_IPSL_2040_2100_585_ca","GR_IPSL_2040_2100_585_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(GR_SpheCoop_src_map_585_IPSL, main = 'GR_IPSL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


#---------------------------------------------------------
# Reunion
#---------------------------------------------------------
#---------------------------------------------------------
# Species range change setup
# Garden route
Sphe_coop_bin_proj_REU_current_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_ens_current\\individual_projections\\Sphaeropteris.cooperi_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_REU_2040_126_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2040_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_REU_2100_126_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2100_126\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

Sphe_coop_bin_proj_REU_2040_585_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2040_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))
Sphe_coop_bin_proj_REU_2100_585_IPSL <- stack(c(
    ca="Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif",
    wm = "Sphaeropteris.cooperi\\proj_REU_IPSL_proj_ens_2100_585\\individual_projections\\Sphaeropteris.cooperi_EMcaByTSS_mergedData_mergedRun_mergedAlgo_TSSbin.tif"
))

#---------------------------------------------------------
# Reunion route ssp126
REU_SRC_current_2040_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_IPSL,
    Sphe_coop_bin_proj_REU_2040_126_IPSL
)
REU_SRC_current_2040_126_IPSL$Compt.By.Models

REU_SRC_current_2100_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_IPSL,
    Sphe_coop_bin_proj_REU_2100_126_IPSL
)
REU_SRC_current_2100_126_IPSL$Compt.By.Models

REU_SRC_2040_2100_126_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_2040_126_IPSL,
    Sphe_coop_bin_proj_REU_2100_126_IPSL
)
REU_SRC_2040_2100_126_IPSL$Compt.By.Models

# Reunion route ssp585
REU_SRC_current_2040_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_IPSL,
    Sphe_coop_bin_proj_REU_2040_585_IPSL
)
REU_SRC_current_2040_585_IPSL$Compt.By.Models

REU_SRC_current_2100_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_current_IPSL,
    Sphe_coop_bin_proj_REU_2100_585_IPSL
)
REU_SRC_current_2100_585_IPSL$Compt.By.Models

REU_SRC_2040_2100_585_IPSL <- BIOMOD_RangeSize(
    Sphe_coop_bin_proj_REU_2040_585_IPSL,
    Sphe_coop_bin_proj_REU_2100_585_IPSL
)
REU_SRC_2040_2100_585_IPSL$Compt.By.Models

#---------------------------------------------------------
# Plot species range change

# Load libraries to change levelplot theme
install.packages("viridis")

# Garden route ssp126
REU_raster_2040_126_IPSL <- REU_SRC_current_2040_126_IPSL$Diff.By.Pixel
REU_raster_2100_126_IPSL <- REU_SRC_current_2100_126_IPSL$Diff.By.Pixel
REU_raster_future_126_IPSL <- REU_SRC_2040_2100_126_IPSL$Diff.By.Pixel

plot(REU_raster_2040_126_IPSL)
plot(REU_raster_2100_126_IPSL)
plot(REU_raster_future_126_IPSL)
REU_layers_126_IPSL<-c(REU_raster_2040_126_IPSL,REU_raster_2100_126_IPSL, REU_raster_future_126_IPSL)
summary(REU_layers_126_IPSL)
REU_SpheCoop_src_map_126_IPSL <- stack(REU_layers_126_IPSL)
summary(REU_SpheCoop_src_map_126_IPSL)
names(REU_SpheCoop_src_map_126_IPSL) <- c("REU_IPSL_current_2040_126_ca", "REU_IPSL_current_2040_126_wm","REU_IPSL_current_2100_126_ca","REU_IPSL_current_2100_126_wm","REU_IPSL_2040_2100_126_ca","REU_IPSL_2040_2100_126_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(REU_SpheCoop_src_map_126_IPSL, main = 'REU_IPSL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )


# Garden route 585
REU_raster_2040_585_IPSL <- REU_SRC_current_2040_585_IPSL$Diff.By.Pixel
REU_raster_2100_585_IPSL <- REU_SRC_current_2100_585_IPSL$Diff.By.Pixel
REU_raster_future_585_IPSL <- REU_SRC_2040_2100_585_IPSL$Diff.By.Pixel

plot(REU_raster_2040_585_IPSL)
plot(REU_raster_2100_585_IPSL)
plot(REU_raster_future_585_IPSL)
REU_layers_585_IPSL<-c(REU_raster_2040_585_IPSL,REU_raster_2100_585_IPSL, REU_raster_future_585_IPSL)
summary(REU_layers_585_IPSL)
REU_SpheCoop_src_map_585_IPSL <- stack(REU_layers_585_IPSL)
summary(REU_SpheCoop_src_map_585_IPSL)
names(REU_SpheCoop_src_map_585_IPSL) <- c("REU_IPSL_current_2040_585_ca", "REU_IPSL_current_2040_585_wm","REU_IPSL_current_2100_585_ca","REU_IPSL_current_2100_585_wm","REU_IPSL_2040_2100_585_ca","REU_IPSL_2040_2100_585_wm")
library(rasterVis)
my.at <- seq(-2.5, 1.5,1)
myColorkey <- list(at = my.at, labels = list(labels = c('lost', 'pres', 'abs', 'gain'),at = my.at[- 1]- 0.5 ))
levelplot(REU_SpheCoop_src_map_585_IPSL, main = 'REU_IPSL_SpheCoop range change', colorkey = myColorkey, layout = c(2, 3), par.settings = BuRdTheme )
