library(raster)
library(ggplot2)

# Load raster data MRI-ESM2-0 ssp585 2041-2070
setwd("C:\\Users\\chris\\Desktop\\Megan\\Meg full git repo")

raster_data1 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data4 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data10 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data11 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data15 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data16 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")
raster_data17 <- raster("input_data\\MRI world data\\MRI SSP585 2041-2070\\CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1.tif")

#---------------------------------------------------------
# MRI Garden Route clipping MRI-ESM2-0 ssp585 2041-2070
#---------------------------------------------------------
# Load mask layer of the Garden Route
mask_data <- shapefile("input_data\\Garden_route_mask.shp")

# Crop the rasters by the global mask
cropped_raster1 <- mask(raster_data1, mask_data)
cropped_raster4 <- mask(raster_data4, mask_data)
cropped_raster10 <- mask(raster_data10, mask_data)
cropped_raster11 <- mask(raster_data11, mask_data)
cropped_raster15 <- mask(raster_data15, mask_data)
cropped_raster16 <- mask(raster_data16, mask_data)
cropped_raster17 <- mask(raster_data17, mask_data)

# Define the output file path and names
output_path <- "input_data\\climate_data\\MRI-ESM2-0\\Garden Route\\ssp585\\2041-2070"

output_name1 <- "CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name4 <- "CHELSA_bio4_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name10 <- "CHELSA_bio10_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name11 <- "CHELSA_bio11_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name15 <- "CHELSA_bio15_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name16 <- "CHELSA_bio16_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"
output_name17 <- "CHELSA_bio17_2041-2070_mri-esm2-0_ssp585_V.2.1_GR.tif"

# Save the cropped rasters to new files
writeRaster(cropped_raster1, filename = file.path(output_path, output_name1), format = "GTiff")
writeRaster(cropped_raster4, filename = file.path(output_path, output_name4), format = "GTiff")
writeRaster(cropped_raster10, filename = file.path(output_path, output_name10), format = "GTiff")
writeRaster(cropped_raster11, filename = file.path(output_path, output_name11), format = "GTiff")
writeRaster(cropped_raster15, filename = file.path(output_path, output_name15), format = "GTiff")
writeRaster(cropped_raster16, filename = file.path(output_path, output_name16), format = "GTiff")
writeRaster(cropped_raster17, filename = file.path(output_path, output_name17), format = "GTiff")

#---------------------------------------------------------
# MRI Reunion clipping MRI-ESM2-0 ssp585 2041-2070
#---------------------------------------------------------
# Load mask layer of Reunion
mask_data <- shapefile("input_data\\REU_adm0.shp") 

# Crop the rasters by the global mask
cropped_raster1 <- mask(raster_data1, mask_data)
cropped_raster4 <- mask(raster_data4, mask_data)
cropped_raster10 <- mask(raster_data10, mask_data)
cropped_raster11 <- mask(raster_data11, mask_data)
cropped_raster15 <- mask(raster_data15, mask_data)
cropped_raster16 <- mask(raster_data16, mask_data)
cropped_raster17 <- mask(raster_data17, mask_data)

# Define the output file path and names
output_path <- "input_data\\climate_data\\MRI-ESM2-0\\Reunion\\585\\2041-2070"

output_name1 <- "CHELSA_bio1_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name4 <- "CHELSA_bio4_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name10 <- "CHELSA_bio10_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name11 <- "CHELSA_bio11_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name15 <- "CHELSA_bio15_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name16 <- "CHELSA_bio16_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"
output_name17 <- "CHELSA_bio16_2041-2070_mri-esm2-0_ssp585_V.2.1_R.tif"

# Save the cropped rasters to new files
writeRaster(cropped_raster1, filename = file.path(output_path, output_name1), format = "GTiff")
writeRaster(cropped_raster4, filename = file.path(output_path, output_name4), format = "GTiff")
writeRaster(cropped_raster10, filename = file.path(output_path, output_name10), format = "GTiff")
writeRaster(cropped_raster11, filename = file.path(output_path, output_name11), format = "GTiff")
writeRaster(cropped_raster15, filename = file.path(output_path, output_name15), format = "GTiff")
writeRaster(cropped_raster16, filename = file.path(output_path, output_name16), format = "GTiff")
writeRaster(cropped_raster17, filename = file.path(output_path, output_name17), format = "GTiff")
