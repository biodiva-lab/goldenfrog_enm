# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 01. Crop Raster

rm(list = ls())
library(terra)
library(sf)

# 1. Load the buffer
path <- "G:/My Drive/UMississippi/amphibian_data"
area <- st_read(file.path(path, "buffer_accessible_area_allpoints.gpkg"))

# 2. List all rasters (including the ones just generated)
base_dir <- file.path(path, "03_var")
rasters <- list.files(base_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

# Remove files that have already been cropped previously to avoid processing twice
rasters <- rasters[!grepl("_crop\\.tif$", rasters)]

# 3. Create a Template (Empty Raster) with the dimensions of your buffer
# This ensures that all cropped rasters have the same grid format at the end
template <- rast(vect(area), res = 0.008333333) # 30 arc-seconds resolution (~1km)

for (f in rasters) {
  message("Processing: ", basename(f))
  
  r <- rast(f)
  
  # Reproject the vector to the raster's CRS
  area_proj <- vect(st_transform(area, crs(r)))
  
  # Step 1: Crop and Mask
  r_crop <- mask(crop(r, area_proj), area_proj)
  
  # Step 2: Alignment (The trick to avoid extent mismatch errors)
  # This forces the raster to have the exact same extent as the template
  r_final <- resample(r_crop, template, method = "near")
  
  # 4. Save
  out_path <- gsub("\\.tif$", "_crop.tif", f)
  writeRaster(r_final, out_path, overwrite = TRUE)
}

message("All rasters have been cropped and standardized to the buffer!")
