# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 02. Data Cleaning

rm(list = ls())

library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(wallace)
library(ggplot2)
library(tidyverse)

# Working directory
path <- "G:/My Drive/UMississippi/amphibian_data"
setwd(path)
setwd("02_occ")
dir()

# Occurrence data
results <- read.csv("coord_atelopus_clean.csv", header = TRUE)
summary(results)

# Excluding NAs
occ_data_na <- results %>% 
  tidyr::drop_na(lon, lat)
summary(occ_data_na)

# Cleaning process - flag problematic points
flags_spatial <- CoordinateCleaner::clean_coordinates(
  x = occ_data_na, 
  species = "species",
  lon = "lon", 
  lat = "lat",
  tests = c(
    # "capitals",   # 10km buffer around capitals
    # "centroids",  # 1km buffer around country centroids
    "duplicates", # Duplicates
    "equal",      # Equal coordinates
    # "institutions", # Buffer around biodiversity institutions
    "seas",       # Sea points
    # "urban",      # Urban areas
    "validity",   # Points outside the coordinate system
    "zeros"       # Zeros and points where lat=lon
  )
)

#' TRUE = clean coordinates
#' FALSE = potentially problematic coordinates

head(flags_spatial)
summary(flags_spatial)

# Exclusion of flagged points
occ_data_tax_date_spa <- occ_data_na %>% 
  dplyr::filter(flags_spatial$.summary == TRUE)
occ_data_tax_date_spa

# Data summary
table(occ_data_na$species)
table(occ_data_tax_date_spa$species)

# Second part of the cleaning process
# Keeping points within our area of interest
setwd(path)
setwd("03_var_30/acc17_infected")

# Select one of the rasters
var_id <- raster::raster("raster_br_res30_wc2.0_bio_30s_01.tif")
var_id

# Extract coordinates from the raster
var_id[!is.na(var_id)] <- raster::cellFromXY(var_id, raster::rasterToPoints(var_id)[, 1:2])

# Filter points that are within the raster coordinates
occ_data_tax_date_spa_oppc <- occ_data_tax_date_spa %>% 
  dplyr::mutate(oppc = raster::extract(var_id, dplyr::select(., lon, lat))) %>% 
  # dplyr::distinct(species, oppc, .keep_all = TRUE) %>% 
  dplyr::filter(!is.na(oppc)) %>% 
  dplyr::add_count(species) %>% 
  dplyr::arrange(species)
occ_data_tax_date_spa_oppc

# Check
table(occ_data_tax_date_spa$species)
table(occ_data_tax_date_spa_oppc$species)

# Export the new data
setwd(path)
dir.create("02_occ", showWarnings = FALSE)
setwd("02_occ")

readr::write_csv(occ_data_tax_date_spa_oppc, "occ_spocc_filtros_taxonomico_data_espatial_oppc.csv")

# End of the cleaning process