# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 06. Binarization

rm(list = ls())

library(raster)
library(tidyverse)
library(biomod2)

# PATH CONFIGURATIONS
path <- "G:/My Drive/UMississippi/amphibian_data"
path_ens <- file.path(path, "05_ensemble_30")
path_bin <- file.path(path, "06_binary_30")
dir.create(path_bin, showWarnings = FALSE, recursive = TRUE)

# 1. Occurrence data
occ <- read.csv(file.path(path, "02_occ/coord_atelopus_clean.csv"))
occ_pts <- occ[, 1:3] # lon, lat, species
occ_pts$presabs <- 1

# 2. List all Ensembles ('ensemble_' pattern)
ens_files <- list.files(path_ens, pattern = "^ensemble_.*Atelopus.*\\.tif$", full.names = TRUE)

if(length(ens_files) == 0) stop("No 'ensemble_' files found in 05_ensemble_30!")

# Identify base scenarios (removing prefix, time suffixes, and species)
# Ex: transforms 'ensemble_acc17_infected_present_Atelopus_zeteki' into 'acc17_infected'
scenarios <- ens_files %>% 
  basename() %>% 
  str_remove("^ensemble_") %>% 
  str_remove("_(present|245|585).*") %>% 
  unique()

species_list <- unique(occ$species)

for(scen in scenarios) {
  for(sp in species_list) {
    
    message("\n>>> FOCUS ON SCENARIO: ", scen, " | SPECIES: ", sp)
    
    # STEP 1: FIND THE PRESENT MAP TO CALIBRATE THE CUTOFF
    pattern_pres <- paste0("ensemble_", scen, "_present_", sp)
    file_pres <- ens_files[grepl(pattern_pres, ens_files)]
    
    if(length(file_pres) == 0) {
      message("⚠️ Present map not found for pattern: ", pattern_pres)
      next
    }
    
    enm_pres <- raster::raster(file_pres)
    
    # STEP 2: GENERATE PSEUDO-ABSENCES (FROM YOUR SCRIPT'S LOGIC)
    set.seed(42)
    pa_points <- dismo::randomPoints(mask = enm_pres, n = nrow(occ_pts[occ_pts$species == sp, ]) * 10) %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(lon = x, lat = y)
    
    pa_points$species <- sp
    pa_points$presabs <- 0
    
    # Combine species presences and absences
    obs <- dplyr::bind_rows(occ_pts[occ_pts$species == sp, ], pa_points)
    
    # Extract suitability values
    thrs_values <- raster::extract(enm_pres, obs[, c("lon", "lat")])
    obs$thrs <- thrs_values
    
    # Clean NAs (points outside the raster)
    obs_clean <- obs %>% filter(!is.na(thrs))
    fit <- obs_clean$thrs
    real <- obs_clean$presabs
    
    # STEP 3: FIND THE BEST CUTOFF (LOOP FROM 0 TO 1)
    message("    Calculating best cutoff via biomod2...")
    n_seq <- seq(0, 1, by = 0.001)
    tss_list <- NULL
    
    for(i in 1:length(n_seq)){
      tss_step <- biomod2::bm_FindOptimStat(
        metric.eval = "TSS",
        obs = real,
        fit = fit,
        threshold = n_seq[i]
      )
      tss_list <- rbind(tss_list, tss_step)
    }
    
    melhor_cutoff <- tss_list$cutoff[which.max(tss_list$best.stat)]
    message("    🎯 Best Cutoff found: ", melhor_cutoff)
    
    # --- STEP 4: BINARIZE ALL TIME PERIODS FOR THIS SCENARIO ---
    # Search for all files matching the scenario and species (present + futures)
    pattern_all <- paste0("ensemble_", scen, ".*_", sp)
    files_to_bin <- ens_files[grepl(pattern_all, ens_files)]
    
    for(f in files_to_bin) {
      r_to_bin <- raster::raster(f)
      
      # Apply binarization
      r_bin <- r_to_bin >= melhor_cutoff
      
      # Output file name (Ex: BIN_acc17_infected_present_MAX_TSS.tif)
      out_name <- basename(f) %>% 
        str_replace("ensemble_", "BIN_") %>% 
        str_replace("\\.tif$", paste0("_MAX_TSS.tif"))
      
      raster::writeRaster(r_bin, 
                          filename = file.path(path_bin, out_name),
                          format = "GTiff",
                          options = c("COMPRESS=DEFLATE"),
                          overwrite = TRUE)
      
      message("    ✅ Generated: ", out_name)
    }
  }
}

message("\n🎉 COMPLETED! All binary maps are located in: 06_binary_30")