# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 05. Ensemble

rm(list = ls())

# Packages
library(raster)
library(tidyverse)
library(vegan)

# Working directory
path <- "G:/My Drive/UMississippi/amphibian_data"
path_in <- file.path(path, "04_algorithms")
path_out <- file.path(path, "05_ensemble")

if(!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)

# Parameters
auc_limit <- 0.75
cenarios <- c("correlative_acc17", "correlative_acc25", 
              "acc17_infected", "acc17_noninfected", 
              "acc25_infected", "acc25_noninfected")
tempos <- c("present", "245", "585")

for(cenario in cenarios) {
  
  message(">>> STARTING ENSEMBLE FOR: ", cenario)
 
  
  dir_rep <- file.path(path_in, cenario, "00_replicas")
  dir_eval <- file.path(path_in, cenario, "01_evaluation", "00_raw")
  
  # List all evaluation spreadsheets generated for this scenario (one per species)
  csv_files <- list.files(dir_eval, pattern = "evaluation_.*\\.csv$", full.names = TRUE)
  
  if(length(csv_files) == 0) {
    message("  [!] No evaluations found for ", cenario, ". Skipping...")
    next
  }
  
  eva_all <- purrr::map_dfr(csv_files, readr::read_csv, show_col_types = FALSE)
  
  for(sp in unique(eva_all$species)) {
    
    # Selecting models with AUC >= 0.75
    eva_i <- eva_all %>% 
      dplyr::filter(species == sp, auc >= auc_limit)
    
    if(nrow(eva_i) == 0) {
      message("  [-] ", sp, " did not reach minimum AUC. Skipping...")
      next
    }
    
    message("  -> Processing species: ", sp)
    
    for(tp in tempos) {
      
      files_to_stack <- c()
      weights_to_use <- c()
      

      # RASTERS AND WEIGHTS SELECTION
     
      if(tp == "present") {
        # In the present, just take the files listed in the evaluation table
        files_to_stack <- file.path(dir_rep, eva_i$file)
        weights_to_use <- (eva_i$auc - 0.5) ^ 2
        
      } else {
        # In the future (245 or 585), we create the 3 GCMs for each good model from the present
        for(row in 1:nrow(eva_i)) {
          alg <- eva_i$algorithm[row]
          rep_num <- eva_i$replica[row]
          rep_str <- ifelse(rep_num < 10, paste0("0", rep_num), as.character(rep_num))
          peso_futuro <- (eva_i$auc[row] - 0.5) ^ 2
          
          for(gcm in c("had", "mir", "mpi")) {
            # Construct the name exactly as the generator script saved it
            f_name <- paste0(cenario, "_", gcm, tp, "_", sp, "_", alg, "_r", rep_str, ".tif")
            f_path <- file.path(dir_rep, f_name)
            
            if(file.exists(f_path)) {
              files_to_stack <- c(files_to_stack, f_path)
              weights_to_use <- c(weights_to_use, peso_futuro)
            }
          }
        }
      }
      
      if(length(files_to_stack) == 0) next
      
      # Importing the models
      enm <- raster::stack(files_to_stack)
      v <- raster::values(enm)
      

      # STANDARDIZATION (Robustness adjustment)

      
      # Normalize each layer individually before combining
      enm_st <- apply(v, 2, function(x) {
        if(all(is.na(x))) return(rep(NA, length(x)))
        
        min_x <- min(x, na.rm = TRUE)
        max_x <- max(x, na.rm = TRUE)
        
        # If min == max (constant layer), return 0 or NA not to influence
        if(max_x == min_x) return(rep(0, length(x)))
        
        return((x - min_x) / (max_x - min_x))
      })
      
      # Filter weights (aucs) to remove failed models (zero range sum)
      # If a layer was all zero, it should not enter the ensemble
      is_valid <- apply(enm_st, 2, function(x) sum(x, na.rm = TRUE) > 0)
      
      if(sum(is_valid) == 0) {
        message("    [!] All layers resulted in zero or constant suitability.")
        next
      }
      
      enm_st <- enm_st[, is_valid, drop = FALSE]
      auc_final <- weights_to_use[is_valid]
      

      # ENSEMBLE

      ens <- enm[[1]]
      ens[] <- apply(enm_st, 1, function(x) {
        if(all(is.na(x))) return(NA)
        sum(x * auc_final, na.rm = TRUE) / sum(auc_final, na.rm = TRUE)
      })
      
      # Exporting the ensemble
      out_name <- paste("ensemble", cenario, tp, sp, sep = "_")
      
      raster::writeRaster(x = ens, 
                          filename = file.path(path_out, paste0(out_name, ".tif")), 
                          format = "GTiff", 
                          options = c("COMPRESS=DEFLATE"), 
                          overwrite = TRUE)
      
      message("    [OK] Saved: ", out_name)
    }
  }
}

message(">>> ENSEMBLES COMPLETED.")
