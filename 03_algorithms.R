# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 03. Algorithms

rm(list = ls())

library(dismo)
library(kernlab)
library(raster)
library(sf)
library(tidyverse)
library(randomForest)
library(maxnet)

# --- PATH CONFIGURATIONS ---
path <- "G:/My Drive/UMississippi/amphibian_data"
path_vars <- file.path(path, "03_var")
path_out <- file.path(path, "04_algorithms")

# Buffer for automated cropping
path_mask <- file.path(path, "buffer_accessible_area_allpoints.gpkg")
shape_mask <- sf::st_read(path_mask)

occ <- readr::read_csv(file.path(path, "02_occ/coord_atelopus_clean.csv"))

replica <- 10 
partition <- .7 

gcms <- c("had", "mir", "mpi")
ssps <- c("245", "585")
scenarios <- c(
  "correlative_acc17", "correlative_acc25",
  "acc17_infected", "acc17_noninfected",
  "acc25_infected", "acc25_noninfected"
)

# Function to automatically load and align rasters
load_aligned_stack <- function(file_list, template_raster) {
  list_ras <- lapply(file_list, raster::raster)
  list_adjusted <- lapply(list_ras, function(r) {
    if (!raster::compareRaster(r, template_raster, stopiffalse = FALSE)) {
      return(raster::resample(r, template_raster, method = "ngb"))
    } else {
      return(r)
    }
  })
  s <- raster::brick(raster::stack(list_adjusted))
  names(s) <- tools::file_path_sans_ext(basename(file_list))
  return(s)
}

# Function that handles prediction for each algorithm
pred_wrapper <- function(model, stack, algorithm_name) {
  if (algorithm_name == "maxent") {
    # 1. template with valid values
    template <- stack[[1]]
    # 2. predict only in pixels that are not NA
    pred_vals <- predict(model, raster::as.data.frame(stack), type = "logistic")
    r <- template
    r[] <- NA
    r[!is.na(template[])] <- pred_vals
    
    return(r)
    
  } else {
    return(dismo::predict(stack, model))
  }
}

for (current_scenario in scnearios) {
  
  # uncomment this if you do not want to run the correlative scenario again
  # if (grepl("correlative", current_scenario)) {
  #   message(">>> Skipping scenario: ", current_scenario)
  #   next
  # }
  
  message("\n=======================================================")
  message(">>> PROCESSING SCENARIO: ", current_scenario)
  message("=======================================================")
  
  dir_rep <- file.path(path_out, current_scenario, "00_replicas")
  dir_eval <- file.path(path_out, current_scenario, "01_evaluation", "00_raw")
  dir_imp <- file.path(path_out, current_scenario, "02_importance") # Folder for the single CSV
  
  dir.create(dir_rep, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_eval, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_imp, recursive = TRUE, showWarnings = FALSE)
  
  is_corr <- grepl("correlative", current_scenario)
  acc_base <- str_extract(current_scenario, "acc\\d+")
  inf_base <- if(is_corr) "infected" else str_extract(current_scenario, "(infected|noninfected)")
  
  baseline_folder <- file.path(path_vars, paste(acc_base, inf_base, sep="_"))
  clim_files <- list.files(baseline_folder, pattern = "\\.tif$", full.names = TRUE, recursive = FALSE)
  
  if (!is_corr) {
    perf_files <- list.files(file.path(baseline_folder, "performance"), pattern = "\\.tif$", full.names = TRUE)
    all_files <- c(clim_files, perf_files)
  } else {
    all_files <- clim_files
  }
  
  var_present <- load_aligned_stack(all_files, raster::raster(all_files[1]))
  names(var_present) <- tools::file_path_sans_ext(basename(all_files))
  
  for(i in occ$species %>% unique) {
    
    final_evaluation <- file.path(dir_eval, paste0("evaluation_", current_scenario, "_", i, ".csv"))
    
    # Check: if evaluation already exists, skip it (to avoid repeating completed runs)
    if (file.exists(final_evaluation)) {
      message("  [!] ", i, " already finished. Skipping...")
      next 
    }
    
    eval_algorithm <- tibble::tibble()
    imp_algorithm <- tibble::tibble() # Importance accumulator
    
    # Sampling (Entire Panama)
    pr_species <- occ %>% filter(species == i) %>% dplyr::select(lon, lat) %>% mutate(id = seq(nrow(.)))
    bg_species <- dismo::randomPoints(mask = var_present, n = nrow(pr_species)*10) %>% as_tibble() %>% rename(lon = x, lat = y) %>% mutate(id = seq(nrow(.)))
    pa_species <- dismo::randomPoints(mask = var_present, n = 10000) %>% as_tibble() %>% rename(lon = x, lat = y) %>% mutate(id = seq(nrow(.)))
    
    for(r in 1:replica){  
      message("\n[SERIES: ", i, "] - Replica: ", r, "/", replica)
      
      pr_sample_train <- pr_species %>% sample_frac(partition) %>% pull(id)
      bg_sample_train <- bg_species %>% sample_frac(partition) %>% pull(id)
      pa_sample_train <- pa_species %>% sample_frac(partition) %>% pull(id)
      
      train <- dismo::prepareData(x = var_present, p = pr_species[pr_species$id %in% pr_sample_train, c("lon","lat")], b = bg_species[bg_species$id %in% bg_sample_train, c("lon","lat")]) %>% na.omit
      test  <- dismo::prepareData(x = var_present, p = pr_species[!pr_species$id %in% pr_sample_train, c("lon","lat")], b = bg_species[!bg_species$id %in% bg_sample_train, c("lon","lat")]) %>% na.omit
      train_pa <- dismo::prepareData(x = var_present, p = pr_species[pr_species$id %in% pr_sample_train, c("lon","lat")], b = pa_species[pa_species$id %in% pa_sample_train, c("lon","lat")]) %>% na.omit
      
      # 2. Fit models
      fit <- list(
        maxent = maxnet::maxnet(p = train$pb, data = train[, -1]),
        Domain = dismo::domain(x = train[train$pb == 1, -1]),
        glm    = glm(formula = pb ~ ., family = binomial(link = "logit"), data = train_pa),
        svm    = kernlab::ksvm(x = pb ~ ., data = train)
      )
      
      for(a in seq(fit)){
        name_alg <- names(fit)[a]
        message("  -> Algorithm: ", toupper(name_alg))
        
        # --- 1. IMPORTANCE (Calculation and Accumulation in Memory) ---
        auc_ref <- dismo::evaluate(p = test[test$pb == 1, -1], a = test[test$pb == 0, -1], model = fit[[a]])@auc
        
        for(var_name in names(var_present)){
          test_perm <- test
          test_perm[, var_name] <- sample(test_perm[, var_name])
          auc_perm <- dismo::evaluate(p = test_perm[test_perm$pb == 1, -1], a = test_perm[test_perm$pb == 0, -1], model = fit[[a]])@auc
          
          imp_algorithm <- bind_rows(imp_algorithm, tibble(
            species = i, replica = r, algorithm = name_alg, 
            variable = var_name, importance = max(0, auc_ref - auc_perm)
          ))
        }
        
        # --- 2. PREDICTIONS + AUTOMATED CROP ---
        p_pres <- pred_wrapper(fit[[a]], var_present, name_alg)
        p_pres_crop <- raster::mask(raster::crop(p_pres, shape_mask), shape_mask)
        
        name_pres <- paste0(current_scenario, "_", i, "_", name_alg, "_r", sprintf("%02d", r))
        raster::writeRaster(p_pres_crop, file.path(dir_rep, paste0(name_pres, ".tif")), format="GTiff", options=c("COMPRESS=DEFLATE"), overwrite=T)
        
        for(gcm in gcms) {
          for(ssp in ssps) {
            # Load future variables
            clim_f <- list.files(file.path(path_vars, paste(gcm, ssp, acc_base, sep="_")), pattern = "\\.tif$", full.names = T)
            if (!is_corr) {
              perf_f <- list.files(file.path(path_vars, paste(acc_base, inf_base, ssp, sep="_")), pattern = "\\.tif$", full.names = T)
              all_f_fut <- c(clim_f, perf_f) # Complete list
            } else { 
              all_f_fut <- clim_f 
            }
            
            # ALIGN WITH PRESENT (var_present is already corrected!)
            var_f <- load_aligned_stack(all_f_fut, var_present[[1]])
            
            names(var_f) <- names(var_present)
            p_fut <- dismo::predict(var_f, fit[[a]])
            p_fut_crop <- raster::mask(raster::crop(p_fut, shape_mask), shape_mask)
            
            name_fut <- paste0(current_scenario, "_", gcm, ssp, "_", i, "_", name_alg, "_r", sprintf("%02d", r))
            raster::writeRaster(p_fut_crop, file.path(dir_rep, paste0(name_fut, ".tif")), format="GTiff", options=c("COMPRESS=DEFLATE"), overwrite=T)
          }
        }
        
        # --- 3. EVALUATION ---
        # --- 4. COMPLETE EVALUATION ---
        eval <- dismo::evaluate(p = test[test$pb == 1, -1], a = test[test$pb == 0, -1], model = fit[[a]])
        
        # Extract the index of the threshold that maximizes Sensitivity + Specificity
        idx <- which.max(eval@TPR + eval@TNR)
        
        # Calculate TSS and capture Threshold
        tss_val   <- eval@TPR[idx] + eval@TNR[idx] - 1
        threshold <- eval@t[idx]
        
        # Save everything in the accumulator
        eval_algorithm <- bind_rows(eval_algorithm, tibble(
          species   = i, 
          replica   = r, 
          algorithm = name_alg, 
          auc       = eval@auc, 
          tss       = tss_val,         # The missing TSS
          threshold = threshold,       # Useful to create the binary map later
          file      = paste0(name_pres, ".tif")
        ))
      } 
    } 
    # Save only TWO CSVs per scenario: one for Evaluation and one for Importance
    readr::write_csv(eval_algorithm, final_evaluation)
    readr::write_csv(imp_algorithm, file.path(dir_imp, paste0("importance_", current_scenario, ".csv")))
    message("  [OK] Consolidated CSVs saved for ", current_scenario)
  } 
}
