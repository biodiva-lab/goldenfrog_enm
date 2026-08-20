# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 00. Mechanistic approach 


# Mapinguari 
# Load packages
library(magraitr)
library(landscapetools)
library(psych)
library(raster)
library(rgdal)
library(rnaturalearth)
library(tidyverse)
library(wesanderson)

# Load temperature rasters
setwd("C:/Users/luisa/Downloads")
tif <- dir(pattern = "tif$")
tif

# Import the rasters and group them into a stack
var <- raster::stack("wc2.1_30s_tmax_MIROC6_ssp245_2041-2060.tif")
names(var)
var <- var[[12]]

br <- rnaturalearth::ne_countries(country = "Panama", returnclass = "sf")
br

# Crop and mask the raster
var_br <- raster::crop(x = var, y = br) %>% 
  raster::mask(mask = br)
var_br
plot(var_br)
names(var_br) <- paste0("tavg_", 12)

# Load ecophysiological and microclimatic data
setwd("G:/Meu Drive/UMississippi/amphibian_data/Mapinguari")

tpc.geral <- read.csv("tpc.geral.csv", sep = ",")
ecofisiofinal <- read.csv("ecofisiofinal.csv", sep = ";")
Thab <- read.csv("Thab.csv", sep = ";")
Toper <- read.csv("Toper.csv", sep = ";")
tpc.geral <- tpc_geral

# Performance models (GAMMs)
perf_Neg_17_SORA <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Neg_17_SORA",])

perf_Neg_25_SORA <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs', k=5), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Neg_25_SORA",])

perf_Pos_17_SORA <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Pos_17_SORA",])

perf_Pos_25_SORA <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Pos_25_SORA",])

perf_Neg_25_AHOS <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Neg_25_AHOS",])

perf_Neg_17_AHOS <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Neg_17_AHOS",])

perf_Pos_17_AHOS <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Pos_17_AHOS",])

perf_Pos_25_AHOS <-
  mgcv::gamm(rate ~ s(temp, bs = 'cs'), 
             random = list(id = ~ 1), 
             data = tpc.geral[tpc.geral$species == "Pos_25_AHOS",])

# Prepare performance predictions
pred_perf <-
  Mapinguari::get_predict(list(perf_Neg_17_SORAFUN = perf_Neg_17_SORA$gam,
                               perf_Neg_25_SORAFUN = perf_Neg_25_SORA$gam,
                               perf_Pos_17_SORAFUN = perf_Pos_17_SORA$gam,
                               perf_Pos_25_SORAFUN = perf_Pos_25_SORA$gam,
                               perf_Neg_17_AHOSFUN = perf_Neg_17_AHOS$gam,
                               perf_Neg_25_AHOSFUN = perf_Neg_25_AHOS$gam,
                               perf_Pos_17_AHOSFUN = perf_Pos_17_AHOS$gam,
                               perf_Pos_25_AHOSFUN = perf_Pos_25_AHOS$gam))

# Calculate performance rasters
perf_raster <-
  Mapinguari::transform_rasters(var_br,
                                perf_Neg_17_SORA=pred_perf$perf_Neg_17_SORAFUN(temp = tavg),
                                perf_Neg_25_SORA=pred_perf$perf_Neg_25_SORAFUN(temp = tavg),
                                perf_Pos_17_SORA=pred_perf$perf_Pos_17_SORAFUN(temp = tavg),
                                perf_Pos_25_SORA=pred_perf$perf_Pos_25_SORAFUN(temp = tavg),
                                perf_Neg_17_AHOS=pred_perf$perf_Neg_17_AHOSFUN(temp = tavg),
                                perf_Neg_25_AHOS=pred_perf$perf_Neg_25_AHOSFUN(temp = tavg),
                                perf_Pos_17_AHOS=pred_perf$perf_Pos_17_AHOSFUN(temp = tavg),
                                perf_Pos_25_AHOS=pred_perf$perf_Pos_25_AHOSFUN(temp = tavg))

raster::writeRaster(perf_raster, "perf_raster_ssp245_2050.tif")

# Save and plot individual layers
for (i in 1:nlayers(perf_raster)) {
  layer_name <- names(perf_raster)[i]
  file_path <- paste0(layer_name, "ssp245_2050.tif")
  writeRaster(perf_raster[[i]], filename = file_path, format = "GTiff", overwrite = TRUE)
  plot(perf_raster[[i]])
}

# Calculate Vtmax and Vtmin for each species
vtmax_ni <- max(ecofisiofinal[ecofisiofinal$Sp == "non-infected",]$Temp.corp, na.rm = TRUE)
vtmin_ni <- min(ecofisiofinal[ecofisiofinal$Sp == "non-infected",]$Temp.corp, na.rm = TRUE)

vtmax_i <- max(ecofisiofinal[ecofisiofinal$Sp == "infected",]$Temp.corp, na.rm = TRUE)
vtmin_i <- min(ecofisiofinal[ecofisiofinal$Sp == "infected",]$Temp.corp, na.rm = TRUE)

# Calculate mean daily temperatures for each area
Thab_day <-
  Thab %>%
  dplyr::group_by(Area, Year, Month, Day) %>% 
  dplyr::summarise(Temp = mean(Temp))

# Format operative temperature table
Toper_C1 <- Toper[1:11]
Toper_C2 <- Toper[c(1:9, 12, 13)]

names(Toper_C1)[10:11] <- c("Temp", "Microhabitat")
names(Toper_C2)[10:11] <- c("Temp", "Microhabitat")

Toper_long <- rbind(Toper_C1, Toper_C2)

# Calculate activity hours
hvtFUN_ni <- function(x) ifelse(x > vtmin_ni & x < vtmax_ni, 1, 0)
hvtFUN_i <- function(x) ifelse(x > vtmin_i & x < vtmax_i, 1, 0)

Toper_long$hvt_ni <- hvtFUN_ni(Toper_long$Temp)
Toper_long$hvt_i <- hvtFUN_i(Toper_long$Temp)

hvt_day <-
  Toper_long %>%
  dplyr::group_by(Area, Datalogger, Microhabitat, Year, Month, Day) %>% 
  dplyr::summarise(hvt_ni = sum(hvt_ni)/30,
                   hvt_i = sum(hvt_i)/30)

hvt_df <- 
  dplyr::left_join(hvt_day, 
                   Thab_day,
                   by = c("Area", "Year", "Month", "Day"))

hvt_df <- hvt_df[complete.cases(hvt_df),]

# Create model relating activity hours to air temperature
hvt_logistic_ni <- nls(hvt_ni + 0.00001 ~ SSlogis(Temp, Asym, xmid, scal), data = hvt_df)
hvt_logistic_i <- nls(hvt_i + 0.00001 ~ SSlogis(Temp, Asym, xmid, scal), data = hvt_df)

pred_hvt <-
  Mapinguari::get_predict(list(hvt_niFUN = hvt_logistic_ni,
                               hvt_iFUN = hvt_logistic_i))

# Calculate activity hours raster
hvt_raster <-
  Mapinguari::transform_rasters(var_br,
                                hvt_ni = pred_hvt$hvt_niFUN(Temp = tavg),
                                hvt_i = pred_hvt$hvt_iFUN(Temp = tavg))

setwd("Mapinguari")
writeRaster(hvt_raster, "hvt_raster_ssp245_2050.tif", overwrite = TRUE)

for (i in 1:nlayers(hvt_raster)) {
  layer_name <- names(hvt_raster)[i]
  file_path <- paste0(layer_name, "ssp245_2050.tif")
  writeRaster(hvt_raster[[i]], filename = file_path, format = "GTiff", overwrite = TRUE)
  plot(hvt_raster[[i]])
}


## NicheMapR

rm(list = ls())

library(NicheMapR)
library(terra)
library(dplyr)
library(ranger)

# ## if we need to create the microclim data
# get.global.climate() 
# setwd("G:/My Drive/UMississippi/amphibian_data")
# loc2 <- read.csv("G:/My Drive/UMississippi/amphibian_data/02_occ/coord_atelopus_clean.csv", header = TRUE) 
# # call the microclimate model, global climate database implementation
# run_microclim <- function(lon, lat) {
#   micro_global(loc = c(lon, lat),
#                timeinterval = 365,  # generates one value per day
#                nyears = 1)          # simulates 1 year
# }

# micro <- mapply(run_microclim, loc2$lon, loc2$lat, SIMPLIFY = FALSE)
# 
# micro_list <- lapply(1:nrow(loc2), function(i) {
#   
#   # 1. extract lat and lon as a numeric vector 
#   coordenada_atual <- c(as.numeric(loc2$lon[i]), as.numeric(loc2$lat[i]))
#   
#   message(sprintf("Running micro_global for point %d of %d...", i, nrow(loc2)))
#   
#   # 2. run the function for each coordinate
#   micro_sim <- tryCatch({
#     micro_global(loc = coordenada_atual, writecsv = 0)
#   }, error = function(e) {
#     message(sprintf("Error at coordinate %d: %s", i, e$message))
#     return(NULL)
#   })
#   
#   return(micro_sim)
# })
# 
# save(micro_list, file = 'micro_atelopus.Rda')


#### Ectotherm

# --- 1. CONFIGURATION AND PARAMETERS ---
setwd("G:/My Drive/UMississippi/amphibian_data")
loc2 <- read.csv("02_occ/coord_atelopus_clean.csv", header = TRUE)
loc2$point <- 1:nrow(loc2)
load('micro_atelopus.Rda') # Loads micro_list

# Morphological and behavioral variables
alpha_min <- 0.80; alpha_max <- 0.85; shape <- 4
Ww_g_i <- 7.4; Ww_g_ni <- 7.4
T_pref_i <- 21.2; T_pref_ni <- 21.2
T_RB_min <- 16; T_B_min <- 16; T_F_min <- 18; T_F_max <- 32
mindepth <- 1; maxdepth <- 2
shade_seek <- 1; burrow <- 0; climb <- 1
nocturn <- 0; crepus <- 0; diurn <- 1
aquabask <- 0; pct_wet <- 1; z.mult <- 1; alpha_sub <- 0.7

# DEB and reproduction variables
V_init <- 3e-9; E.0 <- 1e-3; E_init <- E.0/V_init; E_H_init <- 0; stage <- 0
viviparous <- 0; clutchsize <- 300; photostart <- 1; photofinish <- 2

# CT parameters
CT_max_acc17_ni <- 35.23; CT_min_acc17_ni <- 7.43
CT_max_acc17_i  <- 33.52; CT_min_acc17_i  <- 7.43
CT_max_acc25_ni <- 36.60; CT_min_acc25_ni <- 7.67
CT_max_acc25_i  <- 33.13; CT_min_acc25_i  <- 7.67

# Q10 parameter adapted
# Infected ACC17 2.29345144
# Infected ACC25 4.6404084
# Non-infected ACC17 2.39369045
# Non-infected ACC25 5.6438348

# --- Calculating Metabolic rate parameter 3 (M_3) based on Q10
# the formula is M3 = log10(Q10) / 10

# Parameter M_3
M_3_acc17_ni <- 0.03790680
M_3_acc17_i  <- 0.03604895
M_3_acc25_ni <- 0.06665562
M_3_acc25_i  <- 0.07515743

# Run the ectotherm model
library(raster)
library(pbapply)
library(parallel)
library(NicheMapR)
library(terra)
library(ranger)
library(dplyr)

micro <- micro_list

ecto_list <- lapply(seq_along(micro_list), function(i) {
  micro <<- micro_list[[i]]
  message(paste("Running coordinate", i))
  
  result_acc17_ni <- tryCatch({
    ectotherm(Ww_g = Ww_g_ni, alpha_max = alpha_max, alpha_min = alpha_min, M_3 = M_3_acc17_ni, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max_acc17_ni, CT_min = CT_min_acc17_ni, T_pref = T_pref_ni, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, nocturn = nocturn, diurn = diurn, crepus = crepus, nyears = 1)
  }, error = function(e) {
    message(sprintf("Error at coordinate %d (acc17_ni): %s", i, e$message))
    return(NULL)
  })
  
  result_acc25_ni <- tryCatch({
    ectotherm(Ww_g = Ww_g_ni, alpha_max = alpha_max, alpha_min = alpha_min, M_3 = M_3_acc25_ni, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max_acc25_ni, CT_min = CT_min_acc25_ni, T_pref = T_pref_ni, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, nocturn = nocturn, diurn = diurn, crepus = crepus)
  }, error = function(e) {
    message(sprintf("Error at coordinate %d (acc25_ni): %s", i, e$message))
    return(NULL)
  })
  
  result_acc17_i <- tryCatch({
    ectotherm(Ww_g = Ww_g_i, alpha_max = alpha_max, alpha_min = alpha_min, M_3 = M_3_acc17_i, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max_acc17_i, CT_min = CT_min_acc17_i, T_pref = T_pref_i, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, nocturn = nocturn, diurn = diurn, crepus = crepus)
  }, error = function(e) {
    message(sprintf("Error at coordinate %d (acc17_i): %s", i, e$message))
    return(NULL)
  })
  
  result_acc25_i <- tryCatch({
    ectotherm(Ww_g = Ww_g_i, alpha_max = alpha_max, alpha_min = alpha_min, M_3 = M_3_acc25_i, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max_acc25_i, CT_min = CT_min_acc25_i, T_pref = T_pref_i, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, nocturn = nocturn, diurn = diurn, crepus = crepus)
  }, error = function(e) {
    message(sprintf("Error at coordinate %d (acc25_i): %s", i, e$message))
    return(NULL)
  })
  
  return(list(
    acc17_ni = result_acc17_ni,
    acc25_ni = result_acc25_ni,
    acc17_i  = result_acc17_i,
    acc25_i  = result_acc25_i
  ))
})

# Extract and aggregate data frames
extract_aggregate <- function(ecto_list, scenario, varname) {
  dfs <- lapply(seq_along(ecto_list), function(i) {
    res <- ecto_list[[i]][[scenario]]
    if (!is.null(res)) {
      df <- as.data.frame(res[[varname]])
      df$point <- i
      return(df)
    } else {
      return(NULL)
    }
  })
  combined_df <- do.call(rbind, dfs)
  if (!is.null(combined_df)) {
    combined_df$scenario <- scenario
  }
  return(combined_df)
}

scenarios <- c("acc17_ni", "acc25_ni", "acc17_i", "acc25_i")

enbal_all <- do.call(rbind, lapply(scenarios, function(s) extract_aggregate(ecto_list, s, "enbal")))
environ_all <- do.call(rbind, lapply(scenarios, function(s) extract_aggregate(ecto_list, s, "environ")))

setwd("G:/My Drive/UMississippi/amphibian_data/00_var_30_nichemapr_may29")

write.csv(enbal_all, "enbal_nodeb_may29.csv", row.names = FALSE)
write.csv(environ_all, "environ_nodeb_may29.csv", row.names = FALSE)

#################

library(terra)
library(dplyr)
library(ranger)
setwd("G:/My Drive/UMississippi/amphibian_data")

enbal_df <- read.csv("enbal_nodeb_may29.csv", header = TRUE)
environ_df <- read.csv("environ_nodeb_may29.csv", header = TRUE)

enbal_df <- left_join(enbal_df, loc2, by = "point")
environ_df <- left_join(environ_df, loc2, by = "point")

enbal_df <- enbal_df %>%
  mutate(
    acc = ifelse(grepl("acc17", scenario), "ACC17", "ACC25"),
    status = ifelse(grepl("_ni", scenario), "Non-infected", "Infected")
  )

environ_df <- environ_df %>%
  mutate(
    acc = ifelse(grepl("acc17", scenario), "ACC17", "ACC25"),
    status = ifelse(grepl("_ni", scenario), "Non-infected", "Infected")
  )

### Climate 

make_cov_df <- function(r1, r2) {
  covs <- c(r1, r2)
  names(covs) <- c("bio1", "bio12")
  valid_cells <- which(!is.na(values(r1)))
  cov_values <- terra::extract(covs, valid_cells)
  coords <- terra::xyFromCell(r1, valid_cells) %>% as.data.frame()
  colnames(coords) <- c("lon", "lat")
  cbind(coords, cov_values)
}

# Present
setwd("G:/My Drive/UMississippi/amphibian_data/03_var_30/acc17_infected")
bio1 <- rast("raster_br_res30_wc2.0_bio_30s_01.tif")
bio12 <- rast("raster_br_res30_wc2.0_bio_30s_12.tif")
covariates <- c(bio1, bio12); names(covariates) <- c("bio1", "bio12")
panama_df <- make_cov_df(bio1, bio12)

# RCP 2.45
setwd("G:/My Drive/UMississippi/amphibian_data/03_var_30/mir_245_acc17")
bio1_245 <- rast("raster_br_res30_wc2.0_bio_30s_01.tif")
bio12_245 <- rast("raster_br_res30_wc2.0_bio_30s_12.tif")
covariates_245 <- c(bio1_245, bio12_245); names(covariates_245) <- c("bio1", "bio12")
panama_df_245 <- make_cov_df(bio1_245, bio12_245)

# RCP 5.85
setwd("G:/My Drive/UMississippi/amphibian_data/03_var_30/mir_585_acc17")
bio1_585 <- rast("raster_br_res30_wc2.0_bio_30s_01.tif")
bio12_585 <- rast("raster_br_res30_wc2.0_bio_30s_12.tif")
covariates_585 <- c(bio1_585, bio12_585); names(covariates_585) <- c("bio1", "bio12")
panama_df_585 <- make_cov_df(bio1_585, bio12_585) 

pts_enbal <- terra::vect(enbal_df, geom = c("lon", "lat"), crs = "EPSG:4326")
pts_env   <- terra::vect(environ_df, geom = c("lon", "lat"), crs = "EPSG:4326")

clima_para_enbal <- terra::extract(covariates, pts_enbal)
clima_para_env <- terra::extract(covariates, pts_env)

train_master_enb <- cbind(enbal_df, clima_para_enbal[,-1]) %>% 
  filter(!is.na(ENB), !is.na(bio1)) %>%
  mutate(acc = as.factor(acc), status = as.factor(status))

lookup_ct <- data.frame(
  scenario = c("acc17_ni", "acc17_i", "acc25_ni", "acc25_i"),
  CT_max = c(CT_max_acc17_ni, CT_max_acc17_i, CT_max_acc25_ni, CT_max_acc25_i)
)

train_master_env <- cbind(environ_df, clima_para_env[,-1]) %>% 
  left_join(lookup_ct, by = "scenario") %>%
  mutate(wt = CT_max - TC, tsm = CT_max - TA, 
         acc = as.factor(acc), status = as.factor(status)) %>%
  filter(!is.na(wt), !is.na(bio1))

rf_enb  <- ranger(ENB ~ bio1 + bio12 + acc + status, data = train_master_enb, num.trees = 500)
rf_wt <- ranger(wt ~ bio1 + bio12 + acc + status, data = train_master_env, num.trees = 500)
rf_tsm <- ranger(tsm ~ bio1 + bio12 + acc + status, data = train_master_env, num.trees = 500)

clim_scenarios <- list("Present" = covariates, "245" = covariates_245, "585" = covariates_585)

setwd("G:/My Drive/UMississippi/amphibian_data/00_var_30_nichemapr_may29")

for (name in names(clim_scenarios)) {
  base_df <- as.data.frame(clim_scenarios[[name]], xy = TRUE)
  colnames(base_df) <- c("lon", "lat", "bio1", "bio12")
  
  for (a in c("ACC17", "ACC25")) {
    for (s in c("Infected", "Non-infected")) {
      pred_df <- base_df %>% mutate(acc = as.factor(a), status = as.factor(s))
      
      pred_df$enb   <- predict(rf_enb, data = pred_df)$predictions
      pred_df$wt <- predict(rf_wt, data = pred_df)$predictions
      pred_df$tsm <- predict(rf_tsm, data = pred_df)$predictions
      
      salvar_tif <- function(valores, var_nome) {
        r <- rast(clim_scenarios[[name]][[1]])
        values(r) <- NA
        r[cellFromXY(r, pred_df[,1:2])] <- valores
        writeRaster(r, paste0(var_nome, "_", a, "_", s, "_", name, ".tif"), overwrite = TRUE)
      }
      
      salvar_tif(pred_df$enb, "ENB")
      salvar_tif(pred_df$wt, "wt")
      salvar_tif(pred_df$tsm, "tsm")
      message("Generated: ", a, " | ", s, " | ", name)
    }
  }
}

####################
## PLOTS ##

out_dir <- "G:/My Drive/UMississippi/amphibian_data/00_var_30_nichemapr_may29"

if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

arquivos_finais <- list.files(out_dir, pattern = ".tif$", full.names = TRUE)

lista_df <- lapply(arquivos_finais, function(f) {
  nome <- basename(f)
  partes <- unlist(strsplit(gsub(".tif", "", nome), "_"))
  r <- terra::rast(f)
  df <- as.data.frame(r, xy = FALSE)
  colnames(df)[1] <- "valor"
  
  df <- df %>% dplyr::mutate(
    var     = partes[1],
    acc     = partes[2],
    status  = partes[3],
    cenario = partes[4]
  )
  return(df)
})

df_final <- dplyr::bind_rows(lista_df)

table(df_final$status, df_final$cenario, df_final$var)

ordem_cenarios <- c("Present", "245", "585")

df_summary_enb <- df_final %>%
  filter(var == "ENB") %>%
  group_by(status, cenario, acc) %>%
  summarise(media = mean(valor), sd = sd(valor), .groups = 'drop') %>%
  mutate(cenario = factor(cenario, levels = ordem_cenarios))

df_summary_wt <- df_final %>%
  filter(var == "wt") %>%
  group_by(status, cenario, acc) %>%
  summarise(media = mean(valor), sd = sd(valor), .groups = 'drop') %>%
  mutate(cenario = factor(cenario, levels = ordem_cenarios))

df_summary_tsm <- df_final %>%
  filter(var == "tsm") %>%
  group_by(status, cenario, acc) %>%
  summarise(media = mean(valor), sd = sd(valor), .groups = 'drop') %>%
  mutate(cenario = factor(cenario, levels = ordem_cenarios))

library(ggplot2)

plot_enb <- ggplot(df_summary_enb, aes(x = cenario, y = media, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~acc, scales = "free") +
  theme_classic() +
  labs(title = "Mean Energy Balance (ENB)", y = "J/day", x = "Scenario")

plot_wt <- ggplot(df_summary_wt, aes(x = cenario, y = media, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~acc, scales = "free") +
  theme_classic() +
  labs(title = "Warming Tolerance (WT)", y = "Index", x = "Scenario")

plot_tsm <- ggplot(df_summary_tsm, aes(x = cenario, y = media, fill = status)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), 
                width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~acc, scales = "free") +
  theme_classic() +
  labs(title = "Thermal Safety Margin (TSM)", y = "Index", x = "Scenario")

pasta_graficos <- file.path(out_dir, "09_stats_30")
if(!dir.exists(pasta_graficos)) dir.create(pasta_graficos, recursive = TRUE)

ggsave(file.path(pasta_graficos, "plot_wt.png"), plot = plot_wt, 
       width = 10, height = 6, dpi = 300)

ggsave(file.path(pasta_graficos, "plot_tsm.png"), plot = plot_tsm, 
       width = 10, height = 6, dpi = 300)

ggsave(file.path(pasta_graficos, "plot_enb.png"), plot = plot_enb, 
       width = 10, height = 6, dpi = 300)

