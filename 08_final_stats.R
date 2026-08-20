# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 08. Final Stats

rm(list = ls())
library(terra)
library(sf)
library(tidyverse)
library(nlme)
library(multcomp)
library(tidytext)

path <- "G:/My Drive/UMississippi/amphibian_data"
path_alg <- file.path(path, "04_algorithms_30")
path_ens <- file.path(path, "05_ensemble_30")
path_stats <- file.path(path, "09_stats_30")
path_bin <- file.path(path, "06_binary_30")
dir.create(path_stats, showWarnings = FALSE, recursive = TRUE)

#  1. AUC FILTER 
eval_files <- list.files(path_alg, pattern = "^evaluation_.*\\.csv$", full.names = TRUE, recursive = TRUE)
df_eval <- eval_files %>%
  map_df(~read_csv(.x, show_col_types = FALSE) %>% 
           mutate(scenario = str_extract(.x, "(correlative_acc\\d+|acc\\d+_(infected|noninfected))")))
valid_reps <- df_eval %>% filter(auc > 0.75) %>% dplyr::select(scenario, species, algorithm, replica)

#  2. IMPORTANCE 
message(">>> Generating Variable Importance...")
imp_files <- list.files(path_alg, pattern = "importance", full.names = TRUE, recursive = TRUE)

df_imp <- imp_files %>%
  map_df(function(f) {
    d <- read_csv(f, show_col_types = FALSE)
    scen <- str_extract(f, "(correlative_acc\\d+|acc\\d+_(infected|noninfected))")
    d$scenario <- scen
    d$model_group <- case_when(
      str_detect(scen, "correlative") ~ "Baseline",
      str_detect(scen, "noninfected") ~ "Non-infected",
      TRUE ~ "Infected"
    )
    return(d)
  }) %>%
  inner_join(valid_reps, by = c("scenario", "species", "algorithm", "replica"))

# DEFINING THE ORDER OF THE 6 PANELS (17 ON TOP, 25 ON BOTTOM)
ordem_paineis <- c(
  "correlative_acc17", "acc17_noninfected", "acc17_infected", # Top row
  "correlative_acc25", "acc25_noninfected", "acc25_infected"  # Bottom row
)

df_imp$scenario <- factor(df_imp$scenario, levels = ordem_paineis)
df_imp$model_group <- factor(df_imp$model_group, levels = c("Baseline", "Non-infected", "Infected"))

# Function to clean variable names
limpar_nomes <- function(x) {
  x <- gsub("raster_br_res30_wc2.0_bio_30s_", "Bio", x) # Bio01, Bio02...
  x <- gsub("wt_ACC17_Non.infected", "WT", x)
  x <- gsub("tsm_ACC17_Non.infected", "TSM", x)
  x <- gsub("wt_ACC17_Infected", "WT", x)
  x <- gsub("tsm_ACC17_Infected", "TSM", x)
  x <- gsub("wt_ACC25_Non.infected", "WT", x)
  x <- gsub("tsm_ACC25_Non.infected", "TSM", x)
  x <- gsub("wt_ACC25_Infected", "WT", x)
  x <- gsub("tsm_ACC25_Infected", "TSM", x)
  x <- gsub("ENB_ACC17_Non.infected", "ENB", x)
  x <- gsub("ENB_ACC17_Infected", "ENB", x)
  x <- gsub("ENB_ACC25_Non.infected", "ENB", x)
  x <- gsub("ENB_ACC25_Infected", "ENB", x)
  x <- gsub("ENB_Present", "ENB", x)
  
  return(x)
}

# Apply cleaning
imp_summary <- df_imp %>%
  mutate(variable = limpar_nomes(as.character(variable))) %>%
  group_by(scenario, model_group, variable) %>%
  summarise(mean_imp = mean(importance, na.rm = TRUE), .groups = "drop") %>%
  mutate(variable = tidytext::reorder_within(variable, mean_imp, scenario))

fig_imp <- ggplot(imp_summary, aes(x = variable, y = mean_imp, fill = model_group)) +
  geom_bar(stat = "identity") +
  tidytext::scale_x_reordered() + 
  coord_flip() +
  facet_wrap(~scenario, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Baseline" = "#999999", "Non-infected" = "#1f78b4", "Infected" = "#d95f02")) +
  theme_bw(base_size = 14) + # Increase base font size to 14
  theme(
    axis.text = element_text(size = 12),       # Increase axis text size
    axis.title = element_text(size = 14, face = "bold"), # Increase axis title size
    strip.text = element_text(size = 12, face = "bold"), # Increase facet title size (scenarios)
    legend.text = element_text(size = 12),     # Increase legend text size
    legend.title = element_text(size = 13)     # Increase legend title size
  ) +
  labs(x = "Predictors", y = "Importance (Mean AUC loss)", fill = "Model Type")

ggsave(file.path(path_stats, "Variable_Importance.png"), fig_imp, width = 16, height = 12)

#  3. FOREST PLOT 
message(">>> Generating Forest Plot...")
ens_files <- list.files(path_ens, pattern = "\\.tif$", full.names = TRUE)
all_res <- list()

for(ac in c("acc17", "acc25")){
  for(sce in c("present", "245", "585")){
    p_files <- ens_files[grepl(paste0(ac, ".*", sce), ens_files)]
    if(length(p_files) < 2) next
    s <- terra::rast(p_files)
    names(s) <- case_when(str_detect(names(s), "noninfected") ~ "Mec_NonInf", 
                          str_detect(names(s), "infected") ~ "Mec_Inf", TRUE ~ "Baseline")
    
    v <- terra::spatSample(s, 10000, na.rm = T) %>% mutate(p_id = row_number()) %>%
      pivot_longer(-p_id, names_to = "type", values_to = "val")
    v$type <- factor(v$type, levels = c("Baseline", "Mec_NonInf", "Mec_Inf"))
    
    mod <- lme(val ~ type, random = ~1|p_id, data = v, control = lmeControl(opt = "optim"))
    d_res <- summary(glht(mod, linfct = mcp(type = "Dunnett")))
    ci <- confint(d_res)$confint
    
    # SIGN INVERSION: If Mech < Base, the result should be negative
    all_res[[paste(ac, sce)]] <- data.frame(
      Comp = names(d_res$test$coefficients), 
      Est = as.numeric(d_res$test$coefficients), 
      L = as.numeric(ci[,"lwr"]), 
      U = as.numeric(ci[,"upr"]), 
      Acc = ac, Sce = sce)
  }
}

df_forest <- bind_rows(all_res) %>%
  mutate(Status = factor(ifelse(str_detect(Comp, "Non"), "Non-Infected", "Infected"), levels = c("Non-Infected", "Infected")),
         Season = factor(ifelse(Acc == "acc17", "Wet Season (17°C)", "Dry Season (25°C)"), levels = c("Wet Season (17°C)", "Dry Season (25°C)")),
         Sce = factor(case_when(Sce=="present"~"Present", Sce=="245"~"SSP245", T~"SSP585"), levels = c("Present", "SSP245", "SSP585")))

fig_forest <- ggplot(df_forest, aes(x = Est, y = Sce, color = Status)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = U, xmax = L), height = 0.3, linewidth = 1, position = position_dodge(0.6)) +
  geom_point(size = 4, position = position_dodge(0.6)) + facet_grid(Season ~ .) +
  scale_color_manual(values = c("Non-Infected" = "#1f78b4", "Infected" = "#d95f02")) + 
  theme_bw() + labs(x = "Effect Size (Mechanistic - Baseline)")

ggsave(file.path(path_stats, "ForestPlot.png"), fig_forest, width = 11, height = 7)


#  4. PROTECTED AREAS

message(">>> Generating PAs Analysis...")

# 1. Load Protected Areas (PAs) and Binary Rasters
ucs_sf <- st_read(file.path(path, "07_maps_30/pa_allpoints.gpkg"))
bin_files <- list.files(path_bin, pattern = "MAX_TSS\\.tif$", full.names = TRUE)

# 2. Extract total valid pixels per PA (Denominator)
r_ref <- terra::rast(bin_files[1])
ucs_vect <- terra::vect(ucs_sf)
total_px <- terra::extract(r_ref, ucs_vect, fun = function(x) sum(!is.na(x)))

# 3. Extract suitability (Pixels = 1) for each scenario
results_pa_list <- list()
for(f in bin_files) {
  r <- terra::rast(f)
  vals <- terra::extract(r, ucs_vect, fun = sum, na.rm = TRUE)
  results_pa_list[[basename(f)]] <- vals[, 2]
}

# 4. Organize Dataframe and Calculate Percentage
df_pa_wide <- as.data.frame(results_pa_list)
df_pa_wide <- bind_cols(ucs_sf %>% st_drop_geometry() %>% dplyr::select(NOMBRE), df_pa_wide)

df_slope <- df_pa_wide %>%
  mutate(ID = row_number(), total_pixels_uc = total_px[, 2]) %>%
  pivot_longer(cols = starts_with("BIN_"), names_to = "layer", values_to = "suit_pixels") %>%
  mutate(
    perc = (suit_pixels / total_pixels_uc) * 100,
    # ROWS: Wet season (17°C) on top, Dry season (25°C) on bottom
    acclimation = factor(case_when(
      str_detect(layer, "acc17") ~ "Wet season (17°C)",
      str_detect(layer, "acc25") ~ "Dry season (25°C)"
    ), levels = c("Wet season (17°C)", "Dry season (25°C)")),
    # COLUMNS: Baseline -> Non-infected -> Infected
    model = factor(case_when(
      str_detect(layer, "correlative") ~ "Baseline", 
      str_detect(layer, "noninfected") ~ "Non-infected", 
      TRUE ~ "Infected"
    ), levels = c("Baseline", "Non-infected", "Infected")),
    # X-AXIS: Present -> SSP245 -> SSP585
    scenario = factor(case_when(
      str_detect(layer, "present") ~ "Present", 
      str_detect(layer, "245") ~ "SSP245", 
      TRUE ~ "SSP585"
    ), levels = c("Present", "SSP245", "SSP585"))
  )

# 5. Save Raw Data
write_csv(df_slope, file.path(path_stats, "SlopePlot_data.csv"))

# 6. Generate and Save Plot
fig_slope <- ggplot(df_slope, aes(x = scenario, y = perc, group = NOMBRE, color = NOMBRE)) +
  geom_line(linewidth = 1, alpha = 0.7) + 
  geom_point(size = 2) + 
  facet_grid(acclimation ~ model) + # 17 on top | Baseline on the left
  scale_color_viridis_d() + 
  theme_bw() + 
  labs(title = "Protected Areas Suitability Trends",
       x = "Climate Scenario", 
       y = "Suitable Habitat Area (%)",
       color = "Protected Area") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(path_stats, "SlopePlot.png"), fig_slope, width = 12, height = 8, dpi = 300)


# SUMMARY TABLES (CSV)

message(">>> Saving Summary Tables...")

# 1. VARIABLE IMPORTANCE TABLE (Mean per scenario)
tab_importancia <- df_imp %>%
  group_by(scenario, model_group, variable) %>%
  summarise(
    Mean_Importance = mean(importance, na.rm = TRUE),
    SD_Importance = sd(importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, desc(Mean_Importance))

write_csv(tab_importancia, file.path(path_stats, "Table_Variable_Importance_Full.csv"))

# 2. FOREST PLOT TABLE (LMM Results)
tab_forest <- df_forest %>%
  dplyr::select(Season, Sce, Status, Effect_Size = Est, Lower_CI = L, Upper_CI = U) %>%
  mutate(across(where(is.numeric), ~round(., 4)))

write_csv(tab_forest, file.path(path_stats, "Table_LMM_ForestPlot_Results.csv"))

# 3. PA SUMMARY TABLE (Mean per protection category)
tab_pa_summary <- df_slope %>%
  group_by(acclimation, model, scenario) %>%
  summarise(
    Mean_Suitability_Perc = mean(perc, na.rm = TRUE),
    Min_Suitability = min(perc, na.rm = TRUE),
    Max_Suitability = max(perc, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(tab_pa_summary, file.path(path_stats, "Table_ProtectedAreas_Summary.csv"))

message("✅ TABLES SAVED! Check folder 09_stats_30:")
message("- Table_Variable_Importance_Full.csv")
message("- Table_LMM_ForestPlot_Results.csv")
message("- Table_ProtectedAreas_Summary.csv")