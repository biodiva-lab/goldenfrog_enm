# Disease-mediated thermal limits in frogs
# Diele-Viegas et al.
# Script 04. Evaluation and validation

rm(list = ls())
library(tidyverse)
library(wesanderson)

# Paths
path <- "G:/My Drive/UMississippi/amphibian_data"
path_input <- file.path(path, "04_algorithms_30")


# 1. AUTOMATIC IMPORT
eval_files <- list.files(path_input, 
                         pattern = "^evaluation_.*\\.csv$", 
                         full.names = TRUE, 
                         recursive = TRUE)

if(length(eval_files) == 0) stop("No 'evaluation_' file found!")

message("Files found: ", length(eval_files))

read_eval <- function(f) {
  d <- readr::read_csv(f, show_col_types = FALSE)
  
  # 1. Correct column name (the generator script saves as 'algorithm')
  if ("algorithm" %in% colnames(d)) {
    d <- dplyr::rename(d, model = algorithm)
  }
  
  # 2. Extract exact scenario name from the strict root folder name
  d$folder_scen <- stringr::str_extract(f, "(correlative_acc\\d+|acc\\d+_(infected|noninfected))")
  
  return(d)
}

# Bind everything together and remove any garbage not from target folders
eva_all <- purrr::map_df(eval_files, read_eval) %>% filter(!is.na(folder_scen))


# 2. PROCESSING LOOP

scenarios <- unique(eva_all$folder_scen)

for(scen in scenarios){
  
  message("\n>>> Processing scenario: ", scen)
  eva_scen <- eva_all %>% filter(folder_scen == scen)
  
  path_eval_base <- file.path(path_input, scen, "01_evaluation")
  path_tables    <- file.path(path_eval_base, "01_tables")
  path_plots     <- file.path(path_eval_base, "02_boxplot")
  
  dir.create(path_tables, showWarnings = FALSE, recursive = TRUE)
  dir.create(path_plots,  showWarnings = FALSE, recursive = TRUE)
  
  # --- SPECIES LOOP ---
  for(sp in unique(eva_scen$species)){
    
    message("     Generating metrics for: ", sp)
    
    # 1. SUMMARY TABLE (AUC and TSS)
    eva_table <- eva_scen %>% 
      filter(species == sp) %>% 
      group_by(folder_scen, species, model) %>% 
      summarise(across(any_of(c("auc", "tss", "tss_spec_sens")), 
                       list(mean = ~mean(.x, na.rm = TRUE), 
                            sd = ~sd(.x, na.rm = TRUE)),
                       .names = "{.col}_{.fn}"),
                .groups = "drop") %>% 
      mutate(across(where(is.numeric), ~round(.x, 3)))
    
    write_csv(eva_table, file.path(path_tables, paste0("summary_", sp, ".csv")))
    
    # 2. BOXPLOTS
    eva_long <- eva_scen %>% 
      filter(species == sp) %>% 
      select(model, auc, any_of(c("tss", "tss_spec_sens"))) %>% 
      pivot_longer(cols = -model, names_to = "metric", values_to = "value")
    
    for(m in unique(eva_long$metric)){
      
      # Define cutoff line (0.8 AUC / 0.5 TSS)
      h_line <- ifelse(str_detect(m, "auc"), 0.8, 0.5)
      
      p <- ggplot(eva_long %>% filter(metric == m)) +
        aes(x = model, y = value, fill = model) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 3, alpha = 0.5) +
        scale_fill_manual(values = wes_palette("Darjeeling1", n = length(unique(eva_scen$model)))) +
        geom_hline(yintercept = h_line, linetype = "dashed", color = "red") +
        labs(
          title = paste("Model Performance:", sp),
          subtitle = paste("Scenario:", scen, "| Metric:", toupper(m)),
          x = "Algorithm",
          y = toupper(m)
        ) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        theme_bw() +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold.italic", size = 16),
              strip.text = element_text(face = "bold"))
      
      ggsave(file.path(path_plots, paste0("boxplot_", m, "_", sp, ".png")), 
             plot = p, width = 15, height = 10, units = "cm", dpi = 300)
    }
  }
}

message("\nFinished! AUC and TSS boxplots should now be in the folders.")