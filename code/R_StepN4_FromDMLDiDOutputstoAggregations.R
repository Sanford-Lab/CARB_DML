library(tidyverse);library(splines);library(tictoc);library(reshape2)
library(DoubleML);library(mlr3);library(mlr3learners);library(mlr3tuning)
library(sf);library(paletteer);library(lubridate)
library(ranger);library(scales);library(DRDID);library(Matrix)

# Project-specific information =======

### Project area (in ha) =======

arboc.map.agg <- terra::vect("data/input/OffsetsMap_AGGREGATED.shp")
arboc.map.agg <- terra::project(arboc.map.agg, "epsg:4326")
arboc.map.agg$area <- terra::expanse(arboc.map.agg, unit = "ha")
arboc.map.agg.df <- as.data.frame(arboc.map.agg[,c("ARB_Projec", "area")])

arboc.map.agg.df %>%
  rename(project.name = ARB_Projec) -> arboc.map.agg.df

arboc.map.agg.df$project.name[which(arboc.map.agg.df$project.name=="CAFR5213\r\n")] <- "CAFR5213"

save(list = "arboc.map.agg.df",
     file = "data/output/StepN4_ProjectInfo_Area.RData")

### Other project characteristics =======

mode.giver <- function(x) {
  # Remove NA values from the input vector
  x_clean <- x[!is.na(x)]
  
  # If the vector is empty after removing NAs (i.e., it only contained NAs),
  # return NA of type character to ensure consistency in summarise.
  if (length(x_clean) == 0) {
    return(NA_character_)
  }
  
  # Create a frequency table of the non-NA values
  t <- table(x_clean)
  
  # Return the name of the most frequent value.
  # In case of a tie, which.max() returns the index of the first maximum value.
  names(t)[which.max(t)]
}


fownership.keys <- unlist(list("1" = "Family (Private)",
                               "2" = "Corporate (Private)",
                               "3" = "TIMO/REIT (Private)",
                               "4" = "Other Private",
                               "5" = "Federal (Public)",
                               "6" = "State (Public)",
                               "7" = "Local (Public)",
                               "8" = "Tribal"))

project_characteristics_list <- list()

for (chunkno in 7:7) {
  load(paste0("data/output/StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk", chunkno, ".RData"))
  #the name of this variable will be paste0("carb_data_chunk_", chunkno)
  
  this_chunk_carb_data <- get(paste0("carb_data_chunk_", chunkno))
  
  rm(list = paste0("carb_data_chunk_", chunkno))
  
  print(paste0("DONE BRINGING IN DATA"))
  
  projects_in_this_chunk <- unique(this_chunk_carb_data$projectID)
  
  this_chunk_carb_data %>% 
    filter(Year == START.YEAR-1) %>% 
    #use baseline data
    group_by(projectID, treat) %>% 
    summarise(fownership_mode = mode.giver(fownership),
              average_biomass = mean(biomass, na.rm =T),
              sd_biomass = sd(biomass, na.rm = T),
              nlcd_mode = mode.giver(nlcd),
              distance.to.road.mean = mean(distance.to.road, na.rm = T)
              ) -> this_chunk_carb_data_summarized
  
  this_chunk_carb_data_summarized$fownership_mode_str <- unname(fownership.keys[match(as.character(this_chunk_carb_data_summarized$fownership_mode),
                                                                                           names(fownership.keys))])
  
  this_chunk_carb_data_summarized %>% 
    mutate(priv.type.public = case_when(
      fownership_mode_str == "Corporate (Private)" ~ "Corporate (Private)",
      fownership_mode_str == "Family (Private)" ~ "Family (Private)",
      fownership_mode_str %in% c("TIMO/REIT (Private)", "Other Private") ~ "TIMO/Other private",
      fownership_mode_str %in% c("State (Public)", "Tribal", "Federal (Public)") ~ "Public/Tribal"
    ), priv.or.public = case_when(
      grepl("Private", fownership_mode_str) ~ "Private",
      TRUE ~ "Public"
    )
    ) -> this_chunk_carb_data_summarized
  
  this_chunk_carb_data_summarized %>% 
    mutate(treat = case_when(
      treat == 0 ~ "control",
      treat == 1 ~ "treatment"
    )) %>% 
    pivot_wider(id_cols = "projectID", names_from = "treat", 
                values_from = contains(c("fownership", "biomass", "nlcd", "distance"))) -> this_chunk_carb_data_summarized_wider
  
  
  project_characteristics_list[[as.character(chunkno)]] <- this_chunk_carb_data_summarized_wider
  
  print(paste0("DONE WITH CHARACTERISTICS DATA FOR CHUNK NO. ", chunkno))
  
  rm(this_chunk_carb_data)
}

project_characteristics <- bind_rows(project_characteristics_list)

project_characteristics %>% 
  select(!contains("_NA")) -> project_characteristics

project_type_df <- readxl::read_excel("data/input/OFFSET_ISSUANCE DATA_ProjectType.xlsx")

project_characteristics %>% 
  left_join(project_type_df, by = c("projectID" = "proj.name")) -> project_characteristics

project_characteristics %>% 
  left_join(arboc.map.agg.df %>% rename(area_inha = area),
            by = c("projectID" = "project.name")) -> project_characteristics

save(list = "project_characteristics",
     file = "data/output/StepN4_ProjectInfo_AreawithOtherInfo.RData")

# Going from chunks 1 to 3 ======

load("data/output/StepN4_ProjectInfo_AreawithOtherInfo.RData")

chunk_largeprojs_combinedwithoutput_list <- list()

for (chunkno in 1:3) {
  load(paste0("data/output/StepN3_Showcase_CARB0030andOtherChunk", chunkno, "_Outcome.RData"))
  
  this_chunk_results <- get(paste0("chunk", chunkno, "_results_list"))
  
  rm(list = paste0("chunk", chunkno, "_results_list"))
  
  load(paste0("data/output/StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk", chunkno, ".RData"))
  #the name of this variable will be paste0("carb_data_chunk_", chunkno)
  
  this_chunk_carb_data <- get(paste0("carb_data_chunk_", chunkno))
  
  rm(list = paste0("carb_data_chunk_", chunkno))
  
  print(paste0("DONE BRINGING IN DATA"))
  
  this_chunk_results_combined <- bind_rows(this_chunk_results$MLBased %>% 
                                         map("pscore_and_ell_withIFF") %>%
                                         bind_rows())
  
  this_chunk_largeprojs_combinedwithoutput <- this_chunk_carb_data %>% 
    filter(projectID %in% subset(project_characteristics, area_inha > 1000)$projectID) %>% 
    #only focus on those with sufficiently large treated pixels
    mutate(year.to.treat = as.numeric(Year) - year(DATE.first)) %>% 
    select(-treat) %>%  
    left_join(this_chunk_results_combined %>% 
                filter(projectID %in% subset(project_characteristics, area_inha > 1000)$projectID) %>% 
                select(projectID, cellID, year_val, fold,ghat, phat, treat, outcome.differenced, ellhat, psi1.pre.this.fold, iff.this.fold, att),
              by = c("projectID", "cellID", "year.to.treat" = "year_val")) %>%
    mutate(Year = as.numeric(Year),
           projectID_cellID = dense_rank(paste0(projectID, "-", cellID)),
           projectID_cluster = dense_rank(paste0(projectID, "-", cluster.25km)),
           cellID = as.numeric(cellID)) %>% 
    filter(!is.na(psi1.pre.this.fold) & !is.na(att))
  #drop all cases where the propensity score-trimmed pixels are dropped
  
  rm(this_chunk_carb_data)
  
  chunk_largeprojs_combinedwithoutput_list[[as.character(chunkno)]] <- this_chunk_largeprojs_combinedwithoutput
  
  rm(this_chunk_largeprojs_combinedwithoutput)
  
  print(paste0("DONE WITH CHUNK NO. ", chunkno))
}

chunk_largeprojs_combinedwithoutput <- bind_rows(chunk_largeprojs_combinedwithoutput_list)

manual_did_inference <- function(dml_data, 
                                 cluster_var = "cluster.25km", 
                                 time_var = "year.to.treat", 
                                 alpha = 0.05, 
                                 biter = 1000) {
  
  # 1. SETUP: Unique Cluster IDs
  if("projectID" %in% names(dml_data)){
    dml_data <- dml_data %>%
      mutate(unique_cluster_id = paste0(projectID, "_", .data[[cluster_var]]))
  } else {
    dml_data <- dml_data %>%
      mutate(unique_cluster_id = as.character(.data[[cluster_var]]))
  }
  
  # 2. POINT ESTIMATES & COUNTS (N_obs)
  results_table <- dml_data %>%
    group_by(!!sym(time_var)) %>%
    summarise(
      att = mean(att, na.rm = TRUE),
      n_obs = n(),           # Count of observations (units)
      n_clusters_t = n_distinct(unique_cluster_id), # Count of clusters
      .groups = 'drop'
    ) %>%
    arrange(!!sym(time_var))
  
  # 3. AGGREGATE TO CLUSTER LEVEL (Sum IFFs)
  cluster_iff <- dml_data %>%
    group_by(unique_cluster_id, !!sym(time_var)) %>%
    summarise(psi = sum(iff.this.fold, na.rm = TRUE), .groups = 'drop')
  
  # 4. RESHAPE TO WIDE
  psi_matrix_df <- cluster_iff %>%
    pivot_wider(names_from = !!sym(time_var), 
                values_from = psi, 
                values_fill = 0) 
  
  psi_matrix <- as.matrix(psi_matrix_df %>% select(-unique_cluster_id))
  
  # Match column order
  ordered_times <- as.character(results_table[[time_var]])
  psi_matrix <- psi_matrix[, ordered_times]
  
  # 5. CALCULATE STANDARD ERRORS
  # Formula: SE = sqrt( Sum(Psi_g^2) * (G / G-1) ) / N_obs
  
  sum_sq_psi <- colSums(psi_matrix^2)
  n_clusters <- nrow(psi_matrix)
  
  # Finite Sample Adjustment (G / G-1)
  fpc <- n_clusters / (n_clusters - 1)
  
  # Divide by N_obs, not N_clusters
  se_vec <- sqrt(sum_sq_psi * fpc) / results_table$n_obs 
  
  results_table$se <- se_vec
  
  # 6. SIMULTANEOUS CONFIDENCE BANDS (Multiplier Bootstrap)
  weights <- matrix(sample(c(-1, 1), n_clusters * biter, replace = TRUE), 
                    nrow = n_clusters, ncol = biter)
  
  # Bootstrap Means: (Sum(Psi * w) / N_obs)
  boot_sums <- t(psi_matrix) %*% weights
  boot_means <- boot_sums / results_table$n_obs
  
  # t-statistics
  se_matrix <- matrix(se_vec, nrow = nrow(boot_means), ncol = ncol(boot_means), byrow = FALSE)
  t_stat_matrix <- abs(boot_means) / se_matrix
  
  # Handle 0/0 for empty years
  t_stat_matrix[is.nan(t_stat_matrix) | is.infinite(t_stat_matrix)] <- 0
  
  t_stats <- apply(t_stat_matrix, 2, max, na.rm = TRUE)
  crit_val <- quantile(t_stats, 1 - alpha, na.rm = TRUE)
  
  # 7. FINALIZE
  results_table <- results_table %>%
    mutate(
      c = crit_val,
      lower_pointwise = att - qnorm(1 - alpha/2) * se,
      upper_pointwise = att + qnorm(1 - alpha/2) * se,
      lower_simul     = att - crit_val * se,
      upper_simul     = att + crit_val * se
    )
  
  return(results_table)
}

table(project_characteristics$fownership_mode_str_treatment)

upto3_attgt <- manual_did_inference(chunk_largeprojs_combinedwithoutput)
upto3_family_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                           projectID %in% subset(project_characteristics,
                                                                 fownership_mode_str_treatment %in% "Family (Private)")$projectID))
upto3_corporate_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                                  projectID %in% subset(project_characteristics,
                                                                        fownership_mode_str_treatment %in% "Corporate (Private)")$projectID))

upto3_smallsize_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                                     projectID %in% subset(project_characteristics,
                                                                           area_inha < median(area_inha, na.rm = T))$projectID))

upto3_lowbiomass_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                                     projectID %in% subset(project_characteristics,
                                                                           average_biomass_treatment < median(average_biomass_treatment, na.rm = T))$projectID))

upto3_highbiomass_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                                      projectID %in% subset(project_characteristics,
                                                                            average_biomass_treatment >= median(average_biomass_treatment, na.rm = T))$projectID))

upto3_ifm_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                                     projectID %in% subset(project_characteristics,
                                                                           type %in% "IFM")$projectID))

upto3_ac_attgt <- manual_did_inference(subset(chunk_largeprojs_combinedwithoutput, 
                                               projectID %in% subset(project_characteristics,
                                                                     !type %in% "IFM")$projectID))

ggplot(data = subset(upto3_attgt, abs(year.to.treat) <= 10)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw()

ggplot(data = subset(upto3_family_attgt, abs(year.to.treat) <= 10)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw()


ggplot(data = subset(upto3_corporate_attgt, abs(year.to.treat) <= 10)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw()

ggplot(data = subset(upto3_ac_attgt, abs(year.to.treat) <= 10)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw()

ggplot() +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'low biomass/ha'),
                  data = subset(upto3_lowbiomass_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = -0.2)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'high biomass/ha'),
                  data = subset(upto3_highbiomass_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = 0.2)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw() +
  theme(legend.position = 'bottom')

ggplot() +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'Family-owned'),
                  data = subset(upto3_family_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = -0.2)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'Corporate-owned'),
                  data = subset(upto3_corporate_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = 0.2)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw() +
  theme(legend.position = 'bottom')

ggplot() +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'IFM projects'),
                  data = subset(upto3_ifm_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = -0.2)) +
  geom_pointrange(aes(x = year.to.treat, y = att, ymax = upper_simul,
                      ymin = lower_simul, color = 'Avoided conversion projects'),
                  data = subset(upto3_ac_attgt, abs(year.to.treat) <= 10),
                  position = position_nudge(x = 0.2)) +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw() +
  theme(legend.position = 'bottom')



aggte(first_five_projs_attgt)
