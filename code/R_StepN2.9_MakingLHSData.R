library(tidyverse);library(lubridate);library(Rforestry);library(reshape2)
library(randomForest);library(tictoc);library(sf);library(paletteer)

load("Step2_1_Fixed characteristics_NDVI250m.RData")
crosswalk <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")
load("Step2_BIOMASS_NDVI250m.RData")
load("Step2_BIOMASS_NDVI250m_2018to2023.RData")
load("StepN2.5_NonForestAndTreatedPixels.RData")

biomass.extracted.all <- bind_cols(biomass.extracted.all, biomass.new.df %>% select(-cellID)) %>% relocate(cellID)

colnames(biomass.extracted.all) <- c("cellID", paste0("biomass_", 2000:2023))

biomass.extracted.all %>% 
  pivot_longer(cols = starts_with("biomass_"), names_to = "Year", values_to = "biomass",  
               names_prefix = "biomass_") -> biomass.longer

crosswalk %>% 
  filter(type %in% c("exact.treated", "control.b2k")) -> crosswalk.treated.and.control

rm(crosswalk)
rm(biomass.extracted.all)

biomass.longer %>% 
  left_join(crosswalk.treated.and.control %>% select(proj.id, cellid, type), by = c("cellID" = "cellid")) -> biomass.longer.w.projid

biomass.longer.w.projid %>% 
  arrange(proj.id, cellID, Year) -> biomass.longer.w.projid

biomass.longer.w.projid %>% 
  left_join(proj.years, by = c("proj.id" = "project.name")) -> biomass.longer.w.projid

load("StepN2.5_ProjectSpecificData_b2k_updated.RData")

biomass.longer.w.projid %>% 
  mutate(treat.year = case_when(
    type == "control.b2k" ~ 0,
    type == "exact.treated" ~ year(DATE.first)
  ),
  exit.year = case_when(
    type == "control.b2k" ~ 0,
    type == "exact.treated" ~ year(DATE.last) +1
  )) -> biomass.longer.w.projid

projs <- unique(biomass.longer.w.projid$proj.id)
#138 projects

for (i in 5:7) {
  
  if (i %in% 1:6) {
    projs.in.this.chunk <- projs[(20*(i-1)+1):(20*i)]
  } else {
    projs.in.this.chunk <- projs[121:138]
  }
  
  biomass.longer.w.projid.in.this.chunk <- subset(biomass.longer.w.projid, proj.id %in% projs.in.this.chunk)
  
  biomass.longer.w.projid.wcovs.in.this.chunk <- left_join(biomass.longer.w.projid.in.this.chunk, 
                                                           proj.dats.for.prop.b2k %>% filter(project.ID %in% projs.in.this.chunk) %>% 
                                                             select(-type) %>% 
                                                             mutate(cellID = as.numeric(cellID)), 
                                             by = c("proj.id" = "project.ID", "cellID")) 
  
  rm(biomass.longer.w.projid.in.this.chunk)
  
  biomass.longer.w.projid.wcovs.in.this.chunk %>% 
    relocate(projectID = proj.id, Year, cellID, treat.year, exit.year, biomass,
             starts_with("cluster")) -> biomass.longer.w.projid.wcovs.in.this.chunk
  
  biomass.longer.w.projid.wcovs.in.this.chunk$cellID <- as.character(biomass.longer.w.projid.wcovs.in.this.chunk$cellID)
  
  assign(paste0("carb_data_chunk_", i),
         biomass.longer.w.projid.wcovs.in.this.chunk)
  
  save(list = paste0("carb_data_chunk_", i),
       file = paste0("StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk", i, ".RData"))
  
  rm(list = paste0("carb_data_chunk_", i))
  rm(biomass.longer.w.projid.wcovs.in.this.chunk)
}

rm(list = ls())

# Practice run ====

library(did)

load("StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk1.RData")

carb_data_cafr0030 <- subset(carb_data_chunk_1, projectID == "CAFR0030")

cov.wo.clm.practicerun <- c("clm_DEM", "clm_ned_lf", 
                         "nlcd", 'fownership',
                         "distance.to.road"
                         )

did.onlytreated.chunk1 <- att_gt(yname = "biomass", tname = "Year", idname = "cellID", 
                           gname = "treat.year", 
                           data = carb_data_chunk_1 %>% filter(treat.year!=0) %>% mutate(Year = as.numeric(Year),
                                                                                         cellID = as.numeric(cellID)) %>% 
                             filter(!cellID %in% non.forest.pixels), control_group= "notyettreated", 
                           panel = T, allow_unbalanced_panel = T)

did.plainvanilla <- att_gt(yname = "biomass", tname = "Year", idname = "cellID", 
                           gname = "treat.year", 
                           data = carb_data_cafr0030 %>% mutate(Year = as.numeric(Year),
                                                                cellID = as.numeric(cellID)), 
                           panel = T, allow_unbalanced_panel = T)

ggdid(aggte(did.plainvanilla, type= "dynamic"))

ggdid(aggte(did.onlytreated.chunk1, type = "dynamic"))

did.plainvanilla.forforests <- att_gt(yname = "biomass", tname = "Year", idname = "cellID", 
                           gname = "treat.year", 
                           data = carb_data_cafr0030 %>% mutate(Year = as.numeric(Year),
                                                                cellID = as.numeric(cellID)) %>% 
                             filter(!cellID %in% non.forest.pixels), 
                           panel = T, allow_unbalanced_panel = T)

did.classic.wcovs <- att_gt(yname = "biomass", tname = "Year", idname = "cellID", 
                           gname = "treat.year", 
                           data = carb_data_cafr0030 %>% mutate(Year = as.numeric(Year),
                                                                cellID = as.numeric(cellID)) %>% 
                             select(biomass, Year, cellID, treat.year, cov.wo.clm.practicerun) %>% mutate(fownership = as.factor(fownership),
                                                                                                          nlcd = as.factor(nlcd)),
                           xformla = as.formula(paste0("~", paste(cov.wo.clm.practicerun, collapse = "+"))), 
                           est_method = "ipw",
                           panel = T, allow_unbalanced_panel = T)

ggdid(aggte(did.classic.wcovs, type= "dynamic"))
ggdid(aggte(did.plainvanilla.forforests, type= "dynamic"))

did.classic.wcovs.allchunk1 <- att_gt(yname = "biomass", tname = "Year", 
                                      idname = "cellID", 
                            gname = "treat.year", 
                            data = carb_data_chunk_1 %>% mutate(Year = as.numeric(Year),
                                                                 cellID = as.numeric(cellID)) %>% 
                              select(biomass, Year, cellID, treat.year, cov.wo.clm.practicerun) %>% mutate(fownership = as.factor(fownership),
                                                                                                           nlcd = as.factor(nlcd)),
                            xformla = as.formula(paste0("~", paste(cov.wo.clm.practicerun, collapse = "+"))), 
                            est_method = "ipw",
                            panel = T, allow_unbalanced_panel = T)

ggdid(aggte(did.classic.wcovs.allchunk1, type= "dynamic"))

# 
# 
# evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
# evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
# evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
# evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]
# 
# cov.wo.clm.forminus <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
#                          "nlcd", 'fownership', "distance.to.road",
#                          c(evi.cov.4, evi.cov.5),
#                          'biomass_tminus4', 'biomass_tminus5', "forest.group")
# 
# factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")


# Demo example from doubleML documentation =====

time.periods <- 4
sp <- reset.sim()
sp$te <- 0

set.seed(1814)

# generate dataset with 4 time periods
time.periods <- 4

# add dynamic effects
sp$te.e <- 1:time.periods

# generate data set with these parameters
# here, we dropped all units who are treated in time period 1 as they do not help us recover ATT(g,t)'s.
dta <- build_sim_dataset(sp)

# How many observations remained after dropping the ``always-treated'' units
nrow(dta)
#This is what the data looks like
head(dta)

set.seed(1234)
doubleml_did_linear <- function(y1, y0, D, covariates,
                                ml_g = lrn("regr.lm"),
                                ml_m = lrn("classif.log_reg"),
                                n_folds = 10, n_rep = 1, ...) {
  
  # warning if n_rep > 1 to handle mapping from psi to inf.func
  if (n_rep > 1) {
    warning("n_rep > 1 is not supported.")
  }
  # Compute difference in outcomes
  delta_y <- y1 - y0
  # Prepare data backend
  dml_data = DoubleML::double_ml_data_from_matrix(X = covariates, y = delta_y, d = D)
  # Compute the ATT
  dml_obj = DoubleML::DoubleMLIRM$new(dml_data, ml_g = ml_g, ml_m = ml_m, score = "ATTE", n_folds = n_folds)
  dml_obj$fit()
  att = dml_obj$coef[1]
  # Return results
  inf.func <- dml_obj$psi[, 1, 1]
  output <- list(ATT = att, att.inf.func = inf.func)
  return(output)
}

example_attgt_dml_linear <- att_gt(yname = "Y",
                                   tname = "period",
                                   idname = "id",
                                   gname = "G",
                                   xformla = ~X,
                                   data = dta,
                                   est_method = doubleml_did_linear)


summary(example_attgt_dml_linear)

