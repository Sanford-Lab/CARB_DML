library(tidyverse);library(splines);library(tictoc);library(reshape2)
library(DoubleML);library(mlr3);library(mlr3learners);library(mlr3tuning)
library(sf);library(paletteer);library(lubridate)
library(ranger);library(scales);library(DRDID);library(Matrix)


# Calculation for chunk1: already done in the N3_Showcase code
  #Now we'll do chunks2, 3, ...

source("code/R_StepN3_DMLDiDwithClusters_ShowcasewithCARB0030_customfunctions.R")

load("data/output/StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk1.RData")
#carb_data_chunk_1
load("data/output/StepN2.5_SpatialClusters.RData")
load("data/output/StepN2.5_NonForestAndTreatedPixels.RData")

carb_data_cafr0030 <- subset(carb_data_chunk_1, projectID == "CAFR0030")

evi.cov.2 <- colnames(carb_data_chunk_1)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(carb_data_chunk_1)) & grepl("\\_2", colnames(carb_data_chunk_1)))]
evi.cov.3 <- colnames(carb_data_chunk_1)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(carb_data_chunk_1)) & grepl("\\_3", colnames(carb_data_chunk_1)))]
evi.cov.4 <- colnames(carb_data_chunk_1)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(carb_data_chunk_1)) & grepl("\\_4", colnames(carb_data_chunk_1)))]
evi.cov.5 <- colnames(carb_data_chunk_1)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(carb_data_chunk_1)) & grepl("\\_5", colnames(carb_data_chunk_1)))]

cov.wo.clm.forminus <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         "nlcd", 'fownership', "distance.to.road",
                         c(evi.cov.4, evi.cov.5),
                         'biomass_tminus4', 'biomass_tminus5',
                         "forest.group"
)

factor.covs <- c("fownership", "nlcd", "clm_ned_lf",
                 "forest.group"
)

tic()
toc()

for (chunkno in 5:7) {
  
  load(paste0("data/output/StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk", chunkno, ".RData"))
  #the name of this variable will be paste0("carb_data_chunk_", chunkno)
  
  this_chunk_carb_data <- get(paste0("carb_data_chunk_", chunkno))
  
  rm(list = paste0("carb_data_chunk_", chunkno))
  
  projects_in_this_chunk <- unique(this_chunk_carb_data$projectID)
  
  this_chunk_results_list <- list()
  
  this_chunk_results_list[["MLBased"]] <- list()
  
  for (p in projects_in_this_chunk) {
    
    this_p_data <- subset(this_chunk_carb_data, projectID == p)
    
    this_p_data_cleaned <- dml.datamaker.for.didclust.new.step1(this_p_data, 
                                                                cov.to.use = cov.wo.clm.forminus,  
                                                                nonforest.exclude = T, exclude.treat.pixels = T
    )
    
    this_p_data_cleaned_list <- dml.datamaker.for.didclust.new.step2(dat.to.use = this_p_data_cleaned,
                                                                     n.fold = 3, 
                                                                     cluster.use = T,
                                                                     cluster.cols = "cluster.25km")
    
    if (length(this_p_data_cleaned_list$dfs_for_pscore) <3) {
      print("SKIPPING THIS PROJECT - TREATMENT AREA IS TOO SMALL")
      
      next
    }
    
    this_p_data_dmlrun_result <- dml.runner.for.didclust.new.baselinestr(dat.list.to.use = this_p_data_cleaned_list, cov.to.use = cov.wo.clm.forminus,  
                                                                         g.learner = "classif.rf", m.learner = "regr.rf", trim.use = T)
    
    # this_p_data_dmlrunwithOLSLogit_result <- dml.runner.for.didclust.new.baselinestr(dat.list.to.use = this_p_data_cleaned_list, cov.to.use = cov.wo.clm.forminus,  
    #                                                                      g.learner = "classif.logit", m.learner = "regr.ols", trim.use = T)
    
    this_chunk_results_list[["MLBased"]][[p]] <- this_p_data_dmlrun_result
    
    print(paste0("DONE WITH PROJECT ", p))
  }
  
  assign(paste0("chunk", chunkno, "_results_list"), 
         this_chunk_results_list)
  
  rm(list = "this_chunk_results_list")
  
  tic()
  toc()
  
  save(list = paste0("chunk", chunkno, "_results_list"),
       file = paste0("data/output/StepN3_Showcase_CARB0030andOtherChunk", chunkno, "_Outcome.RData"))
}
