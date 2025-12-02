library(tidyverse);library(splines);library(tictoc);library(reshape2)
library(DoubleML);library(mlr3);library(mlr3learners);library(mlr3tuning)
library(sf);library(paletteer);library(lubridate)
library(ranger);library(scales);library(DRDID)

# Load data and define some variables =======

load("data/output/StepN2.9_CARBData_LHSwTreatYearsClustersandRHS_chunk1.RData")
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

# Data construction, in two functions =======

dml.datamaker.for.didclust.new.step1 <- function(proj.df.from.chunk, 
                                           #year.incr,
                                           # proj.dat.to.use, 
                                           cov.to.use, 
                                           # outcome.to.use,
                                           nonforest.exclude = T, 
                                           exclude.treat.pixels = T) {
  
  #the new data frame already has "treated year" and all the biomass data from prior & post years
    #so no need to use year.incr - which was used to construct data for each specific year
    #instead - step1: final data cleaning
  
  
  # year.of.treat <- year(proj.years[proj.years$project.name == proj.to.look.at, ]$DATE.first)
  year.of.treat <- unique(year(proj.df.from.chunk$DATE.first))
  #point of treatment
  
  dat.to.use <- proj.df.from.chunk[,which(colnames(proj.df.from.chunk) %in% c("cellID", "type", "Year", "treat",
                                                                              "lon", "lat",
                                                                              cov.to.use, "biomass", "biomass_0") | grepl("cluster|biomass_tminus", colnames(proj.df.from.chunk)))] 
  #biomass: biomass at that given year
    #biomass_0: biomass at the start of the year
  #add all cluster columns, as you do not know which one of these would be useful
  
  dat.to.use$forest.or.not <- !dat.to.use$cellID %in% non.forest.pixels 
  #TRUE: is a forest
  
  #dat.to.use[[cluster.cols]] <- proj.dat.to.use[proj.dat.to.use$proj==proj.to.look.at,][[cluster.cols]]
  
  if (nonforest.exclude) {
    dat.to.use <- subset(dat.to.use, forest.or.not)
  }
  
  cellids.included <- unique(dat.to.use$cellID)
  
  dat.to.use$biomass_tminus5.delta <- with(dat.to.use, biomass-biomass_tminus5)
  #difference of that year's biomass with tminus 5 biomass
  
  dat.to.use$biomass_tminus4.delta <- with(dat.to.use, biomass-biomass_tminus4)
  dat.to.use$biomass_tminus3.delta <- with(dat.to.use, biomass-biomass_tminus3)
  dat.to.use$biomass_tminus2.delta <- with(dat.to.use, biomass-biomass_tminus2)
  dat.to.use$biomass_tminus1.delta <- with(dat.to.use, biomass-biomass_tminus1)
  #repeated
  
  dat.to.use$biomass_0.delta <- with(dat.to.use, biomass-biomass_0)
  #difference of that year's biomass with the first year of issuance's biomass
  
  dat.to.use$biomass_tminus5.delta <- with(dat.to.use, biomass-biomass_tminus5)
  #difference of that year's biomass with tminus 5 biomass
  
  if (year.of.treat == 2004) {
    dat.to.use$biomass_tminus5.delta <- 0
    dat.to.use$biomass_tminus4.delta <- 0
  }
  
  if (year.of.treat == 2005) {
    dat.to.use$biomass_tminus5.delta <- 0
  }
  
  # if (outcome.to.use %in% c("conversion.cumul.ratio", "conversion.thisyear.ratio")) {
  #   
  #   if (year.incr > 0) {
  #     lcms.ratio.dat.for.this.cellids <- cells.in.refor.and.ac.projs.unique[cells.in.refor.and.ac.projs.unique$cellID %in% cellids.included,]
  #     
  #     lcms.ratio.dat.for.this.cellids <- lcms.ratio.dat.for.this.cellids[,grepl("cellID|lcms_",colnames(lcms.ratio.dat.for.this.cellids))]
  #     
  #     lcms.years <- as.integer(gsub("lcms_ratio_", "", colnames(lcms.ratio.dat.for.this.cellids)[-1]))
  #     
  #     cols.to.include.from1toyearincr <- lcms.years <= (year.of.treat + year.incr) & lcms.years >= year.of.treat+1
  #     #for example: if the year.of.treat was 2005, and year.incr = 1, then include all columns from 2005 to 2006
  #     
  #     cols.to.include.thisyear <- lcms.years == (year.of.treat + year.incr)
  #     #for example: if the year.of.treat was 2005, and year.incr = 1, then include column from 2006
  #     
  #     lcms.ratio.dat.for.this.cellids.from1toyearincr <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.from1toyearincr))]
  #     lcms.ratio.dat.for.thisyear <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.thisyear))]
  #     
  #     if (ncol(lcms.ratio.dat.for.this.cellids.from1toyearincr)>2) {
  #       conversion.ratio.cumul.fromlcms <- unname(rowSums(lcms.ratio.dat.for.this.cellids.from1toyearincr[, -1], na.rm = T))
  #     } else if (ncol(lcms.ratio.dat.for.this.cellids.from1toyearincr)==2) {
  #       conversion.ratio.cumul.fromlcms <- unname(unlist(lcms.ratio.dat.for.this.cellids.from1toyearincr[, -1]))
  #     }
  #     #whether the pixel has experienced conversion at least once ever since year.of.treat+1 to year.of.treat + year.incr
  #     
  #     conversion.ratio.thisyear <- lcms.ratio.dat.for.thisyear[[2]]
  #     #whether the pixel experienced conversion IN THIS YEAR
  #     
  #     lcms.df.to.fuse <- data.frame(cellID = lcms.ratio.dat.for.this.cellids.from1toyearincr$cellID,
  #                                   conversion.cumul.ratio = as.numeric(conversion.ratio.cumul.fromlcms),
  #                                   conversion.thisyear.ratio = as.numeric(conversion.ratio.thisyear))
  #     
  #     lcms.df.to.fuse$cellID <- as.character(lcms.df.to.fuse$cellID)
  #     
  #     dat.to.use <- left_join(dat.to.use, lcms.df.to.fuse, by = "cellID")
  #   }
  #   
  #   if (year.incr < 0) {
  #     #for year.incr < 0, do not calculate *cumulative* change, but changes in each year
  #     
  #     lcms.ratio.dat.for.this.cellids <- cells.in.refor.and.ac.projs.unique[cells.in.refor.and.ac.projs.unique$cellID %in% cellids.included,]
  #     
  #     lcms.ratio.dat.for.this.cellids <- lcms.ratio.dat.for.this.cellids[,grepl("cellID|lcms_",colnames(lcms.ratio.dat.for.this.cellids))]
  #     
  #     lcms.years <- as.integer(gsub("lcms_ratio_", "", colnames(lcms.ratio.dat.for.this.cellids)[-1]))
  #     
  #     cols.to.include.thisyear <- lcms.years == (year.of.treat + year.incr)
  #     #for example: if the year.of.treat was 2005, and year.incr = -3, then include column from 2002
  #     
  #     lcms.ratio.dat.for.thisyear <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.thisyear))]
  #     
  #     #whether the pixel has experienced conversion at least once ever since year.of.treat+1 to year.of.treat + year.incr
  #     
  #     conversion.ratio.thisyear <- lcms.ratio.dat.for.thisyear[[2]]
  #     #whether the pixel experienced conversion IN THIS YEAR
  #     
  #     lcms.df.to.fuse <- data.frame(cellID = lcms.ratio.dat.for.thisyear$cellID,
  #                                   conversion.cumul.ratio = as.numeric(conversion.ratio.thisyear),
  #                                   conversion.thisyear.ratio = as.numeric(conversion.ratio.thisyear))
  #     
  #     lcms.df.to.fuse$cellID <- as.character(lcms.df.to.fuse$cellID)
  #     
  #     dat.to.use <- left_join(dat.to.use, lcms.df.to.fuse, by = "cellID")
  #   }
  # }
  
  dat.to.use$treat <- as.integer(as.logical(dat.to.use$treat))
  
  if (exclude.treat.pixels) {
    control.but.treated <- which(grepl("control", dat.to.use$type) & dat.to.use$cellID %in% treated.pixels)
    
    if (length(control.but.treated)!=0) {
      dat.to.use <- dat.to.use[-control.but.treated,]
    }
    
  }
  
  if (unique(proj.df.from.chunk$projectID) =="CAFR0001" & any(grepl(paste(evi.cov.4, collapse = "|"), cov.to.use))) {
    dat.to.use[,which(colnames(dat.to.use) %in% evi.cov.5)] <- 0
    
    dat.to.use$biomass_tminus5 <- 0
    
    #CAFR0001 is the only project with NA ndvi.5 data (just make this data completely irrelevant)
  }
  
  dat.to.use <- dat.to.use[complete.cases(dat.to.use),]
  #drop all missing data
  
  for (fc in factor.covs) {
    if (fc %in% colnames(dat.to.use)) {
      dat.to.use[[fc]] <- as.factor(dat.to.use[[fc]])
    }
  }
  
  dat.to.use %>% 
    mutate(year.to.treat = as.numeric(Year) - year.of.treat) %>%
    mutate(projectID = unique(proj.df.from.chunk$projectID)) %>% 
    relocate(projectID, cellID, Year, year.to.treat) -> dat.to.use
  
  return(dat.to.use)
  
}

dml.datamaker.for.didclust.new.step2 <- function(dat.to.use, 
                                                 # dupl.drop = T, 
                                                 n.fold = 3, 
                                                 cluster.use = F, 
                                                 cluster.cols = NA) {
  
  #step 2 - split the data frame into *list*
    #Using cluster-based random split (if cluster.use = T)
      #if cluster.use = F, then full random split
    #list element 1 - "dfs_for_pscore"
    #list element 2 - "dfs_for_outcome"
    #Why two different elements? Because for pscore, only cross-sectional varition is meaningful
      #Either you got treated or not, and all the covariates are predetermined as of treatment
    #So for element 1 - just one entry per cellID
      #Element 2 - Multiple entries per cellID, one for each year
  
  dat.to.use.list <- list()

  dat.to.use %>% 
    distinct(cellID, .keep_all = T) -> dat.to.use.for.pscore
  #make this into cross-sectional
  #and then this will be divvied up using cluster column
  
  ### setting splits for each fold ======
  
  if (cluster.use) {
    #deciding how the clusters will be split
    cols.to.drop <- colnames(dat.to.use)[grepl('cluster', colnames(dat.to.use))]
    
    cols.to.drop <- setdiff(cols.to.drop, cluster.cols)
    
    dat.to.use <- select(dat.to.use, -cols.to.drop)
    #drop the other "cluster" columns to avoid confusion
    
    crosstable.cluster.v.treat <- table(dat.to.use.for.pscore[[cluster.cols]], dat.to.use.for.pscore[["treat"]])
    #use dat.to.use.for.pscore just to make the calculation faster
    
    clusters.w.treat <- rownames(crosstable.cluster.v.treat)[which(crosstable.cluster.v.treat[,2]!=0)]
    #clusters with at least one treated pixel
    clusters.control <- setdiff(rownames(crosstable.cluster.v.treat), clusters.w.treat)
    
    n.fold.to.be.used <- if (n.fold > length(clusters.w.treat)) {
      #if there are more folds asked than the number of clusters with treatment (e.g., only 2 clusters with treatment)
          #Happens for exceptionally small treatment areas (small projects)
        #Then just use the length of clusters with treatment
      length(clusters.w.treat)
    } else {
      n.fold
      #otherwise use the folds requested
    }
    
    sample.treat.clusters <- split(sample(clusters.w.treat, replace =F), 1:n.fold.to.be.used)
    sample.control.clusters <- split(sample(clusters.control, replace = F), 1:n.fold.to.be.used)
    
  } else {
    #if clusters are not being used
    
    dat.to.use <- select(dat.to.use, !contains("cluster"))
    #drop all "cluster" columns
    
    cellIDs.treat <- unique(subset(dat.to.use, treat ==1)$cellID)
    cellIDs.control <- unique(subset(dat.to.use, treat ==0)$cellID)
    
    # dat.to.use.treat.rows <- which(dat.to.use$treat ==1)
    # dat.to.use.control.rows <- which(dat.to.use$treat==0)
    
    sample.treat.cellIDs <- split(sample(cellIDs.treat, replace = F), 1:n.fold)
    sample.control.cellIDs <- split(sample(cellIDs.control, replace = F), 1:n.fold)
  }
  
  if (!cluster.use) {
    n.fold.to.be.used <- n.fold
  }
  
  ### setting the dfs_for_pscore ======
  
  dat.to.use.list[["dfs_for_pscore"]] <- list()
  
  if (cluster.use) {
    
    for (k in 1:n.fold.to.be.used) {
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]] <- list()
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]]$nuisance.fit.data <- rbind(dat.to.use.for.pscore[which(dat.to.use.for.pscore[[cluster.cols]] %in% setdiff(clusters.w.treat, sample.treat.clusters[[k]])),],
                                                                                           dat.to.use.for.pscore[which(dat.to.use.for.pscore[[cluster.cols]] %in% setdiff(clusters.control, sample.control.clusters[[k]])),])
      #data for fitting nuisance parameters
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]]$theta.fit.data <- rbind(dat.to.use.for.pscore[which(dat.to.use.for.pscore[[cluster.cols]] %in% sample.treat.clusters[[k]]),],
                                                                                                   dat.to.use.for.pscore[which(dat.to.use.for.pscore[[cluster.cols]] %in% sample.control.clusters[[k]]),])
      #data that will be used to calculate the theta (use the nuisance parameter fitted in the previous step)
    }
  } else {
    #not using clusters
    
    for (k in 1:n.fold.to.be.used) {
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]] <- list()
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]]$nuisance.fit.data <- rbind(subset(dat.to.use.for.pscore,
                                                                                                  cellID %in%  setdiff(cellIDs.treat, sample.treat.cellIDs[[k]])),
                                                                                           subset(dat.to.use.for.pscore,
                                                                                                  cellID %in%  setdiff(cellIDs.control, sample.control.cellIDs[[k]]))
                                                                                           )
      #data for fitting nuisance parameters
      
      dat.to.use.list[["dfs_for_pscore"]][[paste0("fold.", k)]]$theta.fit.data <- rbind(subset(dat.to.use.for.pscore, cellID %in% sample.treat.cellIDs[[k]]),
                                                                                        subset(dat.to.use.for.pscore, cellID %in% sample.control.cellIDs[[k]]))
      #data that will be used to calculate the theta (use the nuisance parameter fitted in the previous step)
    }
    
  }
  
  print(paste0("DONE MAKING PSCORE DATA!!"))
  
  ### setting the dfs_for_outcome =====
  
  list.of.years.to.treat <- unique(dat.to.use$year.to.treat)
  
  dat.to.use.list[["dfs_for_outcome"]] <- list()
  
  for (y in list.of.years.to.treat) {
    
    dat.to.use.for.this.outcome.y <- subset(dat.to.use, year.to.treat == y)
    
    dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]] <- list()
    
    
    if (cluster.use) {
      
      for (k in 1:n.fold.to.be.used) {
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]] <- list()
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]]$nuisance.fit.data <- rbind(dat.to.use.for.this.outcome.y[which(dat.to.use.for.this.outcome.y[[cluster.cols]] %in% setdiff(clusters.w.treat, sample.treat.clusters[[k]])),],
                                                                                             dat.to.use.for.this.outcome.y[which(dat.to.use.for.this.outcome.y[[cluster.cols]] %in% setdiff(clusters.control, sample.control.clusters[[k]])),])
        #data for fitting nuisance parameters
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]]$theta.fit.data <- rbind(dat.to.use.for.this.outcome.y[which(dat.to.use.for.this.outcome.y[[cluster.cols]] %in% sample.treat.clusters[[k]]),],
                                                                                                     dat.to.use.for.this.outcome.y[which(dat.to.use.for.this.outcome.y[[cluster.cols]] %in% sample.control.clusters[[k]]),])
        #data that will be used to calculate the theta (use the nuisance parameter fitted in the previous step)
      }
    } else {
      #not using clusters
      
      for (k in 1:n.fold.to.be.used) {
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]] <- list()
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]]$nuisance.fit.data <- rbind(subset(dat.to.use.for.this.outcome.y,
                                                                                                    cellID %in%  setdiff(cellIDs.treat, sample.treat.cellIDs[[k]])),
                                                                                             subset(dat.to.use.for.this.outcome.y,
                                                                                                    cellID %in%  setdiff(cellIDs.control, sample.control.cellIDs[[k]]))
        )
        #data for fitting nuisance parameters
        
        dat.to.use.list[["dfs_for_outcome"]][[paste0("year.to.treat_", y)]][[paste0("fold.", k)]]$theta.fit.data <- rbind(subset(dat.to.use.for.this.outcome.y, cellID %in% sample.treat.cellIDs[[k]]),
                                                                                          subset(dat.to.use.for.this.outcome.y, cellID %in% sample.control.cellIDs[[k]]))
        #data that will be used to calculate the theta (use the nuisance parameter fitted in the previous step)
      }
      
    }
  }
  
  return(dat.to.use.list)
}

## Showcasing data construction & successful stratified randomization =======

carb_data_cafr0030_cleaned <- dml.datamaker.for.didclust.new.step1(carb_data_cafr0030, 
                                                                   cov.to.use = cov.wo.clm.forminus,  
                                                                   nonforest.exclude = T, exclude.treat.pixels = T
                                                                   )

carb_data_cafr0030_cleaned_list <- dml.datamaker.for.didclust.new.step2(dat.to.use = carb_data_cafr0030_cleaned,
                                                                        n.fold = 3, 
                                                                        cluster.use = T,
                                                                        cluster.cols = "cluster.25km")

table(carb_data_cafr0030_cleaned_list$dfs_for_pscore$fold.1$nuisance.fit.data$cluster.25km)
table(carb_data_cafr0030_cleaned_list$dfs_for_pscore$fold.1$theta.fit.data$cluster.25km)
#see how these two have no clusters overlapping

carb_data_cafr0030_cleaned_list_woclust <- dml.datamaker.for.didclust.new.step2(dat.to.use = carb_data_cafr0030_cleaned,
                                                                        n.fold = 3, 
                                                                        cluster.use = F,
                                                                        cluster.cols = "cluster.25km")

table(carb_data_cafr0030_cleaned_list_woclust$dfs_for_pscore$fold.1$theta.fit.data$cluster.25km)
table(carb_data_cafr0030_cleaned_list_woclust$dfs_for_pscore$fold.1$nuisance.fit.data$cluster.25km)
#these two have overlapping clusters - because it was not a stratified randomization

# Fitting DMLDiD to get "ell_hat" and "pscore_hat" =========

dml.runner.for.didclust.new.baselinestr <- function(dat.list.to.use, cov.to.use, 
                                                    # outcome.to.use,
                                                    # baseline.string, 
                                                    # diff.with.baseline = T, 
                                                    #difference the tminus_5 with baseline results?
                                                    # diff.with.baseline.vars = c("biomass_tminus5"),
                                                    g.learner = c("classif.rf", "classif.logit", "regr.ols"), 
                                                    m.learner = c("regr.rf", "regr.ols"),
                                                    trim.use = F) {
  
  #Addition 251123:
    # Dropped "difference with baseline" and "difference with baseline" vars thing
    # Seemed very superfluous
  
  #Addition 240810:
  #This function is strictly for DID cluster method
  #Has the option to either use cluster or not, but does NOT allow the classic, non-DID score functions
  #For that, a separate function will be developed
  #This means that I no longer contain diff.use
  #cluster.use is also no longer there, because cluster-based cross fitting.sample splitting is already done beforehand
  
  #A note on outcomes:
  #I do NOT difference the outcomes (outcome in year t [initial year + year.incr] versus outcome in year 0 [e.g., "biomass_start"]) out BEFORE
  #This is because I want to make sure that the outcomes themselves are in Y_{it} terms, even though there are cases when the models are being fitted on Y_{it} - Y_{i0}, as in the exapmle of "ell" models
  
  #On the contrary:
  #If I am going to use "conversion.cumul.ratio" or "conversion.this.year" as the outcome, then I should not be using DMLDiD
  #I should instead be using repeated DMLs
  
  default.mtry.m <- ifelse(grepl("classif", m.learner), floor(sqrt(length(cov.to.use))),
                           floor((length(cov.to.use)+1)/3))
  default.mtry.g <- ifelse(grepl("classif", g.learner), floor(sqrt(length(cov.to.use))),
                           floor((length(cov.to.use)+1)/3))
  
  n.fold.used <- length(dat.list.to.use$dfs_for_pscore)
  
  #### Step 0. Adding differenced
  
  #a step to convert the biomass_0 and biomass_1 in *differences*
  #this step may be removed for non-DID methods
  
  # if (diff.with.baseline) {
  #   # for (k in 1:n.fold.used) {
  #   #   
  #   #   for (var in diff.with.baseline.vars) {
  #   #     dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[var]] <-(dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[var]]- dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[baseline.string]])
  #   #     
  #   #     dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[var]] <-(dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[var]] - dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[baseline.string]])  
  #   #   }
  #   #   
  #   # }
  # # }
  # 
  # vars.to.be.differenced <- c("biomass_tminus4", "biomass_tminus5")
  # 
  # if (any(vars.to.be.differenced %in% cov.to.use)) {
  #   for (k in 1:n.fold.used) {
  #     
  #     for (var in vars.to.be.differenced) {
  #       dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$nuisance.fit.data[[paste0(var, "_diff")]] <-(dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$nuisance.fit.data[[var]]- dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$nuisance.fit.data[[baseline.string]])
  #       
  #       dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$nuisance.fit.data[[paste0(var, "_diff")]] <-(dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$theta.fit.data[[var]] - dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$theta.fit.data[[baseline.string]])  
  #     }
  #     
  #   }
  # }
  
  ## Step 0. Making repositories for propensity scores and other fitted objects =====
  
  psi1.pre.repository <- list()
  # phat.repository <- list()
  D.repository <- list()
  #repositories of these two, for calculating asymptotic variance afterwards
  #we need to store them, because we need theta tilde to calculate the variance
  #and this is only revealed once we have calculated all the g, ell, and subsequently gotten the theta k's
  #note that psi1.pre has "pre" because we need to subtract "theta tilde" from the psi1.pre term to get the actual psi1
  
  thetak.repository <- c()
  
  g.and.ellhat.df.repository <- list()
  
  trimmed.N <- 0

  # for (k in 1:n.fold.used) {
  #   
  #   dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data <- dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[complete.cases(dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data),]
  #   dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data <- dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[complete.cases(dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data),]
  #   
  # }
  
  
  
  ## Step 1. Calculating propensity score (once) =============
  # 
  dat.to.use.for.pscore <- bind_rows(dat.list.to.use$dfs_for_pscore$fold.1$nuisance.fit.data,
                          dat.list.to.use$dfs_for_pscore$fold.1$theta.fit.data)
  # #full cross-sectional data - for reference
  
  pscore_repository <- list()
  
  for (k in 1:n.fold.used) {
    
    crosssectional.dat.this.fold <- dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]
    #
    
    crosssectional.dat.this.fold.nu <- crosssectional.dat.this.fold$nuisance.fit.data
    crosssectional.dat.this.fold.th <- crosssectional.dat.this.fold$theta.fit.data
    
    # dat.this.fold.nu.control <- subset(dat.this.fold.nu, treat == 0)
    #subset to control units (for ell hat)
    
    # dat.this.fold.nu.control$outcome.differenced <- dat.this.fold.nu.control[[outcome.to.use]] - dat.this.fold.nu.control[[baseline.string]]
    
    # dat.this.fold.th$outcome.differenced <- dat.this.fold.th[[outcome.to.use]] - dat.this.fold.th[[baseline.string]]
    
    cov.to.use.temp <- cov.to.use
    
    for (fc in factor.covs) {
      if (!all(c(unique(subset(crosssectional.dat.this.fold.nu, treat==0)[[fc]]) %in% unique(crosssectional.dat.this.fold.th[[fc]]),
                 unique(crosssectional.dat.this.fold.th[[fc]]) %in% unique(subset(crosssectional.dat.this.fold.nu, treat==0)[[fc]])))) {
        # print(paste0("REMOVE FC ", fc))

        cov.to.use.temp <- cov.to.use.temp[-which(cov.to.use.temp %in% fc)]

        # print(paste0("cov.to.use.temp:", paste(cov.to.use.temp, collapse = ", ")))
      }
      #if there is a factor variable that is existent in the theta data but not in the nuisance
      #just drop the factor variable from the covariate set
    }

    # for (fc in intersect(factor.covs, cov.to.use.temp)) {
    #   dat.this.fold.th[[fc]] <- as.factor(dat.this.fold.th[[fc]])
    #   dat.this.fold.nu[[fc]] <- as.factor(dat.this.fold.nu[[fc]])
    # }
    #
    
    phat.this.fold <- mean(crosssectional.dat.this.fold.nu$treat)
    
    if (g.learner == "classif.rf") {
      #RF classifier
      
      # ghat.this.fold.model <- forestry(x = data.frame(crosssectional.dat.this.fold.nu[, cov.to.use.temp]), 
      #                                  y = crosssectional.dat.this.fold.nu$treat,
      #                                  ntree = 500,
      #                                  mtry = default.mtry.g, 
      #                                  splitratio = 1)
      
      ghat.this.fold.model <- ranger(treat ~.,
                                     data = select(crosssectional.dat.this.fold.nu, c("treat",cov.to.use.temp)), 
                                     num.trees = 250,
                                     #reduce tree size to make it less computationally intensive
                                     min.node.size = 1,
                                     mtry = default.mtry.g, 
                                     respect.unordered.factors = 'order',
                                     classification = F)
      
      ghat.this.fold <- predict(ghat.this.fold.model,
                                select(crosssectional.dat.this.fold.th, cov.to.use.temp))$predictions
      
      
    } else if (g.learner == "classif.logit") {
      #Logit classifier
      
      ghat.this.fold.model <- glm(treat ~ .,
                                  data = select(crosssectional.dat.this.fold.nu, c(cov.to.use.temp, "treat")),
                                  family = "binomial")
      
      ghat.this.fold <- predict(ghat.this.fold.model,
                                select(crosssectional.dat.this.fold.th, cov.to.use.temp), type = "response")
      
      
    }
    
    crosssectional.dat.this.fold.th$ghat <- ghat.this.fold
    
    crosssectional.dat.this.fold.th$phat <- mean(crosssectional.dat.this.fold.nu$treat)
    
    pscore_repository[[paste0("fold.", k)]] <- crosssectional.dat.this.fold.th %>% 
      ungroup() %>% 
      select(cellID, ghat, phat, treat)
    
    #removed all parts that have to do with calculating theta's and what not
    
  }
  
  pscore_all <- bind_rows(pscore_repository, .id = "fold")
  
  ## Step 2. Calculating outcome model (for individual years) =================
  
  years <- names(dat.list.to.use$dfs_for_outcome)
  
  all_ell_results <- list()
  
  for (ylabel in years) {
    curr_year_val <- as.numeric(gsub("year.to.treat_", "", ylabel))
    
    print(paste0("RUNNING ELLHAT MODEL FOR YEAR ", curr_year_val))
    
    this_year_ell <- list()
    
    for (k in 1:n.fold.used) {
      
      outcome.thisy.dat.this.fold <- dat.list.to.use$dfs_for_outcome[[ylabel]][[paste0("fold.", k)]]
      #
      
      outcome.thisy.dat.this.fold.nu <- outcome.thisy.dat.this.fold$nuisance.fit.data
      outcome.thisy.dat.this.fold.th <- outcome.thisy.dat.this.fold$theta.fit.data
      
      # dat.this.fold.nu.control <- subset(dat.this.fold.nu, treat == 0)
      #subset to control units (for ell hat)
      
      # dat.this.fold.nu.control$outcome.differenced <- dat.this.fold.nu.control[[outcome.to.use]] - dat.this.fold.nu.control[[baseline.string]]
      
      # dat.this.fold.th$outcome.differenced <- dat.this.fold.th[[outcome.to.use]] - dat.this.fold.th[[baseline.string]]
      
      cov.to.use.temp <- cov.to.use
      
      # for (fc in factor.covs) {
      #   if (!all(c(unique(dat.this.fold.nu.control[[fc]]) %in% unique(dat.this.fold.th[[fc]]),
      #              unique(dat.this.fold.th[[fc]]) %in% unique(dat.this.fold.nu.control[[fc]])))) {
      #     # print(paste0("REMOVE FC ", fc))
      #     
      #     cov.to.use.temp <- cov.to.use.temp[-which(cov.to.use.temp %in% fc)]
      #     
      #     # print(paste0("cov.to.use.temp:", paste(cov.to.use.temp, collapse = ", ")))
      #   }
      #   #if there is a factor variable that is existent in the theta data but not in the nuisance
      #   #just drop the factor variable from the covariate set
      # }
      
      # for (fc in intersect(factor.covs, cov.to.use.temp)) {
      #   dat.this.fold.th[[fc]] <- as.factor(dat.this.fold.th[[fc]])
      #   dat.this.fold.nu[[fc]] <- as.factor(dat.this.fold.nu[[fc]])
      # }
      #
      
      outcome.thisy.dat.this.fold.nu$outcome.differenced <- outcome.thisy.dat.this.fold.nu$biomass_tminus2.delta
      outcome.thisy.dat.this.fold.th$outcome.differenced <- outcome.thisy.dat.this.fold.th$biomass_tminus2.delta
      #biomass_tminus2.delta = for each cellID-year's biomass, biomass_{cellyear} - biomass_{cell, year = t-2}
      
      outcome.thisy.dat.this.fold.nu.control <- subset(outcome.thisy.dat.this.fold.nu, treat ==0)
      
      
      if (m.learner == "regr.rf") {
        
        ellhat.this.fold.model <- ranger(outcome.differenced ~.,
                                         data = select(outcome.thisy.dat.this.fold.nu.control, c("outcome.differenced",cov.to.use.temp)),
                                         num.trees = 250,
                                         respect.unordered.factors = 'order',
                                         mtry = default.mtry.m)
        
        ellhat.this.fold <- predict(ellhat.this.fold.model,
                                    select(outcome.thisy.dat.this.fold.th, cov.to.use.temp))$predictions
        
      } else if (m.learner == "regr.ols") {
        ellhat.this.fold.model <- lm(y ~.,
                                     data = select(outcome.thisy.dat.this.fold.nu.control, c(cov.to.use.temp, y = outcome.differenced)))
        
        ellhat.this.fold <- predict(ellhat.this.fold.model,
                                    select(outcome.thisy.dat.this.fold.th, cov.to.use.temp))
      }
      
      outcome.thisy.dat.this.fold.th$ellhat <- ellhat.this.fold
      
      this_year_ell[[paste0("fold.", k)]] <- outcome.thisy.dat.this.fold.th %>% 
        mutate(year_val = curr_year_val) %>% 
        ungroup() %>% 
        select(cellID, year_val, outcome.differenced, ellhat)
      
      #removed all parts that have to do with calculating theta's and what not
      
    }
    
    this_year_ell_bound <- bind_rows(this_year_ell)
    
    all_ell_results[[ylabel]] <- this_year_ell_bound
    
    
  }
  
  all_ell_results_bounded <- bind_rows(all_ell_results)
  
  pscore_and_ell <- full_join(pscore_all,
                              all_ell_results_bounded,
                              by = "cellID")
  
  
  ## Step 3. Trimming (if requested) & calculating ATT point estimates (for each fold-year) =========
  
  if (trim.use) {
    pscore_and_ell %>% 
      filter(ghat < 0.99 & ghat > 0.01) -> pscore_and_ell_trimmed
    # ellhat.this.fold.trimmed <- ellhat.this.fold[ghat.this.fold < 0.99 & ghat.this.fold > 0.01]
    
    # dat.this.fold.th.trimmed <- dat.this.fold.th[ghat.this.fold < 0.99 & ghat.this.fold > 0.01, ]
  } else {
    pscore_and_ell_trimmed <- pscore_and_ell
  }
  
  pscore_and_ell_trimmed_split <- split(pscore_and_ell_trimmed,
                                        pscore_and_ell_trimmed$year_val)
  
  
  for (ylabel in years) {
    
    curr_year_val <- as.numeric(gsub("year.to.treat_", "", ylabel))
    
    pscore_and_ell_trimmed_df_for_this_y <- subset(pscore_and_ell_trimmed, year_val == curr_year_val)
    
    pscore_and_ell_trimmed_df_for_this_y_list <- list()
    
    for (k in 1:n.fold.used) {
      
      pscore_and_ell_trimmed_df_for_this_y.this.fold <- subset(pscore_and_ell_trimmed_df_for_this_y,
                                                               fold == paste0("fold.", k))
      
      pscore_and_ell_trimmed_df_for_this_y.this.fold$psi1.pre.this.fold <- with(pscore_and_ell_trimmed_df_for_this_y.this.fold,
                                                                                (treat - ghat)/phat/(1-ghat)*(outcome.differenced - ellhat)) 
      
      pscore_and_ell_trimmed_df_for_this_y.this.fold$thetak.this.fold <- mean(pscore_and_ell_trimmed_df_for_this_y.this.fold$psi1.pre.this.fold)
      
      # thetak.repository <- c(thetak.repository, thetak.this.fold)
      
      pscore_and_ell_trimmed_df_for_this_y_list[[paste0("fold.",k)]] <-  pscore_and_ell_trimmed_df_for_this_y.this.fold
      
      print(paste0("DONE WITH FOLD NO. ", k, "/", n.fold.used))
      
    }
    
    pscore_and_ell_trimmed_df_for_this_y <- bind_rows(pscore_and_ell_trimmed_df_for_this_y_list)
    
    pscore_and_ell_trimmed_split[[as.character(curr_year_val)]] <- pscore_and_ell_trimmed_df_for_this_y
    
  }
  
  pscore_and_ell_trimmed_withpsiandthetak <- bind_rows(pscore_and_ell_trimmed_split)
  
  pscore_and_ell_trimmed_withpsiandthetak %>% 
    distinct(fold, year_val, .keep_all = T) %>% 
    select(fold, year_val, thetak.this.fold) -> thetak_repository
  
  pscore_and_ell_trimmed_withpsiandthetak %>% 
    distinct(fold, .keep_all = T) %>% 
    select(fold, phat) -> phat_repository
  
  thetak_repository %>% 
    group_by(year_val) %>% 
    summarise(thetatilde = mean(thetak.this.fold)) -> thetatilde_repository
  
  trimmed.N <- length(unique(pscore_and_ell_trimmed_withpsiandthetak$cellID))
  
  
  ## Step 4. Calculating ATT variance =======
  
  pscore_and_ell_trimmed_withpsiandthetak_split <- split(pscore_and_ell_trimmed_withpsiandthetak,
                                                         pscore_and_ell_trimmed_withpsiandthetak$year_val)
  
  for (ylabel in years) {
    
    curr_year_val <- as.numeric(gsub("year.to.treat_", "", ylabel))
    
    pscore_and_ell_trimmed_withpsiandthetak_for_this_y <- subset(pscore_and_ell_trimmed_withpsiandthetak, 
                                                                 year_val == curr_year_val)
    
    pscore_and_ell_trimmed_withpsiandthetak_for_this_y_list <- list()
    
    for (k in 1:n.fold.used) {
      
      pscore_and_ell_trimmed_withpsiandthetak_for_this_y.this.fold <- subset(pscore_and_ell_trimmed_withpsiandthetak_for_this_y,
                                                               fold == paste0("fold.", k))
      
      psi1.pre.this.fold <-  pscore_and_ell_trimmed_withpsiandthetak_for_this_y.this.fold$psi1.pre.this.fold
      
      psi1.this.fold <- psi1.pre.this.fold - subset(thetatilde_repository, year_val==curr_year_val)$thetatilde
      
      phat.this.fold <- subset(phat_repository, fold == paste0("fold.", k))$phat
      
      Ghat.this.fold <- 0-subset(thetatilde_repository, year_val==curr_year_val)$thetatilde/phat.this.fold
      
      iff.func.this.fold <- psi1.this.fold + Ghat.this.fold*(pscore_and_ell_trimmed_withpsiandthetak_for_this_y.this.fold$treat - phat.this.fold)
      
      # to.be.meaned <- iff.func.this.fold^2
      #the sample mean of this guy is the asymptotic variance of $\sqrt{N}$
      
      pscore_and_ell_trimmed_withpsiandthetak_for_this_y.this.fold$iff.this.fold <- iff.func.this.fold
      
      pscore_and_ell_trimmed_withpsiandthetak_for_this_y_list[[paste0("fold.", k)]] <- pscore_and_ell_trimmed_withpsiandthetak_for_this_y.this.fold
    }
    
    pscore_and_ell_trimmed_withpsiandthetak_for_this_y <- bind_rows(pscore_and_ell_trimmed_withpsiandthetak_for_this_y_list)
    
    
    pscore_and_ell_trimmed_withpsiandthetak_split[[as.character(curr_year_val)]] <- pscore_and_ell_trimmed_withpsiandthetak_for_this_y
    
  }
  
  pscore_and_ell_trimmed_withpsiandthetakandiff <- bind_rows(pscore_and_ell_trimmed_withpsiandthetak_split)
  
  pscore_and_ell_trimmed_withpsiandthetakandiff %>% 
    left_join(thetatilde_repository, by = "year_val") %>%
    rename(att = thetatilde) -> pscore_and_ell_trimmed_withpsiandthetakandiff
  
  pscore_and_ell_trimmed_withpsiandthetakandiff %>% 
    left_join(dat.to.use.for.pscore %>% select(cellID, contains("cluster.")),
              by = "cellID") %>% 
    mutate(projectID = unique(dat.to.use.for.pscore$projectID)) %>% 
    relocate(projectID) -> pscore_and_ell_trimmed_withpsiandthetakandiff
  
  pscore_and_ell_trimmed_withpsiandthetakandiff %>% 
    group_by(projectID, year_val, fold) %>% 
    summarise(mean_iff2 = mean(tcrossprod(iff.this.fold), na.rm = T) 
              #mean of the (iff)(iff') at this fold = $\hat{\Sigma}_{1k} = $
           ) %>%
    group_by(projectID, year_val) %>% 
    summarise(Sigma = mean(mean_iff2)) %>% 
    left_join(pscore_and_ell_trimmed_withpsiandthetakandiff %>% 
                group_by(projectID, year_val) %>% 
                summarise(att = unique(att),
                          cluster_count = n_distinct(cluster.25km)),
              by = c("projectID", "year_val")) %>% 
    mutate(att.se = sqrt(Sigma/cluster_count))-> project_summary
 
  
  final.list <- list(pscore_and_ell_withIFF = pscore_and_ell_trimmed_withpsiandthetakandiff,
                     project_summary = project_summary)
  
  return(final.list)
  
}

# Showcasing usage for CAFR0030 =====

library(did)

carb_data_cafr0030_dmlresult <- dml.runner.for.didclust.new.baselinestr(carb_data_cafr0030_cleaned_list, cov.to.use = cov.wo.clm.forminus,  
                                        g.learner = "classif.rf", m.learner = "regr.rf", trim.use = T)

save(list = "carb_data_cafr0030_dmlresult",
     file = "data/output/StepN3_Showcase_CARB0030_Outcome.RData")

load("data/output/StepN3_Showcase_CARB0030_Outcome.RData")

lookup_dml_estimator <- function(covariates, ...) {
  # 1. Unpack values by COLUMN INDEX (Safest)
  att_vec <- covariates[, 'att']
  inf_vec <- covariates[, 'iff.this.fold']
  
  att <- mean(att_vec[att_vec!=0], na.rm = TRUE)
  #why drop 0? the compute.att_gt function (which is in the did::att_gt)
    #has a wonky feature of adding a bunch of zero-ATT rows
  
  
  # 3. Return BOTH names for compatibility
  # Note: The influence function vector MUST be the same length as the input data
  return(list(
    ATT = att,
    inf.func = inf_vec,
    att.inf.func = inf_vec
  ))
}

carb_data_cafr0030_dml_csapplied <- att_gt(yname = "biomass", tname = "Year", idname = "cellID", 
                                           gname = "treat.year", 
                                           data = carb_data_cafr0030_cleaned %>% select(-treat) %>%  left_join(carb_data_cafr0030_dmlresult$pscore_and_ell_withIFF %>% 
                                                                                                                 select(projectID, cellID, year_val, fold,ghat, phat, treat, outcome.differenced, ellhat, psi1.pre.this.fold, iff.this.fold, att),
                                                                                                               by = c("projectID", "cellID", "year.to.treat" = "year_val")) %>% 
                                             left_join(carb_data_cafr0030 %>% distinct(cellID, .keep_all =T) %>% select(cellID, treat.year),
                                                       by = "cellID") %>% 
                                             mutate(Year = as.numeric(Year),
                                                    cellID = as.numeric(cellID)) %>% 
                                             filter(!is.na(psi1.pre.this.fold) & !is.na(att)),
                                           xformla = ~ att + iff.this.fold + treat,
                                           est_method = lookup_dml_estimator,
                                           anticipation = 1,
                                           base_period = "universal",
                                           control_group = "nevertreated",
                                           #assumed 1-period anticipation
                                           clustervars = "cluster.25km",
                                           bstrap = T,
                                           cband = F,
                                           panel = T, 
                                           allow_unbalanced_panel = T)

ggdid(aggte(carb_data_cafr0030_dml_csapplied, type = "dynamic", clustervars = "cluster.25km"))

# Showcasing usage for all projects in chunk 1 =====

projects_in_chunk1 <- unique(carb_data_chunk_1$projectID)

chunk1_results_list <- list()

chunk1_results_list[["MLBased"]] <- list()

for (p in projects_in_chunk1) {
  
  this_p_data <- subset(carb_data_chunk_1, projectID == p)
  
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
  
  chunk1_results_list[["MLBased"]][[p]] <- this_p_data_dmlrun_result
  
  print(paste0("DONE WITH PROJECT ", p))
}

save(list = "chunk1_results_list",
     file = "data/output/StepN3_Showcase_CARB0030andOtherChunk1_Outcome.RData")

chunk1_results_combined <- bind_rows(chunk1_results_list$MLBased %>% 
                                       map("pscore_and_ell_withIFF") %>%
                                       bind_rows())

chunk1_results_combined2 <- bind_rows(chunk1_results_list$MLBased[1:2] %>% 
                                       map("pscore_and_ell_withIFF") %>%
                                       bind_rows())

carb_data_chunk1_dml_csapplied <- att_gt(yname = "biomass", tname = "Year", 
                                           idname = "projectID_cellID", 
                                           gname = "treat.year", 
                                           data = carb_data_chunk_1 %>% 
                                           filter(projectID %in% names(chunk1_results_list$MLBased[1:2])) %>% 
                                             mutate(year.to.treat = as.numeric(Year) - year(DATE.first)) %>% 
                                             select(-treat) %>%  
                                             left_join(chunk1_results_combined %>% 
                                                         filter(projectID %in% names(chunk1_results_list$MLBased[1:2])) %>% 
                                                         select(projectID, cellID, year_val, fold,ghat, phat, treat, outcome.differenced, ellhat, psi1.pre.this.fold, iff.this.fold, att),
                                                       by = c("projectID", "cellID", "year.to.treat" = "year_val")) %>%
                                             mutate(Year = as.numeric(Year),
                                                    projectID_cellID = dense_rank(paste0(projectID, "-", cellID)),
                                                    projectID_cluster = dense_rank(paste0(projectID, "-", cluster.25km)),
                                                    cellID = as.numeric(cellID)) %>% 
                                             filter(!is.na(psi1.pre.this.fold) & !is.na(att)),
                                           xformla = ~ att + iff.this.fold,
                                           est_method = lookup_dml_estimator,
                                           anticipation = 1,
                                           base_period = "universal",
                                           control_group = "nevertreated",
                                           #assumed 1-period anticipation
                                           clustervars = "projectID_cluster",
                                           bstrap = T,
                                           cband = F,
                                           panel = T, 
                                           allow_unbalanced_panel = T)

carb_data_chunk1_proj1_dml_csapplied <- att_gt(yname = "biomass", tname = "Year", 
                                         idname = "projectID_cellID", 
                                         gname = "treat.year", 
                                         data = carb_data_chunk_1 %>% 
                                           filter(projectID %in% names(chunk1_results_list$MLBased[1])) %>% 
                                           mutate(year.to.treat = as.numeric(Year) - year(DATE.first)) %>% 
                                           select(-treat) %>%  
                                           left_join(chunk1_results_combined %>% 
                                                       filter(projectID %in% names(chunk1_results_list$MLBased[1])) %>% 
                                                       select(projectID, cellID, year_val, fold,ghat, phat, treat, outcome.differenced, ellhat, psi1.pre.this.fold, iff.this.fold, att),
                                                     by = c("projectID", "cellID", "year.to.treat" = "year_val")) %>%
                                           mutate(Year = as.numeric(Year),
                                                  projectID_cellID = dense_rank(paste0(projectID, "-", cellID)),
                                                  projectID_cluster = dense_rank(paste0(projectID, "-", cluster.25km)),
                                                  cellID = as.numeric(cellID)) %>% 
                                           filter(!is.na(psi1.pre.this.fold) & !is.na(att)),
                                         xformla = ~ att + iff.this.fold,
                                         est_method = lookup_dml_estimator,
                                         anticipation = 1,
                                         base_period = "universal",
                                         control_group = "nevertreated",
                                         #assumed 1-period anticipation
                                         clustervars = "projectID_cluster",
                                         bstrap = T,
                                         cband = F,
                                         panel = T, 
                                         allow_unbalanced_panel = T)

carb_data_chunk1_proj2_dml_csapplied <- att_gt(yname = "biomass", tname = "Year", 
                                               idname = "projectID_cellID", 
                                               gname = "treat.year", 
                                               data = carb_data_chunk_1 %>% 
                                                 filter(projectID %in% names(chunk1_results_list$MLBased[2])) %>% 
                                                 mutate(year.to.treat = as.numeric(Year) - year(DATE.first)) %>% 
                                                 select(-treat) %>%  
                                                 left_join(chunk1_results_combined %>% 
                                                             filter(projectID %in% names(chunk1_results_list$MLBased[2])) %>% 
                                                             select(projectID, cellID, year_val, fold,ghat, phat, treat, outcome.differenced, ellhat, psi1.pre.this.fold, iff.this.fold, att),
                                                           by = c("projectID", "cellID", "year.to.treat" = "year_val")) %>%
                                                 mutate(Year = as.numeric(Year),
                                                        projectID_cellID = dense_rank(paste0(projectID, "-", cellID)),
                                                        projectID_cluster = dense_rank(paste0(projectID, "-", cluster.25km)),
                                                        cellID = as.numeric(cellID)) %>% 
                                                 filter(!is.na(psi1.pre.this.fold) & !is.na(att)),
                                               xformla = ~ att + iff.this.fold,
                                               est_method = lookup_dml_estimator,
                                               anticipation = 1,
                                               base_period = "universal",
                                               control_group = "nevertreated",
                                               #assumed 1-period anticipation
                                               clustervars = "projectID_cluster",
                                               bstrap = T,
                                               cband = F,
                                               panel = T, 
                                               allow_unbalanced_panel = T)

first_two_projs <- aggte(carb_data_chunk1_dml_csapplied, type = 'dynamic', na.rm = T)
first_first_proj <- aggte(carb_data_chunk1_proj1_dml_csapplied, type = 'dynamic', na.rm = T)
first_second_proj <- aggte(carb_data_chunk1_proj2_dml_csapplied, type = 'dynamic', na.rm = T)


ggdid(first_two_projs) + labs(title = "CAFR0001 and CAFR0002, combined via Callaway-Sant'Anna") + ylim(c(-7, 25))
ggdid(first_first_proj) + labs(title = "CAFR0001")  + ylim(c(-7, 25))
ggdid(first_second_proj) + labs(title = "CAFR0002") + ylim(c(-7, 25))



