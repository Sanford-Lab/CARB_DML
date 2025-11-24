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
                         "forest.group")

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

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
    relocate(cellID, Year, year.to.treat) -> dat.to.use
  
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
                                                    baseline.string, 
                                                    # diff.with.baseline = T, 
                                                    #difference the tminus_5 with baseline results?
                                                    # diff.with.baseline.vars = c("biomass_tminus5"),
                                                    g.learner = c("classif.rf", "classif.logit", "regr.ols"), 
                                                    m.learner = c("regr.rf", "regr.ols"),
                                                    trim.use = F) {
  
  #Addition 251123:
    # Dropped "difference with baseline" and "difference with baseline" vars thing
    # Will just allow the difference with baseline and the original value itself
  
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
  
  #### Step 0. biomass_0 and biomass_1 differenced with biomass_start
  
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
  # }
  
  vars.to.be.differenced <- c("biomass_tminus4", "biomass_tminus5")
  
  if (any(vars.to.be.differenced %in% cov.to.use)) {
    for (k in 1:n.fold.used) {
      
      for (var in vars.to.be.differenced) {
        dat.list.to.use$dfs_for_pscore[[paste0("fold.", k)]]$nuisance.fit.data[[paste0(var, "_diff")]] <-(dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[var]]- dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[baseline.string]])
        
        dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[var]] <-(dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[var]] - dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[[baseline.string]])  
      }
      
    }
  }
  

  
  for (k in 1:n.fold.used) {
    
    dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data <- dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[complete.cases(dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data),]
    dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data <- dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data[complete.cases(dat.list.to.use[[paste0("fold.", k)]]$theta.fit.data),]
    
  }
  
  dat.to.use <- bind_rows(dat.list.to.use$fold.1$nuisance.fit.data,
                          dat.list.to.use$fold.1$theta.fit.data)
  #full data - this is where we will be adding the fitted g
  
  psi1.pre.repository <- list()
  phat.repository <- list()
  D.repository <- list()
  #repositories of these two, for calculating asymptotic variance afterwards
  #we need to store them, because we need theta tilde to calculate the variance
  #and this is only revealed once we have calculated all the g, ell, and subsequently gotten the theta k's
  #note that psi1.pre has "pre" because we need to subtract "theta tilde" from the psi1.pre term to get the actual psi1
  
  thetak.repository <- c()
  
  g.and.ellhat.df.repository <- list()
  
  trimmed.N <- 0
  
  for (k in 1:n.fold.used) {
    
    dat.this.fold <- dat.list.to.use[[paste0("fold.", k)]]
    
    dat.this.fold.nu <- dat.this.fold$nuisance.fit.data
    dat.this.fold.th <- dat.this.fold$theta.fit.data
    
    dat.this.fold.nu.control <- subset(dat.this.fold.nu, treat == 0)
    #subset to control units (for ell hat)
    
    dat.this.fold.nu.control$outcome.differenced <- dat.this.fold.nu.control[[outcome.to.use]] - dat.this.fold.nu.control[[baseline.string]]
    
    dat.this.fold.th$outcome.differenced <- dat.this.fold.th[[outcome.to.use]] - dat.this.fold.th[[baseline.string]]
    
    cov.to.use.temp <- cov.to.use
    
    for (fc in factor.covs) {
      if (!all(c(unique(dat.this.fold.nu.control[[fc]]) %in% unique(dat.this.fold.th[[fc]]),
                 unique(dat.this.fold.th[[fc]]) %in% unique(dat.this.fold.nu.control[[fc]])))) {
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
    #### Step 1. Calculating nuisance parameters
    
    phat.this.fold <- mean(dat.this.fold.nu$treat)
    
    if (g.learner == "classif.rf") {
      #RF classifier
      
      # ghat.this.fold.model <- forestry(x = data.frame(dat.this.fold.nu[, cov.to.use.temp]), 
      #                                  y = dat.this.fold.nu$treat,
      #                                  ntree = 500,
      #                                  mtry = default.mtry.g, 
      #                                  splitratio = 1)
      
      ghat.this.fold.model <- ranger(treat ~.,
                                     data = select(dat.this.fold.nu, c("treat",cov.to.use.temp)), 
                                     num.trees = 500,
                                     min.node.size = 1,
                                     mtry = default.mtry.g, 
                                     classification = F)
      
      ghat.this.fold <- predict(ghat.this.fold.model,
                                select(dat.this.fold.th, cov.to.use.temp))$predictions
      
      
    } else if (g.learner == "classif.logit") {
      #Logit classifier
      
      ghat.this.fold.model <- glm(treat ~ .,
                                  data = select(dat.this.fold.nu, c(cov.to.use.temp, "treat")),
                                  family = "binomial")
      
      ghat.this.fold <- predict(ghat.this.fold.model,
                                select(dat.this.fold.th, cov.to.use.temp), type = "response")
      
      
    } else if (g.learner == "regr.ols") {
      #LPM
      
      ghat.this.fold.model <- lm(treat ~ .,
                                 data = select(dat.this.fold.nu, c(cov.to.use.temp, "treat")))
      
      ghat.this.fold <- predict(ghat.this.fold.model,
                                select(dat.this.fold.th, cov.to.use.temp))
      
      
    }
    
    
    if (m.learner == "regr.rf") {
      # ellhat.this.fold.model <- forestry(x = data.frame(dat.this.fold.nu.control[, cov.to.use.temp]),
      #                                    y = dat.this.fold.nu.control[["outcome.differenced"]],
      #                                    ntree = 500,
      #                                    mtry = default.mtry.m,
      #                                    splitratio = 1)
      
      ellhat.this.fold.model <- ranger(outcome.differenced ~.,
                                       data = select(dat.this.fold.nu.control, c("outcome.differenced",cov.to.use.temp)),
                                       num.trees = 500,
                                       mtry = default.mtry.m)
      
      ellhat.this.fold <- predict(ellhat.this.fold.model,
                                  select(dat.this.fold.th, cov.to.use.temp))$predictions
      
    } else if (m.learner == "regr.ols") {
      ellhat.this.fold.model <- lm(y ~.,
                                   data = select(dat.this.fold.nu.control, c(cov.to.use.temp, y = outcome.differenced)))
      
      ellhat.this.fold <- predict(ellhat.this.fold.model,
                                  select(dat.this.fold.th, cov.to.use.temp))
    }
    
    
    if ("Year" %in% colnames(dat.to.use)) { # if this is a pooled estimate
      g.and.ellhat.df.repository[[paste0("fold.", k)]] <- data.frame(cellID = dat.this.fold.th$cellID,
                                                                     Year = dat.this.fold.th$Year,
                                                                     ghat = ghat.this.fold,
                                                                     ellhat = ellhat.this.fold)
    } else {
      g.and.ellhat.df.repository[[paste0("fold.", k)]] <- data.frame(cellID = dat.this.fold.th$cellID,
                                                                     ghat = ghat.this.fold,
                                                                     ellhat = ellhat.this.fold)
    }
    
    
    
    #### Step 1-2. Trimming 
    
    if (trim.use) {
      ghat.this.fold.trimmed <- ghat.this.fold[ghat.this.fold < 0.99 & ghat.this.fold > 0.01]
      ellhat.this.fold.trimmed <- ellhat.this.fold[ghat.this.fold < 0.99 & ghat.this.fold > 0.01]
      
      dat.this.fold.th.trimmed <- dat.this.fold.th[ghat.this.fold < 0.99 & ghat.this.fold > 0.01, ]
    } else {
      ghat.this.fold.trimmed <- ghat.this.fold
      ellhat.this.fold.trimmed <- ellhat.this.fold
      dat.this.fold.th.trimmed <- dat.this.fold.th
    }
    
    #### Step 2. Calculating theta_k & (prep for) Step 4. Calculating \Sigma_{1k}
    
    psi1.pre.this.fold <- (dat.this.fold.th.trimmed$treat - ghat.this.fold.trimmed)/phat.this.fold/(1-ghat.this.fold.trimmed)*(dat.this.fold.th.trimmed$outcome.differenced - ellhat.this.fold.trimmed)
    
    psi1.pre.repository[[paste0("fold.", k)]] <- psi1.pre.this.fold
    
    thetak.this.fold <- mean(psi1.pre.this.fold)
    
    thetak.repository <- c(thetak.repository, thetak.this.fold)
    
    phat.repository[[paste0("fold.", k)]] <- mean(dat.this.fold.nu$treat)
    
    D.repository[[paste0("fold.", k)]] <- dat.this.fold.th.trimmed$treat
    
    trimmed.N <- trimmed.N + nrow(dat.this.fold.th.trimmed)
    
    print(paste0("DONE WITH FOLD NO. ", k, "/", n.fold.used))
    
    
  }
  
  #### Step 3. Estimating theta tilde, the final point estimate of ATT
  
  thetatilde <- mean(thetak.repository)
  #the final point estimate of ATT
  
  #### Step 4. Estimating asymptotic variance
  
  Sigma1k.repository <- c()
  
  for (k in 1:n.fold.used) {
    psi1.pre.this.fold <- psi1.pre.repository[[paste0("fold.", k)]]
    
    psi1.this.fold <- psi1.pre.this.fold - thetatilde
    
    phat.this.fold <- phat.repository[[paste0("fold.", k)]]
    
    Ghat.this.fold <- 0-thetatilde/phat.this.fold
    
    to.be.meaned <- (psi1.this.fold + Ghat.this.fold*(D.repository[[paste0("fold.", k)]] - phat.this.fold))^2
    
    Sigma1k.repository <- c(Sigma1k.repository, mean(to.be.meaned))
  }
  
  
  sigma2 <- mean(Sigma1k.repository)/trimmed.N
  
  if ("Year" %in% colnames(dat.to.use)) { # if this is a pooled estimate
    dat.to.use.w.fitted <- left_join(dat.to.use,
                                     bind_rows(g.and.ellhat.df.repository),
                                     by = c("cellID", "Year"))
  } else {
    dat.to.use.w.fitted <- left_join(dat.to.use,
                                     bind_rows(g.and.ellhat.df.repository),
                                     by = "cellID")
  }
  
  
  
  dat.to.use.w.fitted$proj <- proj.to.look.at
  
  dat.to.use.w.fitted %>% 
    relocate(proj) -> dat.to.use.w.fitted
  
  # dat.to.use.w.fitted$lon <- unique.cells.data[match(dat.to.use.w.fitted$cellID, unique.cells.data$cellID),]$lon
  # dat.to.use.w.fitted$lat <- unique.cells.data[match(dat.to.use.w.fitted$cellID, unique.cells.data$cellID),]$lat
  
  t.stat <- thetatilde/sqrt(sigma2)
  
  p.val <- 2*pt(t.stat, df = nrow(dat.to.use.w.fitted), lower.tail = t.stat <=0)
  
  final.list <- list(df = dat.to.use.w.fitted,
                     eff = thetatilde,
                     sigma2 = sigma2,
                     se = sqrt(sigma2),
                     p.val = p.val)
  
  return(final.list)
  
}
