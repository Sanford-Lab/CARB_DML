library(tidyverse);library(splines);library(tictoc);library(reshape2)
library(DoubleML);library(mlr3);library(mlr3learners);library(mlr3tuning)
library(sf);library(paletteer);library(lubridate)
library(ranger);library(scales);library(DRDID)

## ========== NOTES ==============##

# CODES FOR ACTUAL PROJECTS, USING *SPATIAL CLUSTERS*
# VERSION 2024/08/09 (Using 250m pixel)
#Changes from previous version:
#1. Use of DID framework for the analysis, which entails
#1-a. Changing the "dml.datamaker" function, to include the outcome and its changes
#2. Use of clusters

## ===============================##

#### ProjectID Set-up ######

args <- commandArgs(trailingOnly = TRUE)
print(paste0("job ", args[1], " started!"))
project.ID <- "CAFR5090" #the project.ID argument, set externally
year.set <- 1 #the year argument, set externally
###########

# Running DML =====

load("Step2_1_Fixed characteristics_NDVI250m.RData")
load("Step2.5_ProjectSpecificData_b2k_updated.RData")
load("Step2.5_NonForestAndTreatedPixels.RData")
load("Step2_BIOMASS_NDVI250m.RData")
load("Step2.5_OctantsQuadrantsANDDisaggedProjects.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load('Step2_AC_LCMS.RData')
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

projs.in.dat <- unique(proj.dats.for.prop.b2k$proj)
project.ID <- projs.in.dat[args[1]]

colnames(biomass.extracted.all) <- c(paste0("biomass_", 2000:2017), "cellID")

projs.in.dat <- unique(proj.dats.for.prop.b2k$proj)

projs.with.1.possible <- proj.years[which(year(proj.years$DATE.first) <= 2016),]$project.name
projs.with.2.possible <- proj.years[which(year(proj.years$DATE.first) <= 2015),]$project.name
projs.with.3.possible <- proj.years[which(year(proj.years$DATE.first) <= 2014),]$project.name
projs.with.4.possible <- proj.years[which(year(proj.years$DATE.first) <= 2013),]$project.name

## Covariate sets =====

evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]

cov.wo.clm.forminus <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         "nlcd", 'fownership', "distance.to.road",
                         c(evi.cov.4, evi.cov.5),
                         'biomass_tminus4', 'biomass_tminus5', "forest.group")

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")


## DML runner function ======

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

dml.datamaker.for.didclust <- function(proj.to.look.at, year.incr,
                                       proj.dat.to.use, cov.to.use, outcome.to.use,
                                       nonforest.exclude = F, exclude.treat.pixels = F) {
  
  year.of.treat <- year(proj.years[proj.years$project.name == proj.to.look.at, ]$DATE.first)
  #point of treatment
  
  dat.to.use <- proj.dat.to.use[proj.dat.to.use$proj==proj.to.look.at,
                                which(colnames(proj.dat.to.use) %in% c("cellID", "type", "YEAR", "treat",
                                                                       "lon", "lat",
                                                                       cov.to.use, "biomass_0") | grepl("cluster|biomass_tminus", colnames(proj.dat.to.use)))]
  #add all cluster columns, as you do not know which one of these would be useful
  
  dat.to.use$forest.or.not <- !dat.to.use$cellID %in% non.forest.pixels 
  #TRUE: is a forest
  
  #dat.to.use[[cluster.cols]] <- proj.dat.to.use[proj.dat.to.use$proj==proj.to.look.at,][[cluster.cols]]
  
  
  if (nonforest.exclude) {
    dat.to.use <- subset(dat.to.use, forest.or.not)
  }
  
  cellids.included <- unique(dat.to.use$cellID)
  
  biomass.year.of.att.est <- select(subset(biomass.extracted.all, cellID %in% cellids.included), 
                                    c("cellID", paste0("biomass_", year.of.treat + year.incr)))
  # 
  # #in the baseline year
  # biomass.year.of.att.est <- as.data.frame(biomass.extracted.data.table[cellID %in% cellids.included & YEAR==(year.of.treat+year.incr),c(1,3)])
  
  colnames(biomass.year.of.att.est) <- c("cellID", "biomass_OUTCOME")
  
  biomass.year.of.att.est$cellID <- as.character(biomass.year.of.att.est$cellID)
  
  dat.to.use <- left_join(dat.to.use, biomass.year.of.att.est, by ="cellID")
  
  biomass.thisproj <- subset(biomass.extracted.all, cellID %in% cellids.included)
  
  biomass.thisproj %>% 
    pivot_longer(names_to = "Year", cols = starts_with("biomass_"), names_prefix = "", values_to = "biomass") %>% 
    mutate(Year = gsub("biomass_", "", Year)) %>% 
    group_by(cellID) %>% 
    mutate(biomass.delta = biomass - lag(biomass)) %>% 
    ungroup() %>% 
    select(c("cellID", "Year", "biomass.delta")) %>% 
    pivot_wider(id_cols = cellID, names_from = Year,values_from = biomass.delta,
                names_prefix = "biomass_delta.")-> biomass.thisproj
  
  dat.to.use$biomass_tminus5.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat-5)]]
  dat.to.use$biomass_tminus4.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat-4)]]
  dat.to.use$biomass_tminus3.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat-3)]]
  dat.to.use$biomass_tminus2.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat-2)]]
  dat.to.use$biomass_tminus1.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat-1)]]
  dat.to.use$biomass_0.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat)]]
  dat.to.use$biomass_OUTCOME.delta <- biomass.thisproj[match(dat.to.use$cellID, biomass.thisproj$cellID),][[paste0("biomass_delta.", year.of.treat+year.incr)]]
  
  if (year.of.treat == 2004) {
    dat.to.use$biomass_tminus5.delta <- 0
    dat.to.use$biomass_tminus4.delta <- 0
  }
  
  if (year.of.treat == 2005) {
    dat.to.use$biomass_tminus5.delta <- 0
  }
  
  if (outcome.to.use %in% c("conversion.cumul.ratio", "conversion.thisyear.ratio")) {
    
    if (year.incr > 0) {
      lcms.ratio.dat.for.this.cellids <- cells.in.refor.and.ac.projs.unique[cells.in.refor.and.ac.projs.unique$cellID %in% cellids.included,]
      
      lcms.ratio.dat.for.this.cellids <- lcms.ratio.dat.for.this.cellids[,grepl("cellID|lcms_",colnames(lcms.ratio.dat.for.this.cellids))]
      
      lcms.years <- as.integer(gsub("lcms_ratio_", "", colnames(lcms.ratio.dat.for.this.cellids)[-1]))
      
      cols.to.include.from1toyearincr <- lcms.years <= (year.of.treat + year.incr) & lcms.years >= year.of.treat+1
      #for example: if the year.of.treat was 2005, and year.incr = 1, then include all columns from 2005 to 2006
      
      cols.to.include.thisyear <- lcms.years == (year.of.treat + year.incr)
      #for example: if the year.of.treat was 2005, and year.incr = 1, then include column from 2006
      
      lcms.ratio.dat.for.this.cellids.from1toyearincr <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.from1toyearincr))]
      lcms.ratio.dat.for.thisyear <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.thisyear))]
      
      if (ncol(lcms.ratio.dat.for.this.cellids.from1toyearincr)>2) {
        conversion.ratio.cumul.fromlcms <- unname(rowSums(lcms.ratio.dat.for.this.cellids.from1toyearincr[, -1], na.rm = T))
      } else if (ncol(lcms.ratio.dat.for.this.cellids.from1toyearincr)==2) {
        conversion.ratio.cumul.fromlcms <- unname(unlist(lcms.ratio.dat.for.this.cellids.from1toyearincr[, -1]))
      }
      #whether the pixel has experienced conversion at least once ever since year.of.treat+1 to year.of.treat + year.incr
      
      conversion.ratio.thisyear <- lcms.ratio.dat.for.thisyear[[2]]
      #whether the pixel experienced conversion IN THIS YEAR
      
      lcms.df.to.fuse <- data.frame(cellID = lcms.ratio.dat.for.this.cellids.from1toyearincr$cellID,
                                    conversion.cumul.ratio = as.numeric(conversion.ratio.cumul.fromlcms),
                                    conversion.thisyear.ratio = as.numeric(conversion.ratio.thisyear))
      
      lcms.df.to.fuse$cellID <- as.character(lcms.df.to.fuse$cellID)
      
      dat.to.use <- left_join(dat.to.use, lcms.df.to.fuse, by = "cellID")
    }
    
    if (year.incr < 0) {
      #for year.incr < 0, do not calculate *cumulative* change, but changes in each year
      
      lcms.ratio.dat.for.this.cellids <- cells.in.refor.and.ac.projs.unique[cells.in.refor.and.ac.projs.unique$cellID %in% cellids.included,]
      
      lcms.ratio.dat.for.this.cellids <- lcms.ratio.dat.for.this.cellids[,grepl("cellID|lcms_",colnames(lcms.ratio.dat.for.this.cellids))]
      
      lcms.years <- as.integer(gsub("lcms_ratio_", "", colnames(lcms.ratio.dat.for.this.cellids)[-1]))
      
      cols.to.include.thisyear <- lcms.years == (year.of.treat + year.incr)
      #for example: if the year.of.treat was 2005, and year.incr = -3, then include column from 2002
      
      lcms.ratio.dat.for.thisyear <- lcms.ratio.dat.for.this.cellids[,c(1, 1+which(cols.to.include.thisyear))]
      
      #whether the pixel has experienced conversion at least once ever since year.of.treat+1 to year.of.treat + year.incr
      
      conversion.ratio.thisyear <- lcms.ratio.dat.for.thisyear[[2]]
      #whether the pixel experienced conversion IN THIS YEAR
      
      lcms.df.to.fuse <- data.frame(cellID = lcms.ratio.dat.for.thisyear$cellID,
                                    conversion.cumul.ratio = as.numeric(conversion.ratio.thisyear),
                                    conversion.thisyear.ratio = as.numeric(conversion.ratio.thisyear))
      
      lcms.df.to.fuse$cellID <- as.character(lcms.df.to.fuse$cellID)
      
      dat.to.use <- left_join(dat.to.use, lcms.df.to.fuse, by = "cellID")
    }
  }
  
  
  dat.to.use$treat <- as.integer(as.logical(dat.to.use$treat))
  
  if (exclude.treat.pixels) {
    control.but.treated <- which(dat.to.use$type=="control" & dat.to.use$cellID %in% treated.pixels)
    
    if (length(control.but.treated)!=0) {
      dat.to.use <- dat.to.use[-control.but.treated,]
    }
    
  }
  
  if (proj.to.look.at =="CAFR0001" & any(grepl(paste(evi.cov.4, collapse = "|"), cov.to.use))) {
    dat.to.use[,which(colnames(dat.to.use) %in% evi.cov.5)] <- 0
    
    dat.to.use$biomass_tminus5 <- 0
    
    #CAFR0001 is the only project with NA ndvi.5 data (just make this data completely irrelevant)
  }
  
  dat.to.use <- dat.to.use[complete.cases(dat.to.use),]
  
  for (fc in factor.covs) {
    if (fc %in% colnames(dat.to.use)) {
      dat.to.use[[fc]] <- as.factor(dat.to.use[[fc]])
    }
  }
  
  return(dat.to.use)
  
}

dml.datamaker.for.didclust.list <- function(dat.to.use, dupl.drop = T, n.fold = 3, cluster.use = F, cluster.cols = NA) {
  
  if (dupl.drop) {
    dat.to.use <- dat.to.use[!duplicated(dat.to.use$cellID),]
    #drop cellID duplicates
    #normally do this, except for pooled estimations (where we are pooling over multiple years)
  }
  
  dat.to.use.list <- list()
  
  if (cluster.use) {
    #cluster-based cross fitting
    
    cols.to.drop <- colnames(dat.to.use)[grepl('cluster', colnames(dat.to.use))]
    
    cols.to.drop <- setdiff(cols.to.drop, cluster.cols)
    
    dat.to.use <- select(dat.to.use, -cols.to.drop)
    
    crosstable.cluster.v.treat <- table(dat.to.use[[cluster.cols]], dat.to.use[["treat"]])
    
    clusters.w.treat <- rownames(crosstable.cluster.v.treat)[which(crosstable.cluster.v.treat[,2]!=0)]
    #clusters with at least one treated pixel
    clusters.control <- setdiff(rownames(crosstable.cluster.v.treat), clusters.w.treat)
    
    n.fold.to.be.used <- if (n.fold > length(clusters.w.treat)) {
      length(clusters.w.treat)
    } else {
      n.fold
    }
    
    sample.treat.clusters <- split(sample(clusters.w.treat, replace =F), 1:n.fold.to.be.used)
    sample.control.clusters <- split(sample(clusters.control, replace = F), 1:n.fold.to.be.used)
    
    for (k in 1:n.fold.to.be.used) {
      dat.to.use.list[[paste0("fold.", k)]] <- list()
      
      dat.to.use.list[[paste0("fold.", k)]]$nuisance.fit.data <- rbind(dat.to.use[which(dat.to.use[[cluster.cols]] %in% setdiff(clusters.w.treat, sample.treat.clusters[[k]])),],
                                                                       dat.to.use[which(dat.to.use[[cluster.cols]] %in% setdiff(clusters.control, sample.control.clusters[[k]])),])
      #data for fitting nuisance parameters
      
      dat.to.use.list[[paste0("fold.", k)]]$theta.fit.data <- rbind(dat.to.use[which(dat.to.use[[cluster.cols]] %in% sample.treat.clusters[[k]]),],
                                                                    dat.to.use[which(dat.to.use[[cluster.cols]] %in% sample.control.clusters[[k]]),])
    }
    
  } else {
    #random cross fitting
    
    dat.to.use <- select(dat.to.use, !contains("cluster"))
    #drop all "cluster" columns
    
    dat.to.use.treat.rows <- which(dat.to.use$treat ==1)
    dat.to.use.control.rows <- which(dat.to.use$treat==0)
    
    sample.treat <- split(sample(dat.to.use.treat.rows, replace = F), 1:n.fold)
    sample.control <- split(sample(dat.to.use.control.rows, replace = F), 1:n.fold)
    
    for (k in 1:n.fold) {
      dat.to.use.list[[paste0("fold.", k)]] <- list()
      
      
      dat.to.use.list[[paste0("fold.", k)]]$nuisance.fit.data <- rbind(dat.to.use[setdiff(dat.to.use.treat.rows, sample.treat[[k]]),],
                                                                       dat.to.use[setdiff(dat.to.use.control.rows, sample.control[[k]]),])
      #data for fitting nuisance parameters
      
      dat.to.use.list[[paste0("fold.", k)]]$theta.fit.data <- rbind(dat.to.use[sample.treat[[k]],],
                                                                    dat.to.use[sample.control[[k]],])
      #data that will be used to calculate the theta (use the nuisance parameter fitted in the previous step)
    }
  }
  
  return(dat.to.use.list)
}

dml.runner.for.didclust.baselinestr <- function(proj.to.look.at,
                                                dat.list.to.use, cov.to.use, outcome.to.use,
                                                baseline.string, 
                                                diff.with.baseline = T, #difference the tminus_5 with baseline results?
                                                diff.with.baseline.vars = c("biomass_tminus5"),
                                                g.learner = c("classif.rf", "classif.logit", "regr.ols"), 
                                                m.learner = c("regr.rf", "regr.ols"),
                                                trim.use = F) {
  
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
  
  n.fold.used <- length(dat.list.to.use)
  
  #### Step 0. biomass_0 and biomass_1 differenced with biomass_start
  
  #a step to convert the biomass_0 and biomass_1 in *differences*
  #this step may be removed for non-DID methods
  
  if (diff.with.baseline) {
    for (k in 1:n.fold.used) {
      
      for (var in diff.with.baseline.vars) {
        dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[var]] <-(dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[var]]- dat.list.to.use[[paste0("fold.", k)]]$nuisance.fit.data[[baseline.string]])
        
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


# Construct data sets ======

### For biomass ======
dml.dat.forest.didclust <- list()

error.projs.dataset <- c()

for (i in 1:length(projs.in.dat)) {
  
  tryCatch( {
    tic()
    
    proj.to.inspect <- projs.in.dat[i]
    
    print(paste0("STARTING WITH PROJ.", proj.to.inspect, "!!!"))
    
    if (proj.to.inspect %in% red.df$redundant.proj) {
      print(paste0("SKIPPING THIS PROJECT BECAUSE IT'S REDUNADNT!!"))
      
      next
    }
    
    dml.dat.forest.didclust[[proj.to.inspect]] <- list()
    
    years.for.this.proj <- 2017-year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.first)
    # the maximum number of years available for this project
    # for example, if treated in year 2005, then has 2006 ... 2017 --> 12 years to compare
    
    
    if (!proj.to.inspect %in% projs.with.1.possible) { 
      next
    } else if (proj.to.inspect %in% projs.with.1.possible) {
      
      print(paste0("YEARS -3, -2, -1, 0, AND 1..."))
      
      for (y in -3:-1) {
        dml.dat <- dml.datamaker.for.didclust(proj.to.look.at = proj.to.inspect,
                                              year.incr = y,
                                              proj.dat.to.use = proj.dats.for.prop.b2k,
                                              cov.to.use = cov.wo.clm.forminus,
                                              outcome.to.use = "biomass_OUTCOME",
                                              nonforest.exclude = T,
                                              exclude.treat.pixels = T)
        
        dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y))]] <- dml.dat
      }
      
      dml.dat <- dml.datamaker.for.didclust(proj.to.look.at = proj.to.inspect,
                                            year.incr = 0,
                                            proj.dat.to.use = proj.dats.for.prop.b2k,
                                            cov.to.use = cov.wo.clm.forminus,
                                            outcome.to.use = "biomass_OUTCOME",
                                            nonforest.exclude = T,
                                            exclude.treat.pixels = T)
      
      dml.dat.forest.didclust[[proj.to.inspect]]$year0 <- dml.dat
      
      dml.dat <- dml.datamaker.for.didclust(proj.to.look.at = proj.to.inspect,
                                            year.incr = 1,
                                            proj.dat.to.use = proj.dats.for.prop.b2k,
                                            cov.to.use = cov.wo.clm.forminus,
                                            outcome.to.use = "biomass_OUTCOME",
                                            nonforest.exclude = T,
                                            exclude.treat.pixels = T)
      
      dml.dat.forest.didclust[[proj.to.inspect]]$year1 <- dml.dat
      
    }
    
    if (proj.to.inspect %in% projs.with.2.possible) {
      
      for (y in 2:years.for.this.proj) {
        print(paste0("YEAR ", y, "..."))
        
        dml.dat <- dml.datamaker.for.didclust(proj.to.look.at = proj.to.inspect,
                                              year.incr = y,
                                              proj.dat.to.use = proj.dats.for.prop.b2k,
                                              cov.to.use = cov.wo.clm.forminus,
                                              outcome.to.use = "biomass_OUTCOME",
                                              nonforest.exclude = T,
                                              exclude.treat.pixels = T)
        
        dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year", y)]] <- dml.dat
        
        
      }
    }
    
    toc()
    tic()
    
    print(paste0("DONE WITH PROJECT NO. ", i))
    
  }, error = function(e){
    error.projs.dataset <<- c(error.projs.dataset, i)
    
    print(paste0("ERROR IN PROJ NO. ", i, "!! SKIPPING..."))
  }   )
  
  if (i %%5 ==0 | i == length(projs.in.dat)) {
    save(list = c("dml.dat.forest.didclust",
                  "error.projs.dataset"),
         file = paste0("Step3_DML_DIDClust_DataBiomassTminus3_", i, ".RData"))
  }
  
}

# Fit DMLDiD clustered ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")

load("Step2_1_Fixed characteristics_NDVI250m.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step2.5_NonForestAndTreatedPixels.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

projs.in.dat <- unique(proj.dats.for.prop.b2k$proj)

projs.with.1.possible <- proj.years[which(year(proj.years$DATE.first) <= 2016),]$project.name
projs.with.2.possible <- proj.years[which(year(proj.years$DATE.first) <= 2015),]$project.name
projs.with.3.possible <- proj.years[which(year(proj.years$DATE.first) <= 2014),]$project.name
projs.with.4.possible <- proj.years[which(year(proj.years$DATE.first) <= 2013),]$project.name

evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]

cov.wo.clm.forminus <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         "forest.group",
                         "nlcd", 'fownership', "distance.to.road",
                         c(evi.cov.5), 
                         "biomass_tminus5"
                         )

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

dml.fitted.res.trim.forest.didclust <- list()
dml.dat.list.forest.didclust <- list()

error.projs.fit <- c()

for (i in 50:50) {
  
  tryCatch( {tic()
    
    proj.to.inspect <- projs.in.dat[i]
    
    print(paste0("STARTING WITH PROJ.", proj.to.inspect, "!!!"))
    
    if (proj.to.inspect %in% red.df$redundant.proj) {
      print(paste0("SKIPPING THIS PROJECT BECAUSE IT'S REDUNADNT!!"))
      
      next
    }
    
    dml.fitted.res.trim.forest.didclust <- list()
    dml.dat.list.forest.didclust <- list()
    
    dml.fitted.res.trim.forest.didclust[[proj.to.inspect]] <- list()
    
    proj.fitted.models.trim.forest.didclust <- list()
    
    years.for.this.proj <- 2017-year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.first)
    # the maximum number of years available for this project
    # for example, if treated in year 2005, then has 2006 ... 2017 --> 12 years to compare
    
    final.year.of.this.proj <- if (proj.to.inspect %in% red.df$original.proj) {
      later.project <- red.df[red.df$original.proj==proj.to.inspect,]$redundant.proj
      
      min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
          2017)
      
    } else {
       min(year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.last),
                                     2017)
    }
    
    years.to.final.year <- final.year.of.this.proj-year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.first)
    
    
    if (!proj.to.inspect %in% projs.with.1.possible) { 
      next
    } else if (proj.to.inspect %in% projs.with.1.possible) {
      
      print("YEARS -3, -2, -1...")
      
      for (y in -3:-1) {
        dat.to.be.used <- dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y))]]
        
        dat.list.to.use.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                     cluster.use = T, cluster.cols = 'cluster.1km')
        
        dat.list.to.use.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                      cluster.use = T, cluster.cols = 'cluster.25km')
        
        dat.list.to.use.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                     cluster.use = T, cluster.cols = 'cluster.5km')
        
        dat.list.to.use.noclust <-dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use =F)
        
        # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y), ".clust1km")]] <- dat.list.to.use.clust1km
        dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y), ".clust25km")]] <- dat.list.to.use.clust25km
        # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y), ".clust5km")]] <- dat.list.to.use.clust5km
        # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year.minus", abs(y), ".noclust")]] <- dat.list.to.use.noclust
        
        proj.fitted.models.trim.forest.didclust[[paste0("year.minus", abs(y), ".noclust")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                              dat.list.to.use = dat.list.to.use.noclust, 
                                                                                              cov.to.use = cov.wo.clm.forminus, 
                                                                                              outcome.to.use = "biomass_OUTCOME", 
                                                                                              baseline.string = "biomass_tminus5",
                                                                                              diff.with.baseline = F, 
                                                                                              g.learner = "classif.rf",
                                                                                              m.learner = "regr.rf",
                                                                                              trim.use = T
        )
        # 
        # proj.fitted.models.trim.forest.didclust[[paste0("year.minus", abs(y), ".clust1km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
        #                                                                                        dat.list.to.use = dat.list.to.use.clust1km, 
        #                                                                                        cov.to.use = cov.wo.clm.forminus, 
        #                                                                                        outcome.to.use = "biomass_OUTCOME", 
        #                                                                                        baseline.string = "biomass_tminus5",
        #                                                                                        g.learner = "classif.rf",
        #                                                                                        m.learner = "regr.rf",
        #                                                                                        trim.use = T
        # )
        
        proj.fitted.models.trim.forest.didclust[[paste0("year.minus", abs(y), ".clust25km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                dat.list.to.use = dat.list.to.use.clust25km, 
                                                                                                cov.to.use = cov.wo.clm.forminus, 
                                                                                                outcome.to.use = "biomass_OUTCOME", 
                                                                                                baseline.string = "biomass_tminus5",
                                                                                                diff.with.baseline = F, 
                                                                                                g.learner = "classif.rf",
                                                                                                m.learner = "regr.rf",
                                                                                                trim.use = T
        )
        # 
        # proj.fitted.models.trim.forest.didclust[[paste0("year.minus", abs(y), ".clust5km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
        #                                                                                        dat.list.to.use = dat.list.to.use.clust5km, 
        #                                                                                        cov.to.use = cov.wo.clm.forminus, 
        #                                                                                        outcome.to.use = "biomass_OUTCOME", 
        #                                                                                        baseline.string = "biomass_tminus5",
        #                                                                                        g.learner = "classif.rf",
        #                                                                                        m.learner = "regr.rf",
        #                                                                                        trim.use = T
        # )
        
        
      }
      
      print("YEAR 0...")
      
      dat.to.be.used <- dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year", 0)]]
      
      dat.list.to.use.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.1km')
      
      dat.list.to.use.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                    cluster.use = T, cluster.cols = 'cluster.25km')
      
      dat.list.to.use.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.5km')
      
      dat.list.to.use.noclust <-dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                 cluster.use =F)
      
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year0", ".clust1km")]] <- dat.list.to.use.clust1km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year0", ".clust25km")]] <- dat.list.to.use.clust25km
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year0", ".clust5km")]] <- dat.list.to.use.clust5km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("year0", ".noclust")]] <- dat.list.to.use.noclust
      
      proj.fitted.models.trim.forest.didclust[[paste0("year0", ".noclust")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                                 dat.list.to.use = dat.list.to.use.noclust, 
                                                                                                                                 cov.to.use = cov.wo.clm.forminus, 
                                                                                                                                 outcome.to.use = "biomass_OUTCOME", 
                                                                                                                                 baseline.string = "biomass_tminus5",
                                                                                                                    diff.with.baseline = F, 
                                                                                                                                 g.learner = "classif.rf",
                                                                                                                                 m.learner = "regr.rf",
                                                                                                                                 trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("year0", ".clust1km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                                 dat.list.to.use = dat.list.to.use.clust1km, 
      #                                                                                                                 cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                                 outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                                 baseline.string = "biomass_tminus5",
      #                                                                                                                 g.learner = "classif.rf",
      #                                                                                                                 m.learner = "regr.rf",
      #                                                                                                                 trim.use = T
      # )
      
      proj.fitted.models.trim.forest.didclust[[paste0("year0", ".clust25km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                       dat.list.to.use = dat.list.to.use.clust25km, 
                                                                                                                       cov.to.use = cov.wo.clm.forminus, 
                                                                                                                       outcome.to.use = "biomass_OUTCOME", 
                                                                                                                       baseline.string = "biomass_tminus5",
                                                                                                                      diff.with.baseline = F, 
                                                                                                                       g.learner = "classif.rf",
                                                                                                                       m.learner = "regr.rf",
                                                                                                                       trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("year0", ".clust5km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                                 dat.list.to.use = dat.list.to.use.clust5km, 
      #                                                                                                                 cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                                 outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                                 baseline.string = "biomass_tminus5",
      #                                                                                                                 g.learner = "classif.rf",
      #                                                                                                                 m.learner = "regr.rf",
      #                                                                                                                 trim.use = T
      # )
      
      print("YEAR FINAL...") 
      
      
      dat.to.be.used <- dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year", years.to.final.year)]]
      
      dat.list.to.use.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.1km')
      
      dat.list.to.use.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                    cluster.use = T, cluster.cols = 'cluster.25km')
      
      dat.list.to.use.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.5km')
      
      dat.list.to.use.noclust <-dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                 cluster.use =F)
      
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust1km")]] <- dat.list.to.use.clust1km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust25km")]] <- dat.list.to.use.clust25km
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust5km")]] <- dat.list.to.use.clust5km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".noclust")]] <- dat.list.to.use.noclust
      
      proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".noclust")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                    dat.list.to.use = dat.list.to.use.noclust, 
                                                                                                                    cov.to.use = cov.wo.clm.forminus, 
                                                                                                                    outcome.to.use = "biomass_OUTCOME", 
                                                                                                                    baseline.string = "biomass_tminus5",
                                                                                                                    diff.with.baseline = F,
                                                                                                                    g.learner = "classif.rf",
                                                                                                                    m.learner = "regr.rf",
                                                                                                                    trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust1km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust1km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5", 
      #                                                                                                   diff.with.baseline = F,
      #                                                                                                    g.learner = "classif.rf", 
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
      
      proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust25km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                          dat.list.to.use = dat.list.to.use.clust25km, 
                                                                                                          cov.to.use = cov.wo.clm.forminus, 
                                                                                                          outcome.to.use = "biomass_OUTCOME", 
                                                                                                          baseline.string = "biomass_tminus5",
                                                                                                          diff.with.baseline = F,
                                                                                                          g.learner = "classif.rf",
                                                                                                          m.learner = "regr.rf",
                                                                                                          trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust5km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust5km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5",
      #                                                                                                    g.learner = "classif.rf",
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
    }
    
    
    dml.fitted.res.trim.forest.didclust[[proj.to.inspect]] <- proj.fitted.models.trim.forest.didclust
    
    save(list = c("dml.fitted.res.trim.forest.didclust",
                  "error.projs.fit"),
         file = paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
    
    save(list = c("dml.dat.list.forest.didclust"),
         file = paste0("Step3_DML_DIDClustwTminus3_DatList_", i, ".RData"))
    
    toc()
    tic()
    
    print(paste0("DONE WITH PROJECT NO. ", i))
    
    #restore empty lists
    dml.fitted.res.trim.forest.didclust <- list()
    dml.dat.list.forest.didclust <- list()
    
  }, error = function(e){
    error.projs.fit <<- c(error.projs.fit, i)
    
    if (length(proj.fitted.models.trim.forest.didclust)!=0) {
      dml.fitted.res.trim.forest.didclust[[proj.to.inspect]] <<- proj.fitted.models.trim.forest.didclust
      
      save(list = c("dml.fitted.res.trim.forest.didclust",
                    "error.projs.fit"),
           file = paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
      
      save(list = c("dml.dat.list.forest.didclust"),
           file = paste0("Step3_DML_DIDClustwTminus3_DatList_", i, ".RData"))
      
      dml.fitted.res.trim.forest.didclust <- list()
      dml.dat.list.forest.didclust <- list()
    }
    
    
    print(paste0("ERROR IN PROJ NO. ", i, "!! SKIPPING..."))
  }   )
}
#Error with sample.fraction: Only happens when the project area is so tiny that there is just one fold

# Fit DMLDiD clustered, 0 versus year final & 0 versus 5 years ago =====

dml.fitted.res.trim.forest.didclust.v0 <- list()
# dml.dat.list.forest.didclust <- list()

error.projs.fit <- c()

for (i in 50:50) {
  
  tryCatch( {tic()
    
    proj.to.inspect <- projs.in.dat[i]
    
    print(paste0("STARTING WITH PROJ.", proj.to.inspect, "!!!"))
    
    if (proj.to.inspect %in% red.df$redundant.proj) {
      print(paste0("SKIPPING THIS PROJECT BECAUSE IT'S REDUNADNT!!"))
      
      next
    }
    
    dml.fitted.res.trim.forest.didclust.v0 <- list()
    dml.dat.list.forest.didclust <- list()
    
    dml.fitted.res.trim.forest.didclust.v0[[proj.to.inspect]] <- list()
    
    proj.fitted.models.trim.forest.didclust <- list()
    
    years.for.this.proj <- 2017-year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.first)
    # the maximum number of years available for this project
    # for example, if treated in year 2005, then has 2006 ... 2017 --> 12 years to compare
    
    final.year.of.this.proj <- if (proj.to.inspect %in% red.df$original.proj) {
      later.project <- red.df[red.df$original.proj==proj.to.inspect,]$redundant.proj
      
      min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
          2017)
      
    } else {
      min(year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.last),
          2017)
    }
    
    years.to.final.year <- final.year.of.this.proj-year(proj.years[proj.years$project.name == proj.to.inspect,]$DATE.first)
    
    
    if (!proj.to.inspect %in% projs.with.1.possible) { 
      next
    } else if (proj.to.inspect %in% projs.with.1.possible) {
      
      print("YEAR FINAL...") 
    
      dat.to.be.used <- dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year", years.to.final.year)]]
      
      dat.list.to.use.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.1km')
      
      dat.list.to.use.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                    cluster.use = T, cluster.cols = 'cluster.25km')
      
      dat.list.to.use.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.5km')
      
      dat.list.to.use.noclust <-dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                 cluster.use =F)
      
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust1km")]] <- dat.list.to.use.clust1km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust25km")]] <- dat.list.to.use.clust25km
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust5km")]] <- dat.list.to.use.clust5km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".noclust")]] <- dat.list.to.use.noclust
      
      proj.fitted.models.trim.forest.didclust[[paste0("yearfinalv0", ".noclust")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                        dat.list.to.use = dat.list.to.use.noclust, 
                                                                                                                        cov.to.use = cov.wo.clm.forminus, 
                                                                                                                        outcome.to.use = "biomass_OUTCOME", 
                                                                                                                        baseline.string = "biomass_0",
                                                                                                                        diff.with.baseline = F,
                                                                                                                        g.learner = "classif.rf",
                                                                                                                        m.learner = "regr.rf",
                                                                                                                        trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust1km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust1km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5", 
      #                                                                                                   diff.with.baseline = F,
      #                                                                                                    g.learner = "classif.rf", 
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
      
      proj.fitted.models.trim.forest.didclust[[paste0("yearfinalv0", ".clust25km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                          dat.list.to.use = dat.list.to.use.clust25km, 
                                                                                                                          cov.to.use = cov.wo.clm.forminus, 
                                                                                                                          outcome.to.use = "biomass_OUTCOME", 
                                                                                                                          baseline.string = "biomass_0",
                                                                                                                          diff.with.baseline = F,
                                                                                                                          g.learner = "classif.rf",
                                                                                                                          m.learner = "regr.rf",
                                                                                                                          trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust5km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust5km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5",
      #                                                                                                    g.learner = "classif.rf",
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
      
      print("YEAR 0 VERSUS 5...") 
      
      dat.to.be.used <- dml.dat.forest.didclust[[proj.to.inspect]][[paste0("year", 0)]]
      
      dat.list.to.use.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.1km')
      
      dat.list.to.use.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                    cluster.use = T, cluster.cols = 'cluster.25km')
      
      dat.list.to.use.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                   cluster.use = T, cluster.cols = 'cluster.5km')
      
      dat.list.to.use.noclust <-dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used,
                                                                 cluster.use =F)
      
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust1km")]] <- dat.list.to.use.clust1km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust25km")]] <- dat.list.to.use.clust25km
      # dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".clust5km")]] <- dat.list.to.use.clust5km
      dml.dat.list.forest.didclust[[proj.to.inspect]][[paste0("yearfinal", ".noclust")]] <- dat.list.to.use.noclust
      
      proj.fitted.models.trim.forest.didclust[[paste0("year0v5", ".noclust")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                        dat.list.to.use = dat.list.to.use.noclust, 
                                                                                                                        cov.to.use = cov.wo.clm.forminus, 
                                                                                                                        outcome.to.use = "biomass_0", 
                                                                                                                        baseline.string = "biomass_tminus5",
                                                                                                                        diff.with.baseline = F,
                                                                                                                        g.learner = "classif.rf",
                                                                                                                        m.learner = "regr.rf",
                                                                                                                        trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust1km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust1km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5", 
      #                                                                                                   diff.with.baseline = F,
      #                                                                                                    g.learner = "classif.rf", 
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
      
      proj.fitted.models.trim.forest.didclust[[paste0("year0v5", ".clust25km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                                          dat.list.to.use = dat.list.to.use.clust25km, 
                                                                                                                          cov.to.use = cov.wo.clm.forminus, 
                                                                                                                        outcome.to.use = "biomass_0", 
                                                                                                                        baseline.string = "biomass_tminus5",
                                                                                                                          diff.with.baseline = F,
                                                                                                                          g.learner = "classif.rf",
                                                                                                                          m.learner = "regr.rf",
                                                                                                                          trim.use = T
      )
      
      # proj.fitted.models.trim.forest.didclust[[paste0("yearfinal", ".clust5km")]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
      #                                                                                                    dat.list.to.use = dat.list.to.use.clust5km, 
      #                                                                                                    cov.to.use = cov.wo.clm.forminus, 
      #                                                                                                    outcome.to.use = "biomass_OUTCOME", 
      #                                                                                                    baseline.string = "biomass_tminus5",
      #                                                                                                    g.learner = "classif.rf",
      #                                                                                                    m.learner = "regr.rf",
      #                                                                                                    trim.use = T
      # )
    }
    
    
    dml.fitted.res.trim.forest.didclust.v0[[proj.to.inspect]] <- proj.fitted.models.trim.forest.didclust
    
    save(list = c("dml.fitted.res.trim.forest.didclust.v0",
                  "error.projs.fit"),
         file = paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData"))
    
    # save(list = c("dml.dat.list.forest.didclust"),
    #      file = paste0("Step3_DML_DIDClustwTminus3_DatList_", i, ".RData"))
    
    toc()
    tic()
    
    print(paste0("DONE WITH PROJECT NO. ", i))
    
    #restore empty lists
    dml.fitted.res.trim.forest.didclust <- list()
    dml.dat.list.forest.didclust <- list()
    
  }, error = function(e){
    error.projs.fit <<- c(error.projs.fit, i)
    
    if (length(proj.fitted.models.trim.forest.didclust)!=0) {
      dml.fitted.res.trim.forest.didclust[[proj.to.inspect]] <<- proj.fitted.models.trim.forest.didclust
      
      save(list = c("dml.fitted.res.trim.forest.didclust",
                    "error.projs.fit"),
           file = paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData"))
      
      # save(list = c("dml.dat.list.forest.didclust"),
      #      file = paste0("Step3_DML_DIDClustwTminus3_DatList_", i, ".RData"))
      
      dml.fitted.res.trim.forest.didclust <- list()
      dml.dat.list.forest.didclust <- list()
    }
    
    
    print(paste0("ERROR IN PROJ NO. ", i, "!! SKIPPING..."))
  }   )
}
#Error with sample.fraction: Only happens when the project area is so tiny that there is just one fold



# Evaluate results (new version 24/09/18, with multiple baseline years) ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

dml.fitted.result.summarized <- data.frame(proj.id = character(),
                                           fownership.type = character(),
                                           forest.type = character(),
                                           # year.calendar = rep(proj.start.year, 4) + year.integer,
                                           year.start = integer(),
                                           year.final = integer(),
                                           year.str = character(),
                                           # year.str.plusminus = rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                           #                                 as.character(year.integer)), 4),
                                           clust.type = character(),
                                           avg.biomass = numeric(),
                                           sd.biomass = numeric(),
                                           eff = numeric(),
                                           se = numeric(),
                                           trimmed.N = integer(),
                                           treated.N = integer(),
                                           cluster.count.total = integer(),
                                           cluster.count = integer(),
                                           rsquared.outcome = numeric(),
                                           mean.prob.control = numeric(),
                                           mean.prob.treat = numeric(),
                                           mean.prob.control.logit = numeric(),
                                           mean.prob.treat.logit = numeric(),
                                           median.prob.control = numeric(),
                                           median.prob.treat = numeric(),
                                           median.prob.control.logit = numeric(),
                                           median.prob.treat.logit = numeric(),
                                           cor.between.
                                           )

projs.in.dat <- unique(proj.dats.for.prop.b2k$project.ID)

no.dat.projs <- c()

Mode.giver <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for (i in 1:length(projs.in.dat)) {
  
  projs.data.there <- paste0("Step3_DML_DIDClustwTminus3_", i, ".RData") %in% list.files()
  
  if (!projs.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinal", paste0("year.minus", 3:1), "year0")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    avg.biomass <- mean(subset(noclust.dat$df, treat ==1)$biomass_0, na.rm = T)
    sd.biomass <- sd(noclust.dat$df$biomass_0, na.rm = T) #Standard deviation taken over the entire sample
    
    if (!is.null(clust25km.dat)) {
      glm.fitted.25 <- glm(treat~.,
                           select(clust25km.dat$df, c("treat", cov.wo.clm.forminus)), family = binomial())
      
      clust25km.dat$df$logit.p <- glm.fitted.25$fitted.values
    }
    
    if (!is.null(noclust.dat)) {
      
      glm.fitted.noclust <- glm(treat~.,
                                select(noclust.dat$df, c("treat", cov.wo.clm.forminus)), family = binomial())
      
      
      noclust.dat$df$logit.p <- glm.fitted.noclust$fitted.values
    }
    
    
    
    
    dml.fitted.result.summarized <- rbind(dml.fitted.result.summarized,
                                          data.frame(proj.id = rep(proj.here, 2),
                                                     fownership.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$fownership),2),
                                                     forest.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$forest.group),2),
                                                     # year.calendar = rep(proj.start.year, 4) + year.integer,
                                                     year.start = rep(proj.start.year, 2),
                                                     year.final = rep(final.year.of.this.proj, 2),
                                                     year.str = rep(str,2),
                                                     # year.str.plusminus= rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                                     #                                 as.character(year.integer)), 4),
                                                     clust.type = c("noclust", 
                                                                    "clust.25km"),
                                                     avg.biomass = rep(avg.biomass, 2),
                                                     sd.biomass = rep(sd.biomass, 2),
                                                     eff = c(ifelse(is.null(noclust.dat), NA, noclust.dat$eff),
                                                             ifelse(is.null(clust25km.dat), NA, clust25km.dat$eff)),
                                                     se = c(ifelse(is.null(noclust.dat), NA, noclust.dat$se),
                                                            ifelse(is.null(clust25km.dat), NA, clust25km.dat$se)),
                                                     trimmed.N = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                   ifelse(is.null(clust25km.dat), NA, sum(as.integer(clust25km.dat$df$ghat < 0.99 & clust25km.dat$df$ghat > 0.01), na.rm = T))),
                                                     total.N = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                 ifelse(is.null(clust25km.dat), NA, nrow(clust25km.dat$df))),
                                                     treated.N = c(ifelse(is.null(noclust.dat), NA, nrow(subset(noclust.dat$df, treat ==1))),
                                                                 ifelse(is.null(clust25km.dat), NA, nrow(subset(noclust.dat$df, treat ==1)))),
                                                     # cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                     #                   ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km))),
                                                     #                   ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km))),
                                                     #                   ifelse(is.null(clust5km.dat), NA, length(unique(clust5km.dat$df$cluster.5km)))),
                                                     cluster.count.total = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                             ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km)))),
                                                     cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                       ifelse(is.null(clust25km.dat), NA, length(unique(subset(clust25km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.25km)))),
                                                     rsquared.outcome = c(ifelse(is.null(noclust.dat), NA, cor(noclust.dat$df$biomass_OUTCOME - noclust.dat$df$biomass_tminus5,
                                                                                                               noclust.dat$df$ellhat)),
                                                                          ifelse(is.null(clust25km.dat), NA, cor(clust25km.dat$df$biomass_OUTCOME - clust25km.dat$df$biomass_tminus5,
                                                                                                                 clust25km.dat$df$ellhat))
                                                     ),
                                                     mean.prob.control = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==0)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T))),
                                                     mean.prob.treat = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==1)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==1)$ghat, na.rm = T))),
                                                     mean.prob.control.logit = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==0)$logit.p, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$logit.p, na.rm = T))),
                                                     mean.prob.treat.logit = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==1)$logit.p, na.rm = T)),
                                                                         ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==1)$logit.p, na.rm = T))),
                                                     median.prob.control = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==0)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T))),
                                                     median.prob.treat = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==1)$ghat, na.rm = T)),
                                                                         ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==1)$ghat, na.rm = T))),
                                                     median.prob.control.logit = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==0)$logit.p, na.rm = T)),
                                                                                 ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==0)$logit.p, na.rm = T))),
                                                     median.prob.treat.logit = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==1)$logit.p, na.rm = T)),
                                                                               ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==1)$logit.p, na.rm = T)))
                                                     )
    )
    
    
    
  }
  
  projs.0.data.there <- paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData") %in% list.files()
  
  if (!projs.0.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust.v0) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust.v0[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinalv0", "year0v5")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    avg.biomass <- mean(subset(noclust.dat$df, treat ==1)$biomass_0, na.rm = T)
    sd.biomass <- sd(noclust.dat$df$biomass_0, na.rm = T)
    
    if (!is.null(clust25km.dat)) {
      glm.fitted.25 <- glm(treat~.,
                           select(clust25km.dat$df, c("treat", cov.wo.clm.forminus)), family = binomial())
      
      clust25km.dat$df$logit.p <- glm.fitted.25$fitted.values
    }
    
    if (!is.null(noclust.dat)) {
      
      glm.fitted.noclust <- glm(treat~.,
                                select(noclust.dat$df, c("treat", cov.wo.clm.forminus)), family = binomial())
      
      
      noclust.dat$df$logit.p <- glm.fitted.noclust$fitted.values
    }
    
    
    dml.fitted.result.summarized <- rbind(dml.fitted.result.summarized,
                                          data.frame(proj.id = rep(proj.here, 2),
                                                     fownership.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$fownership),2),
                                                     forest.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$forest.group),2),
                                                     # year.calendar = rep(proj.start.year, 4) + year.integer,
                                                     year.start = rep(proj.start.year, 2),
                                                     year.final = rep(final.year.of.this.proj, 2),
                                                     year.str = rep(str,2),
                                                     # year.str.plusminus = rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                                     #                                 as.character(year.integer)), 4),
                                                     clust.type = c("noclust", 
                                                                    "clust.25km"),
                                                     avg.biomass = rep(avg.biomass, 2),
                                                     sd.biomass = rep(sd.biomass, 2),
                                                     eff = c(ifelse(is.null(noclust.dat), NA, noclust.dat$eff),
                                                             ifelse(is.null(clust25km.dat), NA, clust25km.dat$eff)),
                                                     se = c(ifelse(is.null(noclust.dat), NA, noclust.dat$se),
                                                            ifelse(is.null(clust25km.dat), NA, clust25km.dat$se)),
                                                     trimmed.N = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                   ifelse(is.null(clust25km.dat), NA, sum(as.integer(clust25km.dat$df$ghat < 0.99 & clust25km.dat$df$ghat > 0.01), na.rm = T))),
                                                     total.N = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                 ifelse(is.null(clust25km.dat), NA, nrow(clust25km.dat$df))),
                                                     treated.N = c(ifelse(is.null(noclust.dat), NA, nrow(subset(noclust.dat$df, treat ==1))),
                                                                   ifelse(is.null(clust25km.dat), NA, nrow(subset(noclust.dat$df, treat ==1)))),
                                                     # cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                     #                   ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km))),
                                                     #                   ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km))),
                                                     #                   ifelse(is.null(clust5km.dat), NA, length(unique(clust5km.dat$df$cluster.5km)))),
                                                     cluster.count.total = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                             ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km)))),
                                                     cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                       ifelse(is.null(clust25km.dat), NA, length(unique(subset(clust25km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.25km)))),
                                                     rsquared.outcome = c(ifelse(is.null(noclust.dat), NA, cor(noclust.dat$df$biomass_OUTCOME - noclust.dat$df$biomass_tminus5,
                                                                                                               noclust.dat$df$ellhat)),
                                                                          ifelse(is.null(clust25km.dat), NA, cor(clust25km.dat$df$biomass_OUTCOME - clust25km.dat$df$biomass_tminus5,
                                                                                                                 clust25km.dat$df$ellhat))
                                                     ),
                                                     mean.prob.control = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==0)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T))),
                                                     mean.prob.treat = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==1)$ghat, na.rm = T)),
                                                                         ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==1)$ghat, na.rm = T))),
                                                     mean.prob.control.logit = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==0)$logit.p, na.rm = T)),
                                                                                 ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$logit.p, na.rm = T))),
                                                     mean.prob.treat.logit = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==1)$logit.p, na.rm = T)),
                                                                               ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==1)$logit.p, na.rm = T))),
                                                     median.prob.control = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==0)$ghat, na.rm = T)),
                                                                             ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T))),
                                                     median.prob.treat = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==1)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==1)$ghat, na.rm = T))),
                                                     median.prob.control.logit = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==0)$logit.p, na.rm = T)),
                                                                                   ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==0)$logit.p, na.rm = T))),
                                                     median.prob.treat.logit = c(ifelse(is.null(noclust.dat), NA, median(subset(noclust.dat$df, treat==1)$logit.p, na.rm = T)),
                                                                                 ifelse(is.null(clust25km.dat), NA, median(subset(clust25km.dat$df,treat==1)$logit.p, na.rm = T)))
                                          )
    )
    
    
    
  }
  
  print(paste0("Done with proj no. ", i))
  
  rm(dml.fitted.res.trim.forest.didclust)
  
}


dml.fitted.result.summarized$se.f <- dml.fitted.result.summarized$se*sqrt(dml.fitted.result.summarized$total.N)/sqrt(dml.fitted.result.summarized$cluster.count.total)

dml.fitted.result.summarized$t.stat.f <- dml.fitted.result.summarized$eff/dml.fitted.result.summarized$se.f

pval.calc.function <- function(x,y) {2*pt(x, df = y, lower.tail = (x<=0))}

dml.fitted.result.summarized$pval.f <- mapply(pval.calc.function,
                                              dml.fitted.result.summarized$t.stat.f,
                                              dml.fitted.result.summarized$cluster.count)

table(dml.fitted.result.summarized$clust.type,
      dml.fitted.result.summarized$pval.f < 0.05,
      dml.fitted.result.summarized$year.str)

View(subset(dml.fitted.result.summarized, clust.type == "clust.25km" & pval.f < 0.05 & year.str == "yearfinal"))
#note that the ones that are insignificant have very small effect estimates

save(list = "dml.fitted.result.summarized",
     file = "Step3_DMLDiDClustered_summarized.RData")

dml.fitted.result.summarized %>% 
  mutate(mean.prob.overlap.ML = mean.prob.control + 1- mean.prob.treat,
         mean.prob.overlap.nonML = mean.prob.control.logit + 1- mean.prob.treat.logit,
         median.prob.overlap.ML = median.prob.control + 1- median.prob.treat,
         median.prob.overlap.nonML = median.prob.control.logit + 1- median.prob.treat.logit) -> dml.fitted.result.summarized
#mean.prob.control + 1-mean.prob.treat = will be 1 if perfectly balanced, 0.01 + 1 - 0.99 = 0.02 if imbalanced

t.test(subset(dml.fitted.result.summarized,clust.type == "clust.25km" & year.str == "yearfinalv0")$median.prob.overlap.ML,
       subset(dml.fitted.result.summarized,clust.type == "clust.25km" & year.str == "yearfinalv0")$median.prob.overlap.nonML)

ggplot(subset(dml.fitted.result.summarized,clust.type == "clust.25km" & year.str == "yearfinalv0")) +
  geom_point(aes(x = mean.prob.overlap.ML, y = mean.prob.overlap.nonML)) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  paletteer::scale_color_paletteer_d("palettesForR::Named", direction = -1) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0, size = 15, face = 'bold'),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, median.prob.control, mean.prob.treat) %>% 
  pivot_longer(cols = contains('mean.prob'),
               names_to = "group",
               names_prefix = "mean.prob.",
               values_to = "mean.prob") %>% 
  ggpubr::ggdensity(x = "mean.prob", color = "group", fill = "group", add = "mean")

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, median.prob.overlap.ML, median.prob.overlap.nonML) %>% 
  pivot_longer(cols = contains('median.prob.overlap'),
               names_to = "group",
               names_prefix = "median.prob.overlap.",
               values_to = "overlap") %>% 
  ggpubr::gghistogram(x = "overlap", color = "group", fill = "group", add = "mean",
                      bins = 20, alpha = .2) %>% 
  ggpar(xlab = "Propensity score overlap", legend.title = "Model used", palette = c("#0E7175FF", "#FD7901FF"))

ggsave("C:/Users/samue/Dropbox/Apps/Overleaf/Offset_Nature/Figure2_AllProjs_overlap.eps", 
       width = 10, height = 8)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, mean.prob.control, mean.prob.control.logit) %>% 
  pivot_longer(cols = contains('mean.prob'),
               names_to = "group",
               names_prefix = "mean.prob.control",
               values_to = "mean.prob") %>% 
  mutate(group = case_when(
    group == "" ~ "ML model",
    group == ".logit" ~ "Non-ML model")) %>% 
  ggpubr::gghistogram(x = "mean.prob", color = "group", fill = "group", add = "mean",
                      bins = 20, alpha = .2) %>% 
  ggpar(xlab = "Propensity score overlap", legend.title = "Model used", palette = c("#0E7175FF", "#FD7901FF"))


dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, mean.prob.treat, mean.prob.treat.logit) %>% 
  pivot_longer(cols = contains('mean.prob'),
               names_to = "group",
               names_prefix = "mean.prob.treat",
               values_to = "mean.prob") %>% 
  mutate(group = case_when(
    group == "" ~ "ML model",
    group == ".logit" ~ "Non-ML model")) %>% 
  ggpubr::gghistogram(x = "mean.prob", color = "black", fill = "group", add = "mean",
                      alpha = .2, bins = 21) %>% 
  ggpar(xlab = "Average propensity score of FCOP-participating locations", legend.title = "Model used", palette = c("#0E7175FF", "#FD7901FF"))

ggsave("C:/Users/samue/Dropbox/Apps/Overleaf/Offset_Nature/Figure2_AllProjs_pscoreTreat.eps", 
       device = cairo_ps,
       width = 10, height = 8)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, mean.prob.control, mean.prob.treat) %>% 
  pivot_longer(cols = contains('mean.prob'),
               names_to = "group",
               names_prefix = "mean.prob",
               values_to = "mean.prob") %>% 
  mutate(group = case_when(
    group == ".treat" ~ "Participating locations",
    group == ".control" ~ "Non-participating locations")) %>% 
  ggpubr::gghistogram(x = "mean.prob", color = "group", fill = "group", add = "median",
                      bins = 20, alpha = .2) %>% 
  ggpar(xlab = "Average propensity score", legend.title = "", palette = c("#0E7175FF", "#FD7901FF"))

ggsave("C:/Users/samue/Dropbox/Apps/Overleaf/Offset_Nature/Figure2_AllProjs_pscoreML.eps", 
       device = cairo_ps,
       width = 10, height = 8)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, mean.prob.control.logit, mean.prob.treat.logit) %>% 
  pivot_longer(cols = contains('mean.prob'),
               names_to = "group",
               names_prefix = "mean.prob",
               values_to = "mean.prob") %>% 
  mutate(group = case_when(
    group == ".treat.logit" ~ "Participating locations",
    group == ".control.logit" ~ "Non-participating locations")) %>% 
  ggpubr::gghistogram(x = "mean.prob", color = "group", fill = "group", add = "median",
                      bins = 20, alpha = .2) %>% 
  ggpar(xlab = "Average propensity score", legend.title = "", palette = c("#0E7175FF", "#FD7901FF"))

ggsave("C:/Users/samue/Dropbox/Apps/Overleaf/Offset_Nature/Figure2_AllProjs_pscorenonML.eps", 
       device = cairo_ps,
       width = 10, height = 8)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, mean.prob.overlap.ML, mean.prob.overlap.nonML) %>% 
  pivot_longer(cols = contains('mean.prob.overlap'),
               names_to = "group",
               names_prefix = "mean.prob.overlap.",
               values_to = "overlap") %>% 
  ggpubr::ggdensity(x = "overlap", color = "group", fill = "group", add = "mean")

dml.fitted.result.summarized %>% 
  filter(year.str == "yearfinalv0") %>% 
  select(proj.id, clust.type, mean.prob.overlap.ML) %>% 
  ggpubr::ggdensity(x = "mean.prob.overlap.ML", color = "clust.type", fill = "clust.type", add = "mean")


dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str == "yearfinalv0") %>% 
  select(proj.id, contains("mean.prob.control")) %>% 
  pivot_longer(cols = contains('mean.prob.control'),
               names_to = "type",
               names_prefix = "mean.prob.",
               values_to = "mean.prob") %>% 
  ggpubr::ggdensity(x = "mean.prob", color = "type", fill = "type", add = "mean")

ggpubr::gghistogram(subset(dml.fitted.result.summarized, year.str == "yearfinal"),
                    "pval", color = "clust.type", fill = "clust.type")

ggpubr::ggdensity(subset(dml.fitted.result.summarized, year.str == "yearfinal"),
                  "rsquared.outcome", color = "clust.type", fill = "clust.type", add = "mean")

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal"),
                  "mean.prob.control", color = "clust.type", fill = "clust.type", add = 'mean')

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal"),
                  "mean.prob.control", color = "clust.type", fill = "clust.type", add = 'mean')

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal"),
                  "pval.f", color = "clust.type", fill = "clust.type")

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal" & clust.type != "noclust"),
                  "pval.f", color = "clust.type", fill = "clust.type")

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus2`, y = `t_stat_yearfinal`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus3`, y = `t_stat_year.minus1`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & pval.f < 0.05) %>% 
  select(proj.id, eff, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = eff, names_prefix = "eff_") %>% 
  ggplot(aes(x = `eff_year.minus2`, y = `eff_yearfinal`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus3`, y = `t_stat_year.minus1`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()



plot(subset(dml.fitted.result.summarized, !grepl("-", year.str.plusminus)))


custom_log_trans <- function() {
  trans_new("custom_log",
            transform = function (x) ( sign(x)*log(abs(x)+1) ),
            inverse = function (y) ( sign(y)*( exp(abs(y))-1) ),
            domain = c(-Inf,Inf))
}

custom_log_breaks <- function(x) {
  
  log10range <- sign(x)*log10(abs(x))
  
  log10range.int <- c()
  
  for (a in 1:2) {
    log10range.int[a] <- ifelse(log10range[a] < 0, floor(log10range[a]), ceiling(log10range[a]))
  }
  
  ints <- seq(log10range.int[1], log10range.int[2], by = 1)
  
  ints.to.breaks <- c()
  
  for (k in 1:length(ints)) {
    ints.to.breaks[k] <- sign(ints[k])*10^(abs(ints[k]))
  }
  
  return(ints.to.breaks)
}


ggplot(data = subset(dml.res.df, signif & YEAR.INT > 0)) +
  geom_density(aes(x = EFF, fill = priv.type.public), alpha = .3) +
  geom_vline(xintercept = 0, linetype = 'dotted')+
  labs(x = "ATT (Mg biomass/ha)", 
       y = 'Density',
       fill = 'Ownership') +
  # scale_color_manual(values = lll_palette('California')) +
  scale_x_continuous(trans = "custom_log", 
                     breaks = c(-100, -60, -20, -10, 0, 10, 20, 50),
                     minor_breaks = c(-100, -80, -60, -40, -20, -10, 0, 10, 20, 40, 60))+
  scale_fill_nejm() +
  # coord_fixed() +
  theme_bw() + 
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = .5, size = 15, face = 'bold'),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:6), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff))


plot(subset(dml.fitted.result.summarized, !grepl("-", year.str.plusminus)))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_x_discrete(limits = c("-3", "-2", "+1", "+2", "+3", "+4")) +
  scale_y_continuous(trans = "custom_log",
                     breaks = c(-100, -50, -10, 0, 10, 50, 100)) +
  theme_bw()

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "noclust" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff)) +
  scale_x_discrete(limits = c("-3", "-2", "+1", "+2", "+3", "+4")) +
  scale_y_continuous(trans = "custom_log",
                     breaks = c(-100, -50, -10, 0, 10, 50, 100)) +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str %in% paste0("year.minus", c(2,3))) %>% 
  select(proj.id, pval.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = pval.f, names_prefix = "p_val_") %>% 
  ggplot(aes(x = `p_val_-2`, y = `p_val_-3`)) +
  geom_point() +
  geom_vline(xintercept = 0.05, linetype = 'dotted') +
  geom_hline(yintercept = 0.05, linetype = 'dotted') +
  lims(x = c(0, 0.25), y = c(0, 0.25)) +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str %in% paste0("year.minus", c(2,3))) %>% 
  select(proj.id, t.stat.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(y = `t_stat_-2`, x = `t_stat_-3`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_-2`, y = `t_stat_+2`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

plot(subset(dml.fitted.result.summarized, 
            year.str %in% paste0("year.minus", 2) & clust.type == "clust.25km" & pval.f < 0.05)$pval.f,
     subset(dml.fitted.result.summarized, 
            year.str %in% paste0("year.minus", 2) & clust.type == "clust.25km" & pval.f < 0.05)$pval.f)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & !grepl("year.minus", year.str)) %>% 
  group_by(proj.id) %>% 
  summarise(years.participated = n(),
            years.significant = sum(pval.f < 0.05)) -> years.particpiated.vs.sig

dml.fitted.result.summarized$year.final <- paste0("year", years.particpiated.vs.sig[match(dml.fitted.result.summarized$proj.id, years.particpiated.vs.sig$proj.id),]$years.participated)

dml.fitted.result.summarized$is.final.year <- dml.fitted.result.summarized$year.str == dml.fitted.result.summarized$year.final

table(subset(dml.fitted.result.summarized, is.final.year)$clust.type,
      subset(dml.fitted.result.summarized, is.final.year)$pval.f < 0.05)


length(subset(years.particpiated.vs.sig, years.significant >=1)$proj.id)

ggplot(years.particpiated.vs.sig) +
  geom_point(aes(x = years.participated, y = years.significant)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  lims(x = c(0, 13), y = c(0,13))

dml.signif.actual.v.placebo <- subset(dml.fitted.result.summarized, clust.type == "clust.25km" & pval.f < 0.05)
dml.signif.actual.v.placebo$placebo.or.not <- "Actual"

dml.signif.actual.v.placebo <- rbind(dml.signif.actual.v.placebo,
                                     cbind(subset(dml.fitted.result.summarized.placebo, clust.type == "clust.25km" & pval.f < 0.05),
                                           placebo.or.not = "Placebo"))

ggpubr::ggdensity(subset(dml.signif.actual.v.placebo),
                  "eff", color = "placebo.or.not", fill = "placebo.or.not")

# Calculating DRDID classic version ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]

cov.wo.clm.forminus.drdid <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         'fownership', "distance.to.road",
                         c(evi.cov.5), 
                         "biomass_tminus5")
#drop nlcd and forest.group because it causes singularities in OLS

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

drdid.fitted.result.summarized <- data.frame(proj.id = character(),
                                           year.str = character(),
                                           eff = numeric(),
                                           se = numeric())

projs.in.dat <- unique(proj.dats.for.prop.b2k$project.ID)

no.dat.projs <- c()

Mode.giver <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

B <- 500
#number of bootstraps

for (i in 1:length(projs.in.dat)) {
  
  projs.data.there <- paste0("Step3_DML_DIDClustwTminus3_", i, ".RData") %in% list.files()
  
  if (!projs.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinal")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    dat.for.drdid <- select(clust25km.dat$df,
                          c("cellID", "biomass_OUTCOME", "biomass_tminus5", "cluster.25km",
                            "treat", cov.wo.clm.forminus.drdid))
    
    clusters.w.treat <- unique(subset(dat.for.drdid, treat ==1)$cluster.25km)
    clusters.wo.treat <- setdiff(unique(dat.for.drdid$cluster.25km), clusters.w.treat)
    
    dat.for.drdid.split <- split(dat.for.drdid,
                                 dat.for.drdid$cluster.25km)
    
    
    dat.for.drdid %>% 
      pivot_longer(cols = starts_with("biomass_"),
                   values_to = "biomass",
                   names_to = "year",
                   names_prefix = "biomass_") %>% 
      relocate(cellID, biomass, treat) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(year.num = case_when(
        year == "OUTCOME" ~ 1,
        year == "tminus5" ~ 0
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(year.real =  case_when(
      year == "OUTCOME" ~ final.year.of.this.proj,
      year == "tminus5" ~ (proj.start.year - 5)
    )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(treat.num = case_when(
        treat == 1 ~ final.year.of.this.proj,
        treat == 0 ~ 0
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider$cellID <- as.numeric(dat.for.drdid.wider$cellID)
    
    dat.for.drdid.wider$biomass_tminus5 <- dat.for.drdid[match(dat.for.drdid.wider$cellID, dat.for.drdid$cellID),]$biomass_tminus5
      
    # drdid.result <- drdid(yname = "biomass",
    #        tname = "year.num",
    #        idname = "cellID",
    #        dname = "treat",
    #        xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
    #        data = dat.for.drdid.wider,
    #        estMethod = 'trad')
    #traditional for the use of OLS/logit
    
    ipwdid.result <- ipwdid(yname = "biomass",
                           tname = "year.num",
                           idname = "cellID",
                           dname = "treat",
                           xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus, collapse = "+"))),
                           data = dat.for.drdid.wider,
                           normalized = F)
    #normalized=F required for Abadie (2005) method
    
    # ipwdid.result <- att_gt(yname = "biomass",
    #                         tname = "year.real",
    #                         idname = "cellID",
    #                         gname = "treat.num",
    #                         control_group = "notyettreated",
    #                         clustervars = NULL,
    #                         bstrap = T,
    #                         xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
    #                         data = dat.for.drdid.wider)
    #doubly robust, with standard errors clustered
    
    drdid.fitted.result.summarized <- rbind(drdid.fitted.result.summarized,
                                          data.frame(proj.id = proj.here,
                                                     year.str = str,
                                                     eff = ipwdid.result$att,
                                                     se = ipwdid.result$se
                                          )
    )
  }
  
  projs.0.data.there <- paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData") %in% list.files()
  
  if (!projs.0.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust.v0) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust.v0[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinalv0", "year0v5")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    dat.for.drdid <- select(clust25km.dat$df,
                            c("cellID", "biomass_OUTCOME", "biomass_0", "treat", "cluster.25km",
                              cov.wo.clm.forminus.drdid))
    
    dat.for.drdid %>% 
      pivot_longer(cols = c("biomass_OUTCOME", "biomass_0"),
                   values_to = "biomass",
                   names_to = "year",
                   names_prefix = "biomass_") %>% 
      relocate(cellID, biomass, treat) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(year.num = case_when(
        year == "OUTCOME" ~ 1,
        year == "0" ~ 0
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(year.num = case_when(
        year == "OUTCOME" ~ 1,
        year == "0" ~ 0
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(year.real =  case_when(
        year == "OUTCOME" ~ final.year.of.this.proj,
        year == "0" ~ (proj.start.year)
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider %>% 
      mutate(treat.num = case_when(
        treat == 1 ~ final.year.of.this.proj,
        treat == 0 ~ 0
      )) -> dat.for.drdid.wider
    
    dat.for.drdid.wider$cellID <- as.numeric(dat.for.drdid.wider$cellID)
    
    # drdid.result <- drdid(yname = "biomass",
    #                       tname = "year.num", 
    #                       idname = "cellID", 
    #                       dname = "treat",
    #                       xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus, collapse = "+"))),
    #                       data = dat.for.drdid.wider, 
    #                       estMethod = 'trad')
    #OLS/logit
    
    # ipwdid.result <- ipwdid(yname = "biomass",
    #                         tname = "year.num", 
    #                         idname = "cellID", 
    #                         dname = "treat",
    #                         xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus, collapse = "+"))),
    #                         data = dat.for.drdid.wider, 
    #                         normalized = F)
    #normalized=F required for Abadie (2005) method
    
    ipwdid.result <- att_gt(yname = "biomass",
                            tname = "year.num",
                            idname = "cellID",
                            gname = "treat",
                            control_group = "nevertreated",
                            clustervars = "cluster.25km",
                            bstrap = T,
                            xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
                            data = dat.for.drdid.wider)
    
    drdid.fitted.result.summarized <- rbind(drdid.fitted.result.summarized,
                                            data.frame(proj.id = proj.here,
                                                       year.str = str,
                                                       eff = ipwdid.result$att,
                                                       se = ipwdid.result$se
                                            )
    )
    
  }
  
  print(paste0("Done with proj no. ", i))
  
  rm(dml.fitted.res.trim.forest.didclust)
  
}

save(list=c("drdid.fitted.result.summarized"),
     file = "Step3_DRDIDClassic_summarized.RData")


# Calculating DRDID classic version (clustered) ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]

cov.wo.clm.forminus <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         "forest.group",
                         "nlcd", 'fownership', "distance.to.road",
                         c(evi.cov.5), 
                         "biomass_tminus5"
)
#drop nlcd and forest.group because it causes singularities in OLS

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

projs.in.dat <- unique(proj.dats.for.prop.b2k$project.ID)

no.dat.projs <- c()

Mode.giver <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

B <- 500
#number of bootstraps

error.projs.dataset <- c()

drdid.fitted.result.summarized <- data.frame(proj.id = character(),
                                             year.str = character(),
                                             eff = numeric(),
                                             se = numeric())


vec.of.effs.list <- list()

for (i in c(30,79)) {
  #skipped through 30/79, because it takes way too much time
  
  tryCatch( {
  
  projs.data.there <- paste0("Step3_DML_DIDClustwTminus3_", i, ".RData") %in% list.files()
  
  if (!projs.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinal")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    dat.for.drdid <- select(clust25km.dat$df,
                            c("cellID", "biomass_OUTCOME", "biomass_tminus5", "cluster.25km",
                              "treat", cov.wo.clm.forminus))
    
    clusters.w.treat <- unique(subset(dat.for.drdid, treat ==1)$cluster.25km)
    clusters.wo.treat <- setdiff(unique(dat.for.drdid$cluster.25km), clusters.w.treat)
    
    dat.for.drdid.split <- split(dat.for.drdid,
                                 dat.for.drdid$cluster.25km)
    
    vec.of.effs <- c()
    
    tic()
    
    for (b in 1:B) {
      clusters.treat.chosen <- sample(clusters.w.treat, length(clusters.w.treat), replace =T)
      clusters.wotreat.chosen <- sample(clusters.wo.treat, length(clusters.wo.treat), replace =T)
      
      dat.for.drdid.this.time <- bind_rows(dat.for.drdid.split[unique(as.character(c(clusters.treat.chosen, clusters.wotreat.chosen)))])
      
      dat.for.drdid.this.time %>% 
        pivot_longer(cols = starts_with("biomass_"),
                     values_to = "biomass",
                     names_to = "year",
                     names_prefix = "biomass_") %>% 
        relocate(cellID, biomass, treat) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(year.num = case_when(
          year == "OUTCOME" ~ 1,
          year == "tminus5" ~ 0
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(year.real =  case_when(
          year == "OUTCOME" ~ final.year.of.this.proj,
          year == "tminus5" ~ (proj.start.year - 5)
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(treat.num = case_when(
          treat == 1 ~ final.year.of.this.proj,
          treat == 0 ~ 0
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider$cellID <- as.numeric(dat.for.drdid.wider$cellID)
      
      dat.for.drdid.wider$biomass_tminus5 <- dat.for.drdid[match(dat.for.drdid.wider$cellID, dat.for.drdid$cellID),]$biomass_tminus5
      
      # drdid.result <- drdid(yname = "biomass",
      #        tname = "year.num",
      #        idname = "cellID",
      #        dname = "treat",
      #        xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
      #        data = dat.for.drdid.wider,
      #        estMethod = 'trad')
      #traditional for the use of OLS/logit
      
      ipwdid.result <- DRDID::ipwdid(yname = "biomass",
                              tname = "year.num",
                              idname = "cellID",
                              dname = "treat",
                              xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus, collapse = "+"))),
                              data = dat.for.drdid.wider,
                              normalized = F)
      #normalized=F required for Abadie (2005) method
      
      # ipwdid.result <- att_gt(yname = "biomass",
      #                         tname = "year.real",
      #                         idname = "cellID",
      #                         gname = "treat.num",
      #                         control_group = "notyettreated",
      #                         clustervars = NULL,
      #                         bstrap = T,
      #                         xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
      #                         data = dat.for.drdid.wider)
      #doubly robust, with standard errors clustered
      
      vec.of.effs <- c(vec.of.effs,
                       ipwdid.result$ATT)
      
      if (b %% 100 == 0) {
        print(paste0("DONE WITH BOOTSTRAP NO. ", b))
        toc()
        tic()
        
        vec.of.effs.list[[paste0(proj.here, ".", str)]] <- vec.of.effs
      }
    }
    
    drdid.fitted.result.summarized <- rbind(drdid.fitted.result.summarized,
                                            data.frame(proj.id = proj.here,
                                                       year.str = str,
                                                       eff = mean(vec.of.effs, na.rm = T),
                                                       se = sd(vec.of.effs, na.rm = T)
                                            )
    )
    
    
    
  }
  
  projs.0.data.there <- paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData") %in% list.files()
  
  if (!projs.0.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_v0_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust.v0) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust.v0[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinalv0")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    # clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    dat.for.drdid <- select(clust25km.dat$df,
                            c("cellID", "biomass_OUTCOME", "biomass_0", "treat", "cluster.25km",
                              cov.wo.clm.forminus))
    
    dat.for.drdid.split <- split(dat.for.drdid,
                                 dat.for.drdid$cluster.25km)
    
    vec.of.effs <- c()
    
    for (b in 1:B) {
      
      clusters.treat.chosen <- sample(clusters.w.treat, length(clusters.w.treat), replace =T)
      clusters.wotreat.chosen <- sample(clusters.wo.treat, length(clusters.wo.treat), replace =T)
      
      dat.for.drdid.this.time <- bind_rows(dat.for.drdid.split[unique(as.character(c(clusters.treat.chosen, clusters.wotreat.chosen)))])
      
      dat.for.drdid.this.time %>% 
        pivot_longer(cols = c("biomass_OUTCOME", "biomass_0"),
                     values_to = "biomass",
                     names_to = "year",
                     names_prefix = "biomass_") %>% 
        relocate(cellID, biomass, treat) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(year.num = case_when(
          year == "OUTCOME" ~ 1,
          year == "0" ~ 0
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(year.real =  case_when(
          year == "OUTCOME" ~ final.year.of.this.proj,
          year == "0" ~ (proj.start.year - 5)
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider %>% 
        mutate(treat.num = case_when(
          treat == 1 ~ final.year.of.this.proj,
          treat == 0 ~ 0
        )) -> dat.for.drdid.wider
      
      dat.for.drdid.wider$cellID <- as.numeric(dat.for.drdid.wider$cellID)
      
      # drdid.result <- drdid(yname = "biomass",
      #        tname = "year.num",
      #        idname = "cellID",
      #        dname = "treat",
      #        xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
      #        data = dat.for.drdid.wider,
      #        estMethod = 'trad')
      #traditional for the use of OLS/logit
      
      ipwdid.result <- DRDID::ipwdid(yname = "biomass",
                                     tname = "year.num",
                                     idname = "cellID",
                                     dname = "treat",
                                     xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus, collapse = "+"))),
                                     data = dat.for.drdid.wider,
                                     normalized = F)
      #normalized=F required for Abadie (2005) method
      
      # ipwdid.result <- att_gt(yname = "biomass",
      #                         tname = "year.real",
      #                         idname = "cellID",
      #                         gname = "treat.num",
      #                         control_group = "notyettreated",
      #                         clustervars = NULL,
      #                         bstrap = T,
      #                         xformla = as.formula(paste0("~", paste(cov.wo.clm.forminus.drdid, collapse = "+"))),
      #                         data = dat.for.drdid.wider)
      #doubly robust, with standard errors clustered
      
      vec.of.effs <- c(vec.of.effs,
                       ipwdid.result$ATT)
      
      if (b %% 100 == 0) {
        print(paste0("DONE WITH BOOTSTRAP NO. ", b))
        toc()
        tic()
        
        vec.of.effs.list[[paste0(proj.here, ".", str)]] <- vec.of.effs
      }
    }
    
    drdid.fitted.result.summarized <- rbind(drdid.fitted.result.summarized,
                                            data.frame(proj.id = proj.here,
                                                       year.str = str,
                                                       eff = mean(vec.of.effs, na.rm = T),
                                                       se = sd(vec.of.effs, na.rm = T)
                                            )
    )
    
    rm(vec.of.effs)
    
  }
  
  }, error = function(e){
    error.projs.dataset <<- c(error.projs.dataset, i)
    
    print(paste0("ERROR IN PROJ NO. ", i, "!! SKIPPING..."))
  }) 
  
  print(paste0("Done with proj no. ", i))
  
  rm(dml.fitted.res.trim.forest.didclust)
  
  save(list=c("drdid.fitted.result.summarized", "vec.of.effs.list", "error.projs.dataset"),
       file = "Step3_DRDIDClassicClustered_summarized.RData")
  
}

drdid.fitted.result.summarized$treated.N <- dml.fitted.result.summarized[match(drdid.fitted.result.summarized$proj.id, dml.fitted.result.summarized$proj.id),]$treated.N

save(list=c("drdid.fitted.result.summarized"),
     file = "Step3_DRDIDClassic_summarized.RData")


# Fit DML for deltas, DiD ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")

load("Step2_1_Fixed characteristics_NDVI250m.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step2.5_NonForestAndTreatedPixels.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

projs.in.dat <- unique(proj.dats.for.prop.b2k$proj)

projs.with.1.possible <- proj.years[which(year(proj.years$DATE.first) <= 2016),]$project.name
projs.with.2.possible <- proj.years[which(year(proj.years$DATE.first) <= 2015),]$project.name
projs.with.3.possible <- proj.years[which(year(proj.years$DATE.first) <= 2014),]$project.name
projs.with.4.possible <- proj.years[which(year(proj.years$DATE.first) <= 2013),]$project.name

evi.cov.2 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_2", colnames(proj.dats.for.prop.b2k)))]
evi.cov.3 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_3", colnames(proj.dats.for.prop.b2k)))]
evi.cov.4 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_4", colnames(proj.dats.for.prop.b2k)))]
evi.cov.5 <- colnames(proj.dats.for.prop.b2k)[which(grepl("evi.*ns.coef|ns.coef*evi", colnames(proj.dats.for.prop.b2k)) & grepl("\\_5", colnames(proj.dats.for.prop.b2k)))]

cov.wo.clm.forminus.delta <- c("clm_DEM", "clm_ned_lf", "clm_soilCLAY", "clm_soilpH", "clm_soilORGC",
                         "forest.group",
                         "nlcd", 'fownership', "distance.to.road",
                         c(evi.cov.5), "biomass_tminus5.delta"
)

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

# load("Step3_DMLDiDClustered_summarized.RData")
# 
# proj.with.positive.final.att <- subset(dml.fitted.result.summarized, clust.type == "clust.25km" & pval.f < 0.05 & year.str == "yearfinal" & eff > 0)$proj.id

for (i in c(1:length(projs.in.dat))) {
  
  tryCatch( {tic()
    
    proj.to.inspect <- projs.in.dat[i]
    
    if (proj.to.inspect %in% red.df$redundant.proj) {
      print(paste0("SKIPPING THIS PROJECT BECAUSE IT'S REDUNADNT!!"))
      #if this is a redundant project, then its "true" start year is not the PLACEBO.YEAR
      #it's sometime before this --> so you must drop this
      next
    }
    
    print(paste0("STARTING WITH PROJ.", proj.to.inspect, "!!!"))
    
    dat.to.be.used.list <- dml.dat.forest.didclust[[proj.to.inspect]]
    
    minus.years <- names(dat.to.be.used.list)[grepl("year.minus|year0", names(dat.to.be.used.list))]
    
    dat.to.be.used.minus <- bind_rows(dat.to.be.used.list[[minus.years[1]]])
    
    dat.to.be.used.minus$Year <- minus.years[1]
    
    for (my in minus.years[2:length(minus.years)]) {
      dat.to.be.used.minus.thismy <- dat.to.be.used.list[[my]]
      dat.to.be.used.minus.thismy$Year <- my
      
      dat.to.be.used.minus <- bind_rows(dat.to.be.used.minus,
                                        dat.to.be.used.minus.thismy)
    }
    
    dat.to.be.used.minus$cluster.1km.year <- paste0(dat.to.be.used.minus$cluster.1km, "-", dat.to.be.used.minus$Year)
    dat.to.be.used.minus$cluster.25km.year <- paste0(dat.to.be.used.minus$cluster.25km, "-", dat.to.be.used.minus$Year)
    dat.to.be.used.minus$cluster.5km.year <- paste0(dat.to.be.used.minus$cluster.5km, "-", dat.to.be.used.minus$Year)
    
    plus.years <- setdiff(names(dat.to.be.used.list), minus.years)
    
    dat.to.be.used.plus <- bind_rows(dat.to.be.used.list[[plus.years[1]]])
    
    dat.to.be.used.plus$Year <- plus.years[1]
    
    for (py in plus.years[2:length(plus.years)]) {
      dat.to.be.used.plus.thispy <- dat.to.be.used.list[[py]]
      dat.to.be.used.plus.thispy$Year <- py
      
      dat.to.be.used.plus <- bind_rows(dat.to.be.used.plus,
                                       dat.to.be.used.plus.thispy)
    }
    
    dat.to.be.used.plus$cluster.1km.year <- paste0(dat.to.be.used.plus$cluster.1km, "-", dat.to.be.used.plus$Year)
    dat.to.be.used.plus$cluster.25km.year <- paste0(dat.to.be.used.plus$cluster.25km, "-", dat.to.be.used.plus$Year)
    dat.to.be.used.plus$cluster.5km.year <- paste0(dat.to.be.used.plus$cluster.5km, "-", dat.to.be.used.plus$Year)
    
    dml.fitted.res.trim.forest.didclust.delta <- list()
    
    dml.fitted.res.trim.forest.didclust.delta[[proj.to.inspect]] <- list()
    
    proj.fitted.models.trim.forest.didclust <- list()
    
    # dat.list.to.use.minus.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used.minus,
    #                                                                     cluster.use = T, cluster.cols = 'cluster.1km.year', dupl.drop = F)
    
    dat.list.to.use.minus.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used.minus,
                                                                       cluster.use = T, cluster.cols = 'cluster.25km.year', dupl.drop = F)
    
    # dat.list.to.use.minus.clust5km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used.minus,
    #                                                                   cluster.use = T, cluster.cols = 'cluster.5km.year', dupl.drop = F)
    
    
    # dat.list.to.use.plus.clust1km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used.plus,
    #                                                                    cluster.use = T, cluster.cols = 'cluster.1km.year', dupl.drop = F)
    dat.list.to.use.plus.clust25km <- dml.datamaker.for.didclust.list( dat.to.use = dat.to.be.used.plus,
                                                                  cluster.use = T, cluster.cols = 'cluster.25km.year', dupl.drop = F)
    
    
    # proj.fitted.models.trim.forest.didclust[["yearminus.clust1km"]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
    #                                                                                                    dat.list.to.use = dat.list.to.use.minus.clust1km, 
    #                                                                                                    cov.to.use = cov.wo.clm.forminus.delta, 
    #                                                                                                    outcome.to.use = "biomass_OUTCOME.delta", 
    #                                                                                                    baseline.string = "biomass_tminus5.delta",
    #                                                                                                    diff.with.baseline = F, 
    #                                                                                                    g.learner = "classif.rf",
    #                                                                                                    m.learner = "regr.rf",
    #                                                                                                    trim.use = T
    # )
    
    proj.fitted.models.trim.forest.didclust[["yearminus.clust25km"]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                            dat.list.to.use = dat.list.to.use.minus.clust25km, 
                                                                                                            cov.to.use = cov.wo.clm.forminus.delta, 
                                                                                                            outcome.to.use = "biomass_OUTCOME.delta", 
                                                                                                            baseline.string = "biomass_tminus5.delta",
                                                                                                            diff.with.baseline = F, 
                                                                                                            g.learner = "classif.rf",
                                                                                                            m.learner = "regr.rf",
                                                                                                            trim.use = T
    )
    
    # proj.fitted.models.trim.forest.didclust[["yearplus.clust1km"]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
    #                                                                                                        dat.list.to.use = dat.list.to.use.plus.clust1km, 
    #                                                                                                        cov.to.use = cov.wo.clm.forminus.delta, 
    #                                                                                                        outcome.to.use = "biomass_OUTCOME.delta", 
    #                                                                                                        baseline.string = "biomass_tminus5.delta",
    #                                                                                                        diff.with.baseline = F, 
    #                                                                                                        g.learner = "classif.rf",
    #                                                                                                        m.learner = "regr.rf",
    #                                                                                                        trim.use = T
    # )
    
    proj.fitted.models.trim.forest.didclust[["yearplus.clust25km"]] <- dml.runner.for.didclust.baselinestr(proj.to.look.at = proj.to.inspect, 
                                                                                                            dat.list.to.use = dat.list.to.use.plus.clust25km, 
                                                                                                            cov.to.use = cov.wo.clm.forminus.delta, 
                                                                                                            outcome.to.use = "biomass_OUTCOME.delta", 
                                                                                                            baseline.string = "biomass_tminus5.delta",
                                                                                                            diff.with.baseline = F, 
                                                                                                            g.learner = "classif.rf",
                                                                                                            m.learner = "regr.rf",
                                                                                                            trim.use = T
    )
    
    
    dml.fitted.res.trim.forest.didclust.delta[[proj.to.inspect]] <- proj.fitted.models.trim.forest.didclust
    
    
    save(list = c("dml.fitted.res.trim.forest.didclust.delta"),
         file = paste0("Step3_DML_DIDClust_wDelta_", i, ".RData"))
    
    toc()
    tic()
    
    print(paste0("DONE WITH PROJECT NO. ", i))
    
    #restore empty lists
    dml.fitted.res.trim.forest.didclust.delta <- list()
    
  }, error = function(e){
    
    if (length(proj.fitted.models.trim.forest.didclust)!=0) {
      dml.fitted.res.trim.forest.didclust.delta[[proj.to.inspect]] <<- proj.fitted.models.trim.forest.didclust
      
      save(list = c("dml.fitted.res.trim.forest.didclust.delta"),
           file = paste0("Step3_DML_DIDClust_wDelta_", i, ".RData"))
      
      dml.fitted.res.trim.forest.didclust.delta <- list()
    }
    
    
    print(paste0("ERROR IN PROJ NO. ", i, "!! SKIPPING..."))
  }   )
  
}

# Evaluate results ======

load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

dml.fitted.delta.summarized <- data.frame(proj.id = character(),
                                                   # year.calendar = integer(),
                                                   year.str = character(),
                                                   # year.str.plusminus = character(),
                                                   clust.type = character(),
                                                   eff = numeric(),
                                                   se = numeric(),
                                                   trimmed.N = integer(),
                                                   cluster.count.total = integer(),
                                                   cluster.count = integer(),
                                                   pval = numeric(),
                                                   rsquared.outcome = numeric(),
                                                   mean.prob.control = numeric())


for (i in 1:111) {
  
  projs.data.there <- paste0("Step3_DML_DIDClust_wDelta_", i, ".RData") %in% list.files()
  
  if (!projs.data.there) {
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClust_wDelta_", i, ".RData"))
  #load the data
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust.delta)
  
  proj.dat <- dml.fitted.res.trim.forest.didclust.delta[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  for (str in c("yearplus", "yearminus")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    # clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    # noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    dml.fitted.delta.summarized <- rbind(dml.fitted.delta.summarized,
                                                  data.frame(proj.id = rep(proj.here, 2),
                                                             # year.calendar = rep(proj.start.year, 2) + year.integer,
                                                             year.str = rep(str,2),
                                                             # year.str.plusminus = rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                                             #                                 as.character(year.integer)), 4),
                                                             clust.type = c("clust.1km",
                                                                            "clust.25km"),
                                                             eff = c(ifelse(is.null(clust1km.dat), NA, clust1km.dat$eff),
                                                                     ifelse(is.null(clust25km.dat), NA, clust25km.dat$eff)),
                                                             se = c(ifelse(is.null(clust1km.dat), NA, clust1km.dat$se),
                                                                    ifelse(is.null(clust25km.dat), NA, clust25km.dat$se)),
                                                             trimmed.N = c(ifelse(is.null(clust1km.dat), NA, sum(as.integer(clust1km.dat$df$ghat < 0.99 & clust1km.dat$df$ghat > 0.01), na.rm = T)),
                                                                           ifelse(is.null(clust25km.dat), NA, sum(as.integer(clust25km.dat$df$ghat < 0.99 & clust25km.dat$df$ghat > 0.01), na.rm = T))),
                                                             total.N = c(ifelse(is.null(clust1km.dat), NA, nrow(clust1km.dat$df)),
                                                                         ifelse(is.null(clust25km.dat), NA, nrow(clust25km.dat$df))),
                                                             # cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                             #                   ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km))),
                                                             #                   ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km))),
                                                             #                   ifelse(is.null(clust5km.dat), NA, length(unique(clust5km.dat$df$cluster.5km)))),
                                                             cluster.count.total = c(ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km.year))),
                                                                                     ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km.year)))),
                                                             cluster.count = c(ifelse(is.null(clust1km.dat), NA, length(unique(subset(clust1km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.1km.year))),
                                                                               ifelse(is.null(clust25km.dat), NA, length(unique(subset(clust25km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.25km.year)))),
                                                             pval = c(ifelse(is.null(clust1km.dat), NA, clust1km.dat$p.val),
                                                                      ifelse(is.null(clust25km.dat), NA, clust25km.dat$p.val)),
                                                             rsquared.outcome = c(ifelse(is.null(clust1km.dat), NA, cor(clust1km.dat$df$biomass_OUTCOME.delta - clust1km.dat$df$biomass_tminus5.delta,
                                                                                                                        clust1km.dat$df$ellhat)),
                                                                                  ifelse(is.null(clust25km.dat), NA, cor(clust25km.dat$df$biomass_OUTCOME.delta - clust25km.dat$df$biomass_tminus5.delta,
                                                                                                                         clust25km.dat$df$ellhat))),
                                                             mean.prob.control = c(ifelse(is.null(clust1km.dat), NA, mean(subset(clust1km.dat$df, treat==0)$ghat, na.rm =T)),
                                                                                   ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T))
                                                             ))
    
    )
  }
  
  print(paste0("Done with proj no. ", i))
  
  rm(dml.fitted.res.trim.forest.didclust.delta)
  
}

dml.fitted.delta.summarized$se.f <- dml.fitted.delta.summarized$se*sqrt(dml.fitted.delta.summarized$total.N)/sqrt(dml.fitted.delta.summarized$cluster.count.total)

dml.fitted.delta.summarized$t.stat.f <- dml.fitted.delta.summarized$eff/dml.fitted.delta.summarized$se.f

pval.calc.function <- function(x,y) {2*pt(x, df = y, lower.tail = (x<=0))}

dml.fitted.delta.summarized$pval.f <- mapply(pval.calc.function,
                                                      dml.fitted.delta.summarized$t.stat.f,
                                                      dml.fitted.delta.summarized$trimmed.N)

dml.fitted.delta.summarized <- subset(dml.fitted.delta.summarized, 
                                      clust.type == "clust.25km")

table(dml.fitted.delta.summarized$pval.f < 0.05,
      dml.fitted.delta.summarized$year.str,
      dml.fitted.delta.summarized$eff > 0)

table(dml.fitted.delta.summarized$clust.type,
      dml.fitted.delta.summarized$pval.f < 0.05,
      dml.fitted.delta.summarized$year.str)

View(subset(dml.fitted.delta.summarized, clust.type == "clust.25km" & pval.f < 0.05))

ggpubr::ggdensity(subset(dml.fitted.delta.summarized, year.str == "yearplus"),
                    "mean.prob.control", color = "clust.type", fill = "clust.type")

save(list = "dml.fitted.delta.summarized",
     file = "Step3_DMLDiDClustered_Delta_summarized.RData")




# Evaluate results (deprecated 24/09/18) ======

load("Step3_DML_DIDClust_DataBiomassTminus3_111.RData")
load("Step2.5_ProjectSpecificData_b2k.RData")
load("Step3_REDUNDANT_PROJECTS.RData")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

factor.covs <- c("fownership", "nlcd", "clm_ned_lf", "forest.group")

dml.fitted.result.summarized <- data.frame(proj.id = character(),
                                           fownership.type = character(),
                                           forest.type = character(),
                                           # year.calendar = rep(proj.start.year, 4) + year.integer,
                                           year.start = integer(),
                                           year.final = integer(),
                                           year.str = character(),
                                           # year.str.plusminus = rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                           #                                 as.character(year.integer)), 4),
                                           clust.type = character(),
                                           eff = numeric(),
                                           se = numeric(),
                                           trimmed.N = integer(),
                                           cluster.count.total = integer(),
                                           cluster.count = integer(),
                                           pval = numeric(),
                                           rsquared.outcome = numeric(),
                                           mean.prob.control = numeric())

projs.in.dat <- unique(proj.dats.for.prop.b2k$project.ID)

no.dat.projs <- c()

Mode.giver <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

for (i in 1:length(projs.in.dat)) {
  
  projs.data.there <- paste0("Step3_DML_DIDClustwTminus3_", i, ".RData") %in% list.files()
  
  if (!projs.data.there) {
    no.dat.projs <- c(no.dat.projs, projs.in.dat[i]) 
    next
    #if this project was skipped, move onto the next
  }
  
  load(paste0("Step3_DML_DIDClustwTminus3_", i, ".RData"))
  #load the data
  
  # proj.here <- names(dml.fitted.res.trim.forest.didclust)[i]
  
  proj.here <- names(dml.fitted.res.trim.forest.didclust) 
  #a single project within each file
  
  proj.dat <- dml.fitted.res.trim.forest.didclust[[proj.here]]
  
  proj.start.year <- year(proj.years[proj.years$project.name == proj.here, ]$DATE.first)
  
  final.year.of.this.proj <- if (proj.here %in% red.df$original.proj) {
    later.project <- red.df[red.df$original.proj==proj.here,]$redundant.proj
    
    min(year(proj.years[proj.years$project.name == later.project,]$DATE.last),
        2017)
    
  } else {
    min(year(proj.years[proj.years$project.name == proj.here,]$DATE.last),
        2017)
  }
  
  for (str in c("yearfinal", paste0("year.minus", 3:1), "year0")) {
    if (!any(grepl(str, names(proj.dat)))) {
      next
    }
    
    # year.integer <- as.integer(ifelse(grepl("minus", str), 0-as.integer(gsub("year.minus", "", str)),
    #                                   gsub("year", "", str)))
    
    clust1km.dat <- proj.dat[[paste0(str, ".clust1km")]]
    clust25km.dat <- proj.dat[[paste0(str, ".clust25km")]]
    clust5km.dat <- proj.dat[[paste0(str, ".clust5km")]]
    noclust.dat <- proj.dat[[paste0(str, ".noclust")]]
    
    avg.biomass <- mean(subset(noclust.dat$df, treat ==1)$biomass_0, na.rm = T)
    
    dml.fitted.result.summarized <- rbind(dml.fitted.result.summarized,
                                          data.frame(proj.id = rep(proj.here, 4),
                                                     fownership.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$fownership),4),
                                                     forest.type = rep(Mode.giver(subset(noclust.dat$df, treat ==1)$forest.group),4),
                                                     # year.calendar = rep(proj.start.year, 4) + year.integer,
                                                     year.start = rep(proj.start.year, 4),
                                                     year.final = rep(final.year.of.this.proj, 4),
                                                     year.str = rep(str,4),
                                                     # year.str.plusminus = rep(ifelse(year.integer > 0, paste0("+",year.integer),
                                                     #                                 as.character(year.integer)), 4),
                                                     clust.type = c("noclust", "clust.1km",
                                                                    "clust.25km", "clust.5km"),
                                                     avg.biomass = rep(avg.biomass, 4),
                                                     eff = c(ifelse(is.null(noclust.dat), NA, noclust.dat$eff),
                                                             ifelse(is.null(clust1km.dat), NA, clust1km.dat$eff),
                                                             ifelse(is.null(clust25km.dat), NA, clust25km.dat$eff),
                                                             ifelse(is.null(clust5km.dat), NA, clust5km.dat$eff)),
                                                     se = c(ifelse(is.null(noclust.dat), NA, noclust.dat$se),
                                                            ifelse(is.null(clust1km.dat), NA, clust1km.dat$se),
                                                            ifelse(is.null(clust25km.dat), NA, clust25km.dat$se),
                                                            ifelse(is.null(clust5km.dat), NA, clust5km.dat$se)),
                                                     trimmed.N = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                   ifelse(is.null(clust1km.dat), NA, sum(as.integer(clust1km.dat$df$ghat < 0.99 & clust1km.dat$df$ghat > 0.01), na.rm = T)),
                                                                   ifelse(is.null(clust25km.dat), NA, sum(as.integer(clust25km.dat$df$ghat < 0.99 & clust25km.dat$df$ghat > 0.01), na.rm = T)),
                                                                   ifelse(is.null(clust5km.dat), NA, sum(as.integer(clust5km.dat$df$ghat < 0.99 & clust5km.dat$df$ghat > 0.01), na.rm = T))),
                                                     total.N = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                 ifelse(is.null(clust1km.dat), NA, nrow(clust1km.dat$df)),
                                                                 ifelse(is.null(clust25km.dat), NA, nrow(clust25km.dat$df)),
                                                                 ifelse(is.null(clust5km.dat), NA, nrow(clust5km.dat$df))),
                                                     # cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                     #                   ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km))),
                                                     #                   ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km))),
                                                     #                   ifelse(is.null(clust5km.dat), NA, length(unique(clust5km.dat$df$cluster.5km)))),
                                                     cluster.count.total = c(ifelse(is.null(noclust.dat), NA, nrow(noclust.dat$df)),
                                                                             ifelse(is.null(clust1km.dat), NA, length(unique(clust1km.dat$df$cluster.1km))),
                                                                             ifelse(is.null(clust25km.dat), NA, length(unique(clust25km.dat$df$cluster.25km))),
                                                                             ifelse(is.null(clust5km.dat), NA, length(unique(clust5km.dat$df$cluster.5km)))),
                                                     cluster.count = c(ifelse(is.null(noclust.dat), NA, sum(as.integer(noclust.dat$df$ghat < 0.99 & noclust.dat$df$ghat > 0.01), na.rm = T)),
                                                                       ifelse(is.null(clust1km.dat), NA, length(unique(subset(clust1km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.1km))),
                                                                       ifelse(is.null(clust25km.dat), NA, length(unique(subset(clust25km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.25km))),
                                                                       ifelse(is.null(clust5km.dat), NA, length(unique(subset(clust5km.dat$df, ghat < 0.99 & ghat > 0.01)$cluster.5km)))),
                                                     pval = c(ifelse(is.null(noclust.dat), NA, noclust.dat$p.val),
                                                              ifelse(is.null(clust1km.dat), NA, clust1km.dat$p.val),
                                                              ifelse(is.null(clust25km.dat), NA, clust25km.dat$p.val),
                                                              ifelse(is.null(clust5km.dat), NA, clust5km.dat$p.val)),
                                                     rsquared.outcome = c(ifelse(is.null(noclust.dat), NA, cor(noclust.dat$df$biomass_OUTCOME - noclust.dat$df$biomass_tminus5,
                                                                                                               noclust.dat$df$ellhat)),
                                                                          ifelse(is.null(clust1km.dat), NA, cor(clust1km.dat$df$biomass_OUTCOME - clust1km.dat$df$biomass_tminus5,
                                                                                                                clust1km.dat$df$ellhat)),
                                                                          ifelse(is.null(clust25km.dat), NA, cor(clust25km.dat$df$biomass_OUTCOME - clust25km.dat$df$biomass_tminus5,
                                                                                                                 clust25km.dat$df$ellhat)),
                                                                          ifelse(is.null(clust5km.dat), NA, cor(clust5km.dat$df$biomass_OUTCOME - clust5km.dat$df$biomass_tminus5,
                                                                                                                clust5km.dat$df$ellhat))
                                                     ),
                                                     mean.prob.control = c(ifelse(is.null(noclust.dat), NA, mean(subset(noclust.dat$df, treat==0)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust1km.dat), NA, mean(subset(clust1km.dat$df, treat==0)$ghat, na.rm =T)),
                                                                           ifelse(is.null(clust25km.dat), NA, mean(subset(clust25km.dat$df,treat==0)$ghat, na.rm = T)),
                                                                           ifelse(is.null(clust5km.dat), NA, mean(subset(clust5km.dat$df, treat==0)$ghat, na.rm = T))
                                                     ))
    )
    
    
    
  }
  
  print(paste0("Done with proj no. ", i))
  
  rm(dml.fitted.res.trim.forest.didclust)
  
}


dml.fitted.result.summarized$se.f <- dml.fitted.result.summarized$se*sqrt(dml.fitted.result.summarized$total.N)/sqrt(dml.fitted.result.summarized$cluster.count.total)

dml.fitted.result.summarized$t.stat.f <- dml.fitted.result.summarized$eff/dml.fitted.result.summarized$se.f

pval.calc.function <- function(x,y) {2*pt(x, df = y, lower.tail = (x<=0))}

dml.fitted.result.summarized$pval.f <- mapply(pval.calc.function,
                                              dml.fitted.result.summarized$t.stat.f,
                                              dml.fitted.result.summarized$cluster.count)

table(dml.fitted.result.summarized$clust.type,
      dml.fitted.result.summarized$pval.f < 0.05,
      dml.fitted.result.summarized$year.str)

View(subset(dml.fitted.result.summarized, clust.type == "clust.25km" & pval.f < 0.05 & year.str == "yearfinal"))
#note that the ones that are insignificant have very small effect estimates

save(list = "dml.fitted.result.summarized",
     file = "Step3_DMLDiDClustered_summarized.RData")

ggpubr::gghistogram(subset(dml.fitted.result.summarized, year.str == "yearfinal"),
                    "pval", color = "clust.type", fill = "clust.type")

ggpubr::ggdensity(subset(dml.fitted.result.summarized, year.str == "yearfinal"),
                  "rsquared.outcome", color = "clust.type", fill = "clust.type", add = "mean")

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal"),
                  "mean.prob.control", color = "clust.type", fill = "clust.type", add = 'mean')

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal"),
                  "pval.f", color = "clust.type", fill = "clust.type")

ggpubr::ggdensity(subset(dml.fitted.result.summarized,  year.str == "yearfinal" & clust.type != "noclust"),
                  "pval.f", color = "clust.type", fill = "clust.type")

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus2`, y = `t_stat_yearfinal`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus3`, y = `t_stat_year.minus1`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & pval.f < 0.05) %>% 
  select(proj.id, eff, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = eff, names_prefix = "eff_") %>% 
  ggplot(aes(x = `eff_year.minus2`, y = `eff_yearfinal`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_year.minus3`, y = `t_stat_year.minus1`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()



plot(subset(dml.fitted.result.summarized, !grepl("-", year.str.plusminus)))


custom_log_trans <- function() {
  trans_new("custom_log",
            transform = function (x) ( sign(x)*log(abs(x)+1) ),
            inverse = function (y) ( sign(y)*( exp(abs(y))-1) ),
            domain = c(-Inf,Inf))
}

custom_log_breaks <- function(x) {
  
  log10range <- sign(x)*log10(abs(x))
  
  log10range.int <- c()
  
  for (a in 1:2) {
    log10range.int[a] <- ifelse(log10range[a] < 0, floor(log10range[a]), ceiling(log10range[a]))
  }
  
  ints <- seq(log10range.int[1], log10range.int[2], by = 1)
  
  ints.to.breaks <- c()
  
  for (k in 1:length(ints)) {
    ints.to.breaks[k] <- sign(ints[k])*10^(abs(ints[k]))
  }
  
  return(ints.to.breaks)
}


ggplot(data = subset(dml.res.df, signif & YEAR.INT > 0)) +
  geom_density(aes(x = EFF, fill = priv.type.public), alpha = .3) +
  geom_vline(xintercept = 0, linetype = 'dotted')+
  labs(x = "ATT (Mg biomass/ha)", 
       y = 'Density',
       fill = 'Ownership') +
  # scale_color_manual(values = lll_palette('California')) +
  scale_x_continuous(trans = "custom_log", 
                     breaks = c(-100, -60, -20, -10, 0, 10, 20, 50),
                     minor_breaks = c(-100, -80, -60, -40, -20, -10, 0, 10, 20, 40, 60))+
  scale_fill_nejm() +
  # coord_fixed() +
  theme_bw() + 
  theme(legend.position = 'bottom',
        plot.title = element_text(hjust = .5, size = 15, face = 'bold'),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:6), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff))


plot(subset(dml.fitted.result.summarized, !grepl("-", year.str.plusminus)))

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "clust.25km" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff)) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  scale_x_discrete(limits = c("-3", "-2", "+1", "+2", "+3", "+4")) +
  scale_y_continuous(trans = "custom_log",
                     breaks = c(-100, -50, -10, 0, 10, 50, 100)) +
  theme_bw()

ggplot(subset(dml.fitted.result.summarized, 
              year.str %in% c(paste0("year", 1:4), paste0("year.minus", 2:3)) & clust.type == "noclust" & pval.f < 0.05)) +
  geom_boxplot(aes(x = year.str.plusminus, y = eff)) +
  scale_x_discrete(limits = c("-3", "-2", "+1", "+2", "+3", "+4")) +
  scale_y_continuous(trans = "custom_log",
                     breaks = c(-100, -50, -10, 0, 10, 50, 100)) +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str %in% paste0("year.minus", c(2,3))) %>% 
  select(proj.id, pval.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = pval.f, names_prefix = "p_val_") %>% 
  ggplot(aes(x = `p_val_-2`, y = `p_val_-3`)) +
  geom_point() +
  geom_vline(xintercept = 0.05, linetype = 'dotted') +
  geom_hline(yintercept = 0.05, linetype = 'dotted') +
  lims(x = c(0, 0.25), y = c(0, 0.25)) +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & year.str %in% paste0("year.minus", c(2,3))) %>% 
  select(proj.id, t.stat.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(y = `t_stat_-2`, x = `t_stat_-3`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km") %>% 
  select(proj.id, t.stat.f, year.str.plusminus) %>% 
  pivot_wider(id_cols = proj.id, names_from = year.str.plusminus, values_from = t.stat.f, names_prefix = "t_stat_") %>% 
  ggplot(aes(x = `t_stat_-2`, y = `t_stat_+2`)) +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 'dotted') +
  geom_point() +
  coord_fixed( ratio=1) +
  geom_abline(intercept =0, slope = 1, linetype = 'dashed') +
  theme_bw()

plot(subset(dml.fitted.result.summarized, 
            year.str %in% paste0("year.minus", 2) & clust.type == "clust.25km" & pval.f < 0.05)$pval.f,
     subset(dml.fitted.result.summarized, 
            year.str %in% paste0("year.minus", 2) & clust.type == "clust.25km" & pval.f < 0.05)$pval.f)

dml.fitted.result.summarized %>% 
  filter(clust.type == "clust.25km" & !grepl("year.minus", year.str)) %>% 
  group_by(proj.id) %>% 
  summarise(years.participated = n(),
            years.significant = sum(pval.f < 0.05)) -> years.particpiated.vs.sig

dml.fitted.result.summarized$year.final <- paste0("year", years.particpiated.vs.sig[match(dml.fitted.result.summarized$proj.id, years.particpiated.vs.sig$proj.id),]$years.participated)

dml.fitted.result.summarized$is.final.year <- dml.fitted.result.summarized$year.str == dml.fitted.result.summarized$year.final

table(subset(dml.fitted.result.summarized, is.final.year)$clust.type,
      subset(dml.fitted.result.summarized, is.final.year)$pval.f < 0.05)


length(subset(years.particpiated.vs.sig, years.significant >=1)$proj.id)

ggplot(years.particpiated.vs.sig) +
  geom_point(aes(x = years.participated, y = years.significant)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  lims(x = c(0, 13), y = c(0,13))

dml.signif.actual.v.placebo <- subset(dml.fitted.result.summarized, clust.type == "clust.25km" & pval.f < 0.05)
dml.signif.actual.v.placebo$placebo.or.not <- "Actual"

dml.signif.actual.v.placebo <- rbind(dml.signif.actual.v.placebo,
                                     cbind(subset(dml.fitted.result.summarized.placebo, clust.type == "clust.25km" & pval.f < 0.05),
                                           placebo.or.not = "Placebo"))

ggpubr::ggdensity(subset(dml.signif.actual.v.placebo),
                  "eff", color = "placebo.or.not", fill = "placebo.or.not")


### Update on confidence intervals (use percentile CI, rather than SD) - dperecated 24/10/05 =======

projs.there <- unique(gsub(".yearfinal|.yearfinalv0", "", names(vec.of.effs.list)))

projs.there <- projs.there[!grepl(".year0v5", projs.there)]

drdid.fitted.result.summarized <- data.frame(proj.id = character(),
                                             year.str = character(),
                                             eff = numeric(),
                                             eff.lb = numeric(),
                                             eff.ub = numeric())

for (yearstr in c("yearfinal", "yearfinalv0")) {
  
  for (p in projs.there) {
    
    drdid.dat <- vec.of.effs.list[[paste0(p, ".", yearstr)]]
    
    drdid.fitted.result.summarized <- bind_rows(drdid.fitted.result.summarized,
                                                data.frame(proj.id = p,
                                                           year.str = yearstr,
                                                           eff = mean(drdid.dat, na.rm = T),
                                                           eff.lb = quantile(drdid.dat, 0.025, na.rm = T),
                                                           eff.ub = quantile(drdid.dat, 0.975, na.rm = T)))
  }
  
  
}
