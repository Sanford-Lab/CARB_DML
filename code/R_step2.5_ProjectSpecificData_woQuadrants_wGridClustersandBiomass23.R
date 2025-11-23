library(tidyverse);library(lubridate);library(Rforestry);library(reshape2)
library(randomForest);library(tictoc);library(sf);library(paletteer)

#With biomass_2, biomass_3 instead of biomass_1, biomass_0
#and without quadrants

# load("Step2_1_Fixed characteristics_ALL.RData")
load("Step2_1_Fixed characteristics_NDVI250m.RData")
crosswalk <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv")
load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")
load("Step2_BIOMASS_NDVI250m.RData")

for (yyyy in 2000:2018) {
  load(paste0("Step2_2_NDVI_SUMMARIZED_", yyyy, ".RData"))
}

projs.to.look.at <- proj.years[year(proj.years$DATE.first)<=2016,]$project.name

nlcd.years <- c(2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019)

colnames(biomass.extracted.all) <- c(paste0("biomass_", 2000:2017), "cellID")

data.maker <- function(proj.name, treat.type, control.type, flag.qc = FALSE) {
  
  # function that prepares the data for ATT/ATE analysis
  
  # proj.name: CARB project to assess
  # treat.type: the treatment cell type to use
  # control.type: control group cell type to use
  
  ##---Fixed characteristics--##
  
  point.of.treat <- proj.years[proj.years$project.name == proj.name, ]$DATE.first
  #point of treatment
  
  nlcd.to.look.at <- nlcd.years[nlcd.years <= (year(point.of.treat)-4)][length(nlcd.years[nlcd.years <= (year(point.of.treat)-4)])]
  # the closest NLCD cover
  
  if (year(point.of.treat)==2004) {
    nlcd.to.look.at <- 2001
  }
  
  nlcds.to.remove <- nlcd.years[!nlcd.years %in% nlcd.to.look.at]
  
  treated.cells <- crosswalk[crosswalk$proj.id==proj.name & crosswalk$type==treat.type, ]
  control.cells <- crosswalk[crosswalk$proj.id==proj.name & crosswalk$type==control.type, ]
  
  treated.cellid <- unique.cells.data[match(treated.cells$cellid, unique.cells.data$cellID), "cellID"]
  control.cellid <- unique.cells.data[match(control.cells$cellid, unique.cells.data$cellID), "cellID"]
  
  treated.fix <- unique.cells.data[unique.cells.data$cellID %in% treated.cellid, ] #fixed characteristics for treated cells
  control.fix <- unique.cells.data[unique.cells.data$cellID %in% control.cellid, ] #fixed characteristics for control cells
  
  treated.fix <- treated.fix[,-which(colnames(treated.fix) %in% paste0("nlcd_", nlcds.to.remove))]
  control.fix <- control.fix[,-which(colnames(control.fix) %in% paste0("nlcd_", nlcds.to.remove))]
  
  colnames(treated.fix)[which(colnames(treated.fix)==paste0("nlcd_", nlcd.to.look.at))] <- "nlcd"
  colnames(control.fix)[which(colnames(control.fix)==paste0("nlcd_", nlcd.to.look.at))] <- "nlcd"
  
  ##--Time-varying characteristics--##
  ## Default level of aggregation: month 
  # we do not do na.rm, as NA may not be random
  
  years.minus2and3 <- c(year(point.of.treat)-2, year(point.of.treat)-3)
  years.minus4and5 <- c(year(point.of.treat)-4, year(point.of.treat)-5)
  
  if (year(point.of.treat) == 2004) { #case for CAFR0001
    years.minus4and5[2] <- 2000
    #turn it into 2000
  }
  
  ndvi.dat.tminus2 <- get(paste0("ndvi.evi.summ.dat.", years.minus2and3[1]))
  ndvi.dat.tminus3 <- get(paste0("ndvi.evi.summ.dat.", years.minus2and3[2]))
  ndvi.dat.tminus4 <- get(paste0("ndvi.evi.summ.dat.", years.minus4and5[1]))
  ndvi.dat.tminus5 <- get(paste0("ndvi.evi.summ.dat.", years.minus4and5[2]))
  
  ndvi.dat.for.preceding.tminus2to5 <- rbind(ndvi.dat.tminus2[ndvi.dat.tminus2$cellID %in% c(treated.cellid, control.cellid),],
                                             ndvi.dat.tminus3[ndvi.dat.tminus3$cellID %in% c(treated.cellid, control.cellid),],
                                             ndvi.dat.tminus4[ndvi.dat.tminus4$cellID %in% c(treated.cellid, control.cellid),],
                                             ndvi.dat.tminus5[ndvi.dat.tminus5$cellID %in% c(treated.cellid, control.cellid),])
  
  ndvi.dat.for.preceding.tminus2to5$YEAR <- abs(ndvi.dat.for.preceding.tminus2to5$YEAR-year(point.of.treat))
  # example: if project started on 2004-05-01, then 2002 (2), 2001 (3), 2000 (4), 1999 (5)
  
  biomass.tminus0 <- select(subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid)),
                            c("cellID", paste0("biomass_", year(point.of.treat))))
  
  colnames(biomass.tminus0) <- c("cellID", "biomass_0")
  
  biomass.tminus0$biomass_tminus1 <- subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid))[[paste0("biomass_", year(point.of.treat)-1)]]
  biomass.tminus0$biomass_tminus2 <- subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid))[[paste0("biomass_", year(point.of.treat)-2)]]
  biomass.tminus0$biomass_tminus3 <- subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid))[[paste0("biomass_", year(point.of.treat)-3)]]
  biomass.tminus0$biomass_tminus4 <- subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid))[[paste0("biomass_", year(point.of.treat)-4)]]
  biomass.tminus0$biomass_tminus5 <- subset(biomass.extracted.all, cellID %in% c(treated.cellid, control.cellid))[[paste0("biomass_", year(point.of.treat)-5)]]

  
  if (nrow(ndvi.dat.for.preceding.tminus2to5)==0) {
    #do nothing
  } else {
    if (flag.qc) { #if we are removing QC-flagged
      
      ndvi.dat.for.preceding.tminus2to5 %>% 
        filter(cellID %in% treated.cellid & flagged) %>% 
        select(!flagged)-> ndvi.evi.treat
      
      ndvi.dat.for.preceding.tminus2to5 %>% 
        filter(cellID %in% control.cellid & flagged) %>% 
        select(!flagged)-> ndvi.evi.control
      
      ndvi.evi.treat <- ndvi.evi.treat[!duplicated(ndvi.evi.treat[c("cellID", "YEAR")]),]
      ndvi.evi.control <- ndvi.evi.control[!duplicated(ndvi.evi.control[c("cellID", "YEAR")]),]
      
    } else { #we are keeping QC-flagged entries
      
      ndvi.dat.for.preceding.tminus2to5 %>% 
        filter(cellID %in% treated.cellid & !flagged) %>% 
        select(!flagged) -> ndvi.evi.treat
      
      ndvi.dat.for.preceding.tminus2to5 %>% 
        filter(cellID %in% control.cellid & !flagged) %>% 
        select(!flagged)-> ndvi.evi.control
      
      ndvi.evi.treat <- ndvi.evi.treat[!duplicated(ndvi.evi.treat[c("cellID", "YEAR")]),]
      ndvi.evi.control <- ndvi.evi.control[!duplicated(ndvi.evi.control[c("cellID", "YEAR")]),]
    }
    
    ndvi.evi.control <- dcast(melt(ndvi.evi.control, id.vars=c("cellID", "YEAR")), cellID~variable+YEAR)
    
    if (nrow(ndvi.evi.treat)==0) {
      ndvi.evi.treat <- data.frame(matrix(ncol = ncol(ndvi.evi.control), nrow = 0))
      colnames(ndvi.evi.treat) <- colnames(ndvi.evi.control)
      ndvi.evi.treat$cellID <- as.character(ndvi.evi.treat$cellID)
    } else {
      ndvi.evi.treat <- dcast(melt(ndvi.evi.treat, id.vars=c("cellID", "YEAR")), cellID~variable+YEAR)  
    }
    
    # from long to wide
    
    biomass.treat <- subset(biomass.tminus0, cellID %in% treated.cellid)
    biomass.control <- subset(biomass.tminus0, cellID %in% control.cellid)
    
    biomass.control$treat <- 0
    biomass.treat$treat <- 1
    
    biomass.treat$type <- "treated"
    biomass.control$type <- "control"
    
    treated.fix$cellID <- as.character(treated.fix$cellID)
    control.fix$cellID <- as.character(control.fix$cellID)
    
    biomass.treat$cellID <- as.character(biomass.treat$cellID)
    biomass.control$cellID <- as.character(biomass.control$cellID)
    
    ##--Final dataset--#
    
    final.dataset <- bind_rows(left_join(left_join(control.fix, ndvi.evi.control, by = 'cellID'), biomass.control, by = "cellID"), 
                               left_join(left_join(treated.fix, ndvi.evi.treat, by = 'cellID'), biomass.treat, by = "cellID"))
    
    final.dataset$proj <- proj.name
    
    final.dataset %>% 
      select(-proj.from.crosswalk) %>% #this was used for sanity checks in the fixed characteristics dataset
      relocate(project.ID = proj) -> final.dataset
    
    return(final.dataset)
  }
}

proj.dats.for.prop.b2k <- data.maker(projs.to.look.at[1], "exact.treated", "control.b2k", flag.qc = F)
proj.dats.for.prop.b2k.qcT <- data.maker(projs.to.look.at[1], "exact.treated", "control.b2k", flag.qc = T)
proj.dats.for.prop.pmethod <- data.maker(projs.to.look.at[1], "exact.treated", "control.pmethod", flag.qc = F)
proj.dats.for.prop.b1k <- data.maker(projs.to.look.at[1], "exact.treated", "control.b1k", flag.qc = F)

for (i in 2:length(projs.to.look.at)) {
  proj.dats.for.prop.b2k <- bind_rows(proj.dats.for.prop.b2k,
                                      data.maker(projs.to.look.at[i], "exact.treated", "control.b2k", flag.qc = F))
  proj.dats.for.prop.b2k.qcT <- bind_rows(proj.dats.for.prop.b2k.qcT,
                                          data.maker(projs.to.look.at[i], "exact.treated", "control.b2k", flag.qc = T))
  proj.dats.for.prop.pmethod <- bind_rows(proj.dats.for.prop.pmethod,
                                          data.maker(projs.to.look.at[i], "exact.treated", "control.pmethod", flag.qc = F))
  proj.dats.for.prop.b1k <- bind_rows(proj.dats.for.prop.b1k,
                                      data.maker(projs.to.look.at[i], "exact.treated", "control.b1k", flag.qc = F))
  
  print(paste0("DONE WITH PROJECT ", i))
}

proj.dats.for.prop.b2k$nlcd <- as.factor(proj.dats.for.prop.b2k$nlcd)
proj.dats.for.prop.b2k.qcT$nlcd <- as.factor(proj.dats.for.prop.b2k.qcT$nlcd)
proj.dats.for.prop.pmethod$nlcd <- as.factor(proj.dats.for.prop.pmethod$nlcd)
proj.dats.for.prop.b1k$nlcd <- as.factor(proj.dats.for.prop.b1k$nlcd)

save(list=c('proj.dats.for.prop.b2k'),
     file = "Step2.5_ProjectSpecificData_b2k_updated.RData")

save(list=c('proj.dats.for.prop.b2k.qcT'),
     file = "Step2.5_ProjectSpecificData_b2kqcT_updated.RData")

save(list=c('proj.dats.for.prop.pmethod'),
     file = "Step2.5_ProjectSpecificData_pmethod_updated.RData")

save(list=c('proj.dats.for.prop.b1k'),
     file = "Step2.5_ProjectSpecificData_b1k_updated.RData")

rm(list=ls())

# Identifying treatment and non-forest pixels =====
treated.pixels <- unique(subset(proj.dats.for.prop.b1k,treat ==1)$cellID)
non.forest.pixels <- union(union(unique(subset(proj.dats.for.prop.b1k,!nlcd %in% c("41", "42", "43",
                                                                                   "51", "52", "90", "95"))$cellID),
                                 unique(subset(proj.dats.for.prop.b2k,!nlcd %in% c("41", "42", "43",
                                                                                   "51", "52", "90", "95"))$cellID)),
                           unique(subset(proj.dats.for.prop.pmethod,!nlcd %in% c("41", "42", "43",
                                                                                 "51", "52", "90", "95"))$cellID))

save(list = c("treated.pixels", "non.forest.pixels"),
     file = "Step2.5_NonForestAndTreatedPixels.RData")

# Identifying spatial clusters ======
load("Step2.5_ProjectSpecificData_b1k_updated.RData")

### Checking if the projects are disjoint (repeated) =====
arboc.map.projs.simpl <- terra::vect("ARBOC Map/OffsetsMapForestSimplified.shp")

terra::plot(arboc.map.projs.simpl[20,]) #example of project with 'clustered' areas

sample.buffer <- terra::buffer(arboc.map.projs.simpl[20,], 10^4)

terra::plot(sample.buffer) #the buffers do not intersect

disagged.projs <- terra::disagg(sample.buffer)

clustered.projs <- c()

arboc.map.projs.simpl.buff <- terra::buffer(arboc.map.projs.simpl, 10^4)

for (i in 1:nrow(arboc.map.projs.simpl)) {
  
  buffered.proj <- arboc.map.projs.simpl.buff[i,]
  #buffer the project
  
  buffered.proj.disagg <- terra::disagg(buffered.proj)
  
  if (nrow(buffered.proj.disagg)==1) {
    print(paste0("NO CLUSTERING FOR PROJECT NO. ", i, "!! MOVING ON :)"))
    
  } else {
    
    print(paste0('CLUSTERS FOUND FOR PROJECT NO. ', i, "!!!"))
    
    clustered.projs <- c(clustered.projs, arboc.map.projs.simpl[i,]$ARB_Projec)
    
    disagged.projs <- c(disagged.projs,
                        buffered.proj.disagg)
    
  }
  
  
}

disagged.projs.list <- c()

for (i in 1:length(disagged.projs)) {
  disagged.projs.list <- c(disagged.projs.list,
                           unique(disagged.projs[[i]]$ARB_Projec))
}

projs.in.dat <- unique(proj.dats.for.prop.b1k$proj)

spatial.cluster.list <- list()

for (i in 1:length(projs.in.dat)) {
  
  proj.to.assess <- projs.in.dat[i]
  
  if (proj.to.assess %in% disagged.projs.list) {
    print(paste0("THIS PROJECT HAS DISJOINTED BORDERS: WILLL BE ASSESSED BELOW"))
    
    next
  }
  
  proj.df <- proj.dats.for.prop.b1k[proj.dats.for.prop.b1k$project.ID==proj.to.assess,]
  
  proj.vect <- terra::vect(proj.df, geom = c("lon", "lat"), crs = "epsg:4269", keepgeom = T)
  
  proj.bbox <- terra::ext(terra::buffer(proj.vect, width = 10000, capstyle = "square"))
  #buffer by 10km
  
  proj.bbox.grid.1km <- terra::as.polygons(terra::rast(extent = proj.bbox, resolution = 0.01, crs = terra::crs(proj.vect)))
  #0.01 degree: Around 1km
  
  proj.bbox.grid.25km <- terra::as.polygons(terra::rast(extent = proj.bbox, resolution = 0.025, crs = terra::crs(proj.vect)))
  
  proj.bbox.grid.5km <- terra::as.polygons(terra::rast(extent = proj.bbox, resolution = 0.05, crs = terra::crs(proj.vect)))
  
  proj.bbox.grid.1km.extracted <- terra::extract(proj.bbox.grid.1km, proj.vect)
  proj.bbox.grid.25km.extracted <- terra::extract(proj.bbox.grid.25km, proj.vect)
  proj.bbox.grid.5km.extracted <- terra::extract(proj.bbox.grid.5km, proj.vect)
  
  proj.df$cluster.1km <- paste0(proj.bbox.grid.1km.extracted[,2])
  proj.df$cluster.25km <- paste0(proj.bbox.grid.25km.extracted[,2])
  proj.df$cluster.5km <- paste0(proj.bbox.grid.5km.extracted[,2])
  
  spatial.cluster.list[[proj.to.assess]] <- proj.df
  
  print(paste0("DONE WITH NO. ", i))
  
  
}

for (i in 1:length(disagged.projs.list)) {
  
  #for projects with disjointed borders
  
  proj.to.assess <- disagged.projs.list[i]
  
  proj.df <- proj.dats.for.prop.b1k[proj.dats.for.prop.b1k$project.ID==proj.to.assess,]
  
  proj.vect <- terra::vect(proj.df, geom = c("lon", "lat"), crs = "epsg:4269", keepgeom = T)
  
  if (nrow(proj.df)==0)  {
    print(paste0("THIS PROJECT NOT IN THE DATASET (TOO RECENT)! MOVING ON..."))
    
    next
  }
  
  proj.disagged.vect <- disagged.projs[i]
  #the shapefile belonging to this project, disaggregated (if there are three disjointed buffers, then three)
  
  num.of.disjoints <- length(proj.disagged.vect)
  
  disjoint.df.list <- list()
  
  for (k in 1:num.of.disjoints) {
    
    proj.disagged.vect.for.this.cluster <- proj.disagged.vect[k]
    
    pixels.belonging.to.this.cluster <- which(!is.na(terra::extract(terra::buffer(proj.disagged.vect.for.this.cluster, 5000), proj.vect)$ARB_Projec))
    
    proj.df.belonging.to.this.cluster <- proj.df[pixels.belonging.to.this.cluster,]
    proj.vect.belonging.to.this.cluster <- proj.vect[pixels.belonging.to.this.cluster,]
    
    if (nrow(proj.vect.belonging.to.this.cluster[proj.vect.belonging.to.this.cluster$treat == 1,]) == 0) {
      print(paste0('STOP!!!! THERE IS NOT ENOUGH TREATMENTS HERE'))
      
      break
    }
    
    proj.bbox.for.this.cluster <- terra::ext(terra::buffer(proj.vect.belonging.to.this.cluster, width = 10000, capstyle = "square"))
    #buffer by 10000m
    
    proj.bbox.grid.1km <- terra::as.polygons(terra::rast(extent = proj.bbox.for.this.cluster, resolution = 0.01, crs = terra::crs(proj.bbox)))
    #0.01 degree: Around 1km
    
    proj.bbox.grid.25km <- terra::as.polygons(terra::rast(extent = proj.bbox.for.this.cluster, resolution = 0.025, crs = terra::crs(proj.bbox)))
    
    proj.bbox.grid.5km <- terra::as.polygons(terra::rast(extent = proj.bbox.for.this.cluster, resolution = 0.05, crs = terra::crs(proj.bbox)))
    
    proj.bbox.grid.1km.extracted <- terra::extract(proj.bbox.grid.1km, proj.vect.belonging.to.this.cluster)
    proj.bbox.grid.25km.extracted <- terra::extract(proj.bbox.grid.25km, proj.vect.belonging.to.this.cluster)
    proj.bbox.grid.5km.extracted <- terra::extract(proj.bbox.grid.5km, proj.vect.belonging.to.this.cluster)
    
    proj.df.belonging.to.this.cluster$cluster.1km <- paste0(LETTERS[k], "-", proj.bbox.grid.1km.extracted[,2])
    proj.df.belonging.to.this.cluster$cluster.25km <- paste0(LETTERS[k], "-", proj.bbox.grid.25km.extracted[,2])
    proj.df.belonging.to.this.cluster$cluster.5km <- paste0(LETTERS[k], "-", proj.bbox.grid.5km.extracted[,2])
    
    disjoint.df.list[[k]] <- proj.df.belonging.to.this.cluster
    
    
  }
  
  disjoint.df <- bind_rows(disjoint.df.list)
  
  spatial.cluster.list[[proj.to.assess]] <- disjoint.df
  
  print(paste0("DONE WITH NO. ", i))
  
}

spatial.cluster <- bind_rows(spatial.cluster.list)

spatial.cluster$proj.and.cellID <- paste0(spatial.cluster$proj, "-", spatial.cluster$cellID)

save(list = "spatial.cluster",
     file = "Step2.5_SpatialClusters.RData")

### Adding clusters to the data set ======

load("Step2.5_SpatialClusters.RData")
load("Step2.5_ProjectSpecificData_b1k_updated.RData")
load("Step2.5_ProjectSpecificData_b2k_updated.RData")
load("Step2.5_ProjectSpecificData_b2kqcT_updated.RData")

for (a in c('proj.dats.for.prop.b1k', 'proj.dats.for.prop.b2k', 
            'proj.dats.for.prop.b2k.qcT')) {
  
  dat <- get(a)
  
  dat$proj.and.cellID <- paste0(dat$proj, "-", dat$cellID)
  
  dat %>% 
    select(!contains("cluster")) -> dat
  
  dat$cluster.1km <- spatial.cluster[match(dat$proj.and.cellID, spatial.cluster$proj.and.cellID),]$cluster.1km
  dat$cluster.25km <- spatial.cluster[match(dat$proj.and.cellID, spatial.cluster$proj.and.cellID),]$cluster.25km
  dat$cluster.5km <- spatial.cluster[match(dat$proj.and.cellID, spatial.cluster$proj.and.cellID),]$cluster.5km
  
  dat <- dat[,-which(colnames(dat)=="proj.and.cellID")]
  
  assign(a, dat)
  
  rm(dat)
  
}

save(list=c('proj.dats.for.prop.b2k'),
     file = "Step2.5_ProjectSpecificData_b2k_updated.RData")

save(list=c('proj.dats.for.prop.b2k.qcT'),
     file = "Step2.5_ProjectSpecificData_b2kqcT_updated.RData")

save(list=c('proj.dats.for.prop.b1k'),
     file = "Step2.5_ProjectSpecificData_b1k_updated.RData")


proj.dats.for.prop.b1k %>% 
  filter(project.ID=="CAFR5089") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(cluster.25km))) + scale_color_paletteer_d("palettesForR::Named") + geom_sf() + theme_bw() + theme(legend.position = "none")

proj.dats.for.prop.b1k %>% 
  filter(project.ID=="CAFR5089") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(cluster.5km))) + scale_color_paletteer_d("palettesForR::Named") + geom_sf() + theme_bw() + theme(legend.position = "none")

proj.dats.for.prop.b1k %>% 
  filter(project.ID=="CAFR0030") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(cluster.25km))) + scale_color_paletteer_d("palettesForR::Named") + geom_sf() + theme_bw() + theme(legend.position = "none")

proj.dats.for.prop.b1k %>% 
  filter(project.ID=="CAFR0030") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(cluster.5km))) + scale_color_paletteer_d("palettesForR::Named") + geom_sf() + theme_bw() + theme(legend.position = "none")


proj.dats.for.prop.b2k %>% 
  filter(project.ID=="CAFR0080") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(quadrants.treatBRUTEMED))) + scale_color_paletteer_d("khroma::okabeito") + geom_sf() + theme_bw()

proj.dats.for.prop.b2k %>% 
  filter(project.ID=="CAFR0080") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(octants.treatBRUTEMED))) + scale_color_paletteer_d("khroma::okabeito") + geom_sf() + theme_bw()

proj.dats.for.prop.b2k %>% 
  filter(project.ID=="CAFR0030") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(octants.treatBRUTEMED))) + scale_color_paletteer_d("khroma::okabeito") + geom_sf() + theme_bw()

proj.dats.for.prop.b2k %>% 
  filter(project.ID=="CAFR0030") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(quadrants.treatBRUTEMED))) + scale_color_paletteer_d("khroma::okabeito") + geom_sf() + theme_bw()

proj.dats.for.prop.b2k %>% 
  filter(project.ID=="CAFR0028") %>%  
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = factor(octants.treatBRUTEMED))) + scale_color_paletteer_d("khroma::okabeito") + geom_sf() + theme_bw()
