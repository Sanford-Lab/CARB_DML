library(tidyverse);library(lubridate);library(exactextractr);library(sf)
library(tictoc);library(splines);library(paletteer)

cell.df <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv") # crosswalk file
cell.df.unique <- cell.df[!duplicated(cell.df$cellid),]
cell.df.unique.shp <- terra::vect("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_ALL_NDVI250m.shp")

# cell.df.unique.buffered.shp <- terra::vect("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m.shp")

cell.df.unique.buffered.shp.sf <- st_read("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m.shp")

sample.gee.dat <- read_csv("F:/OFFSET_250m/ClmData/clm_DEM_250m_1.csv")

for (i in 2:41) {
  new.dat <- read_csv(paste0('F:/OFFSET_250m/ClmData/clm_DEM_250m_', i, ".csv"))
  
  if (all(new.dat$FID == (1:nrow(new.dat)-1))) {
    #check if GEE has messed up the FID orders or not
    sample.gee.dat <- bind_rows(sample.gee.dat,
                                new.dat)
  } else {
    paste0('ERROR!! SOMETHING WRONG WITH BULK NO. ', i)
    
    break
  }
  
  
}

# see that the coordinates or the same
# i.e., sample.gee.dat[1,".geo"] is the same as the coordinates of terra::geom(cell.df.unique.shp[1,])

unique.cells.data <- data.frame(cellID = cell.df.unique$cellid)

number2binary <- function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

# Covariates: Non-climate, fixed characteristics =======

## Soil and DEM =======

fixed.chr <- list.files("F:/Offset_250m/ClmData")[grepl('clm', list.files("F:/Offset_250m/ClmData"))]

fixed.chr <- gsub(paste(paste0("\\_", 1:41, ".csv"), collapse = "|"),
                  "",
                  fixed.chr)

fixed.chr <- fixed.chr[!duplicated(fixed.chr)]

fixed.chr <- fixed.chr[!fixed.chr %in% c('clm_nlcd_2001_250m', 'clm_nlcd_2004_250m',
                                         'clm_nlcd_2006_250m', 'clm_nlcd_2008_250m',
                                         'clm_nlcd_2011_250m', 'clm_nlcd_2013_250m',
                                         'clm_nlcd_2016_250m', 'clm_nlcd_2019_250m')]

for (datname in fixed.chr) {
  
  dat <- read_csv(paste0("F:/Offset_250m/ClmData/", datname, "_1.csv"))
  
  
  for (i in 2:41) {
    new.dat <- read_csv(paste0("F:/Offset_250m/ClmData/", datname, "_", i, ".csv"))
    
    if (all(new.dat$FID == (1:nrow(new.dat)-1))) {
      #check if GEE has messed up the FID orders or not
      
      dat <- bind_rows(dat, new.dat)
      
    } else {
      paste0('ERROR!! SOMETHING WRONG WITH BULK NO. ', i)
      
      break
    }
  }
  
  if (nrow(dat) != nrow(cell.df.unique)) {
    
    print(paste0("ERROR!!! SOMETHING WRONG WITH FILE ", datname))
    
    #these fixed characteristics are supposed to have the same nrow as the shapefile that we used for GEE query
    
    break
    
  }
  
  col.position <- which(grepl("mean|b100|mode", colnames(dat)))
  #use 1m soil property and mean DEM
  
  name.of.this.var <- gsub("_250m", "", datname)
  
  unique.cells.data[[name.of.this.var]] <- dat[[col.position]]
  
  print(paste0("DONE WITH DATA ", datname, "! FIRST QUERY: ",
               dat[1,col.position]))
  
}

rm(dat)

table(unique.cells.data$clm_ned_lf)
#see that this column is a factor column

unique.cells.data$clm_ned_lf <- as.character(unique.cells.data$clm_ned_lf)

## NLCD data ======

nlcd.legend <- list("11" = "Open.Water", "12" = "P.Ice.Snow",
                    "21" = "Developed.Open", "22" = "Developed.LowI",
                    "23" = "Developed.MedI", "24" = "Developed.HighI",
                    "31" = "Barren.Land", "41" = "Decidious.Forest", "42" = "Evergreen.Forest",
                    "43" = "Mixed.Forest", "51" = "Dwarf.Shrub", "52" = "Shrub.scrub",
                    "71" = "Grassland", "72" = "Sedge", "73" = "Lichens", "74" = "Moss",
                    "81" = "Pasture", "82" = "Cultivated.crops", "90" = "Woody.wetlands", 
                    "95" = "Emergent.Herb.Wetlands",
                    "99" = "Unknown")

for (yyyy in c(2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019)) {
  
  nlcd.dat <- read_csv(paste0("F:/Offset_250m/ClmData/clm_nlcd_", yyyy, "_250m_1.csv"))
  
  for (i in 2:41) {
    new.dat <- read_csv(paste0("F:/Offset_250m/ClmData/clm_nlcd_", yyyy, "_250m_", i, ".csv"))
    
    if (all(new.dat$FID == (1:nrow(new.dat)-1))) {
      #check if GEE has messed up the FID orders or not
      
      nlcd.dat <- bind_rows(nlcd.dat, new.dat)
      
    } else {
      paste0('ERROR!! SOMETHING WRONG WITH BULK NO. ', i)
      
      break
    }
    
  }
  
  unique.cells.data[[paste0("nlcd_", yyyy)]] <- as.character(nlcd.dat$mode)
  
  # unique.cells.data[[paste0("nlcd_expl_", yyyy)]] <- unname(unlist(nlcd.legend[match(unique.cells.data$nlcd, names(nlcd.legend),
  #                                                                                    nomatch = 21)])) 
  #return 99 if does not match to anything (because 99 is the 21st entry)
  
  
}

unique.cells.data$lon <- cell.df.unique[match(unique.cells.data$cellID, cell.df.unique$cellid),]$lon
unique.cells.data$lat <- cell.df.unique[match(unique.cells.data$cellID, cell.df.unique$cellid),]$lat
unique.cells.data$proj.from.crosswalk <- cell.df.unique[match(unique.cells.data$cellID, cell.df.unique$cellid),]$proj.id

### sanity check for NLCD data ======
cafr0001.fixed <- unique.cells.data[unique.cells.data$proj.from.crosswalk=="CAFR0001",]

cafr0001.fixed %>% st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = clm_DEM)) + geom_sf()
cafr0001.fixed %>% st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = clm_soilORGC)) + geom_sf()
cafr0001.fixed %>% st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = nlcd_2001)) + geom_sf()
#continuous/smooth transitions

save(list = "unique.cells.data", file = 'Step2_1_Fixed characteristics_NDVI250m.RData')

## Tree characteristics =====

tree.type.dat <- terra::rast("E:/OFFSET_DATA/Forest_Type/conus_foresttype.img")
tree.type.dat <- terra::project(tree.type.dat, "epsg:4326")

terra::activeCat(tree.type.dat) <- "ForestType" # Activate "ForestType" for this column

treetype.extracted <- exact_extract(tree.type.dat, cell.df.unique.buffered.shp.sf[1,], fun = "mode")

nrow(cell.df.unique) # the number of FIDs

tic()

for (p in 1:(nrow(cell.df.unique) %/% 5000)) {
  
  treetype.extracted.for.5000 <- exact_extract(tree.type.dat, 
                                               cell.df.unique.buffered.shp.sf[(5000*(p-1) +2):(5000*p+1), ], 
                                              fun = "mode", max_cells_in_memory = 1e+08)
  
  treetype.extracted <- c(treetype.extracted,
                          treetype.extracted.for.5000)
  
  print(paste0("DONE WITH ", p, "TH BUNCH!"))
  
  toc()
  tic()
  
  if (p %% 100 == 0) { # every 500,000th cellID
    save(list = "treetype.extracted",
         file = paste0("F:/Offset_250m/Step2_TREETYPE_UNTIL", p*5000, "_UPDATED.RData"))
  }
  
}

tree.extracted.rest <- exact_extract(tree.type.dat, 
                                     cell.df.unique.buffered.shp.sf[(5000*(p) +2):nrow(cell.df.unique.buffered.shp.sf),], 
                                        fun = "mode")

treetype.extracted <- c(treetype.extracted, tree.extracted.rest)

unique.cells.data$forest.type <- as.character(treetype.extracted)

unique.cells.data$forest.group <- as.factor(paste0(substr(unique.cells.data$forest.type,1,2), "0"))
#larger category than forest.type

save(list = "unique.cells.data", file = 'Step2_1_Fixed characteristics_NDVI250m.RData')

## Forest ownership (2023.06.02) =======

load("Step2_1_Fixed characteristics_NDVI250m.RData")

fownership <- terra::rast("E:/OFFSET_DATA_ALL/FOREST_OWERNSHIP_MAP/year2017/Data/forest_own1/forest_own1.tif")

cell.df.unique.buffered.shp.sf <- st_transform(cell.df.unique.buffered.shp.sf, terra::crs(fownership))
# turn into fownership coors

nrow(cell.df.unique.buffered.shp) # the number of FIDs

tic()

fownership.extracted <- exact_extract(fownership, cell.df.unique.buffered.shp.sf[1:10**5,], fun = "mode")

for (i in 1:40) {
  fownership.extracted <- c(fownership.extracted,
                            exact_extract(fownership, cell.df.unique.buffered.shp.sf[(i*(10**5)+1):((i+1)*(10**5)),], fun = "mode"))  
  
  print(paste0("DONE WITH BUNCH NO.", i))
}

fownership.extracted <- c(fownership.extracted,
                          exact_extract(fownership, cell.df.unique.buffered.shp.sf[(41*(10**5)+1):nrow(cell.df.unique.buffered.shp.sf),], 
                                        fun = "mode"))  

unique.cells.data$fownership <- fownership.extracted

unique.cells.data$fownership <- as.character(unique.cells.data$fownership)

unique.cells.data %>% filter(proj.from.crosswalk=="CAFR0002") %>%  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = fownership)) + geom_sf()
#sanity check

# save(list = "unique.cells.data", file = 'Step2_1_Fixed characteristics_NLCDUPDATED_ALL.RData')
# this was an error that caused the addition of the "adding fownership" section in R_step2.5_projectspecificdata.R
# rectified as follows

save(list = "unique.cells.data", file = 'Step2_1_Fixed characteristics_NDVI250m.RData')

# Outcome variable: biomass ======

# biomass.dat <- terra::rast("E:/OFFSET_DATA_ALL/Biomass/lt-stem_biomass_nbcd_v0.1_conus_median_2000_2017.tif")

biomass.dat <- terra::rast("D:/OFFSET_DATA_ALL/Biomass/lt-stem_biomass_nbcd_v0.1_conus_median_2000_2017.tif")

cell.df.unique.buffered.shp.sf <- st_transform(cell.df.unique.buffered.shp.sf, terra::crs(biomass.dat))

biomass.extracted <- exact_extract(biomass.dat, cell.df.unique.buffered.shp.sf[1,], fun = "mean")

nrow(cell.df.unique) # the number of FIDs

tic()

for (p in 1:(4157943 %/% 10000)) {
  
  biomass.extracted.for.10000 <- exact_extract(biomass.dat, 
                                              cell.df.unique.buffered.shp.sf[(10000*(p-1) +2):(10000*p+1) ,], 
                                              fun = "mean", max_cells_in_memory = 1e+08)
  
  biomass.extracted <- as.data.frame(data.table::rbindlist(list(biomass.extracted, biomass.extracted.for.10000)))
  
  print(paste0("DONE WITH ", p, "TH BUNCH!"))
  
  toc()
  tic()
  
  if (p %% 10 == 0) { # every 100,000th FID
    save(list = "biomass.extracted",
         file = paste0("Step2_BIOMASS_UNTIL", p*10000, "_UPDATED.RData"))
  }
  
}

biomass.extracted.rest <- exact_extract(biomass.dat, 
                                        cell.df.unique.buffered.shp.sf[(10000*p +2):nrow(cell.df.unique.buffered.shp.sf),], 
                                        fun = "mean")

biomass.extracted.all <- as.data.frame(data.table::rbindlist(list(biomass.extracted, biomass.extracted.rest)))

biomass.extracted.all$cellID <- cell.df.unique$cellid

save(list = "biomass.extracted.all",
     file = paste0("Step2_BIOMASS_NDVI250m.RData"))

load("Step2_BIOMASS_ALL.RData")

#Outcome variable: Avoided conversion (for AC projects) =====

load("STEP1_AC_LONLATS_NDVI250m.RData")

projs.in.ac <- unique(cells.in.refor.and.ac.projs.unique$proj.id)

yyyy <- 2000

lcms.dat.for.this.yyyy <- read_csv(paste0('F:/Offset_250m/ACData/lcms_', yyyy, '_share_1.csv'))  

for (i in 2:8) {
  lcms.dat.for.this.yyyy <- rbind(lcms.dat.for.this.yyyy,
                                  read_csv(paste0('F:/Offset_250m/ACData/lcms_', yyyy, "_share_", i, '.csv')))
}

lcms.dat.ratio <- data.frame(cellID = lcms.dat.for.this.yyyy$cellID)

lcms.dat.ratio[["lcms_ratio_2000"]] <- lcms.dat.for.this.yyyy$share

for (yyyy in 2001:2022) {
  
  lcms.dat.for.this.yyyy <- read_csv(paste0('F:/Offset_250m/ACData/lcms_', yyyy, '_share_1.csv'))  
  
  for (i in 2:8) {
    lcms.dat.for.this.yyyy <- rbind(lcms.dat.for.this.yyyy,
                                    read_csv(paste0('F:/Offset_250m/ACData/lcms_', yyyy, "_share_", i, '.csv')))
  }
  
  
  if (all(lcms.dat.ratio$cellID == lcms.dat.for.this.yyyy$cellID)) {
    lcms.dat.ratio[[paste0("lcms_ratio_", yyyy)]] <- lcms.dat.for.this.yyyy$share
  } else {
    print(paste0('ERROR!!! SOMETHING WRONG WITH DATA IN YEAR ', yyyy))
    
    break
  }
}

cells.in.refor.and.ac.projs.unique %>% 
  rename(cellID = cellid) -> cells.in.refor.and.ac.projs.unique

cells.in.refor.and.ac.projs.unique <- left_join(cells.in.refor.and.ac.projs.unique,
                                                lcms.dat.ratio, 
                                                by = "cellID")

cells.in.refor.and.ac.projs.unique %>% 
  select(-c("...1")) %>% 
  rename(lon.lcms = lon, lat.lcms = lat)-> cells.in.refor.and.ac.projs.unique

cells.in.refor.and.ac.projs.unique %>% 
  select(-c(type, vect.id, proj.id)) -> cells.in.refor.and.ac.projs.unique

save(list = c("cells.in.refor.and.ac.projs.unique", "projs.in.ac"),
     file = "Step2_AC_LCMS.RData")

# Past and present NDVI trends =======

cell.df <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv")
#crosswalk file

cell.df.split <- split(cell.df, cell.df$proj.id)

rm(cell.df)

arboc.dat <- readxl::read_excel("OFFSET_ISSUANCE DATA.xlsx", sheet = "ARB Offset Credit Issuance")

arboc.dat$project.name <-str_split_fixed(arboc.dat$`CARB Issuance ID`, "-",2)[,1] 

arboc.dat %>% 
  arrange(project.name, desc(`Reporting Period End Date`)) %>% 
  filter(!duplicated(project.name) & `Project Type` =='Forest') -> arboc.dat.last

arboc.dat %>%
  arrange(project.name,`Reporting Period Start Date`) %>% 
  filter(!duplicated(project.name) & `Project Type` =='Forest') -> arboc.dat.first

arboc.dat.first$DATE.first <- as.Date(as.integer(arboc.dat.first$`Reporting Period Start Date`), origin = "1899-12-30")
arboc.dat.last$DATE.last <- as.Date(as.integer(arboc.dat.last$`Reporting Period End Date`), origin = "1899-12-30")

proj.years <- merge(dplyr::select(arboc.dat.first, project.name, DATE.first), 
                    dplyr::select(arboc.dat.last, project.name, DATE.last),
                    by = "project.name")

proj.years$START.YEAR <- year(proj.years$DATE.first) -2
proj.years$END.YEAR <- year(proj.years$DATE.last) +1

proj.years$START.YEAR.prec2 <- proj.years$START.YEAR -2
#2 years preceding (to assess pre-treatment trend)

proj.years$START.YEAR.prec2 <- proj.years$START.YEAR -2

cellids.for.each.proj <- list()

for (proj in names(cell.df.split)) {
  
  cellids.for.this.proj <- unique(cell.df.split[[proj]]$cellid)
  
  cellids.for.each.proj[[proj]] <- cellids.for.this.proj
  
}

cellids.for.each.year <- list()

for (yyyy in 2000:2021) {
  
  projs.with.this.year <- proj.years[year(proj.years$DATE.first) %in% c(yyyy+2, yyyy+3, yyyy+4, yyyy+5), ]$project.name
  
  #Example: A project started on 2007-XX-XX
    #Treatment period: 2007
    #Baseline: 2006
    #Two preceding years: 2005, 2004
    #Two preceding years' two preceding years: 2003, 2004 / 2002, 2003 --> 2002, 2003, 2004
    #This project hence belongs to years 2002, 2003, 2004, 2005
  
  #Year 2005 would have projects that started on 
    #2007-XX-XX (2005, 2004, 2003, 2002)
    #2008-XX-XX (2006, 2005, 2004, 2003)
    #2009-XX-XX (2007, 2006, 2005, 2004)
    #2010-XX-XX (2008, 2007, 2006, 2005)
  
  cellids.this.year <- c()
  
  for (p in projs.with.this.year) {
    cellids.this.year <- c(cellids.this.year, cellids.for.each.proj[[p]])
  }
  
  cellids.this.year <- cellids.this.year[!duplicated(cellids.this.year)]
  
  cellids.for.each.year[[as.character(yyyy)]] <- cellids.this.year
  
}


ndvi.evi.summ.dat <- data.frame(cellID = character(),
                                YEAR = integer(),
                                flagged = logical(),
                                NDVI.max = numeric(),
                                EVI.max = numeric(),
                                NDVI.min = numeric(),
                                EVI.min = numeric(),
                                NDVI.max.doy = integer(),
                                NDVI.min.date = integer(),
                                EVI.max.date = integer(),
                                EVI.min.date = integer(),
                                ndvi.ns.coef0 = numeric(),
                                ndvi.ns.coef1 = numeric(),
                                ndvi.ns.coef2 = numeric(),
                                ndvi.ns.coef3 = numeric(),
                                evi.ns.coef0 = numeric(),
                                evi.ns.coef1 = numeric(),
                                evi.ns.coef2 = numeric(),
                                evi.ns.coef3 = numeric())

ndvi.evi.summ.giver <- function(data.to.use, flagged.tf) {
  
  cellid <- unique(data.to.use$cellID)
  year <- unique(data.to.use$YEAR)
  ndvi.max.where <- max(which(data.to.use$NDVI == max(data.to.use$NDVI)))
  ndvi.min.where <- min(which(data.to.use$NDVI == min(data.to.use$NDVI)))
  evi.min.where <- min(which(data.to.use$EVI == min(data.to.use$EVI)))
  evi.max.where <- max(which(data.to.use$EVI == max(data.to.use$EVI)))
  
  ndvi.spline.fitted <- lm(NDVI~ns(doy, 3), data = data.to.use)
  evi.spline.fitted <- lm(EVI~ns(doy, 3), data = data.to.use)
  
  return(data.frame(cellID = cellid,
                    YEAR = year,
                    flagged = flagged.tf,
                    NDVI.max = data.to.use[ndvi.max.where,]$NDVI,
                    EVI.max = data.to.use[evi.max.where,]$EVI,
                    NDVI.min = data.to.use[ndvi.min.where,]$NDVI,
                    EVI.min = data.to.use[evi.min.where,]$EVI,
                    NDVI.max.doy = data.to.use[ndvi.max.where,]$doy,
                    NDVI.min.date = data.to.use[ndvi.min.where,]$doy,
                    EVI.max.date = data.to.use[evi.max.where,]$doy,
                    EVI.min.date = data.to.use[evi.min.where,]$doy,
                    ndvi.ns.coef0 = unname(coefficients(ndvi.spline.fitted)[1]),
                    ndvi.ns.coef1 = unname(coefficients(ndvi.spline.fitted)[2]),
                    ndvi.ns.coef2 = unname(coefficients(ndvi.spline.fitted)[3]),
                    ndvi.ns.coef3 = unname(coefficients(ndvi.spline.fitted)[4]),
                    evi.ns.coef0 = unname(coefficients(evi.spline.fitted)[1]),
                    evi.ns.coef1 = unname(coefficients(evi.spline.fitted)[2]),
                    evi.ns.coef2 = unname(coefficients(evi.spline.fitted)[3]),
                    evi.ns.coef3 = unname(coefficients(evi.spline.fitted)[4])))
  
}

cellids.wo.ndvi <- list()

save(list= c("cellids.for.each.proj", "cellids.for.each.year", "proj.years"),
     file = "Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")

load('Step2_1_Fixed characteristics_NDVI250m.RData')

for (yyyy in 2014:2019) {
  
  print(paste0("Staring year ", yyyy))
  
  cellids.to.look.at <- cellids.for.each.year[[as.character(yyyy)]]
  
  if (length(cellids.to.look.at)==0) {
    next
  } else {
    print(paste0("ROW COUNT FOR THIS YEAR: ", length(cellids.to.look.at)))
  }
  
  dat <- read_csv(paste0("F:/Offset_250m/NdviEvi/", paste0("ndvi_evi_", yyyy, "_ALL_1.csv")))
  
  cellids.for.this.bunch <- unique.cells.data$cellID[1:10^5]
  
  dat$cellID <- cellids.for.this.bunch[dat$FID+1]
  
  dat <- subset(dat, cellID %in% cellids.to.look.at)
  
  for (i in 2:41) {
    
    new.dat <- read_csv(paste0("F:/Offset_250m/NdviEvi/", paste0("ndvi_evi_", yyyy, "_ALL_", i, ".csv")))
    
    if (i %in% 2:40) {
      cellids.for.this.bunch <- unique.cells.data$cellID[((10^5)*(i-1)+1) : ((10^5)*i)]
    } else if (i==41)  {
      cellids.for.this.bunch <- unique.cells.data$cellID[(10^5*41+1):nrow(unique.cells.data)]
    }
    
    if (all(unique(new.dat$FID) %in% c(1:nrow(new.dat)-1))   ){
      #check if GEE has not added random FIDs
      
      new.dat$cellID <- cellids.for.this.bunch[new.dat$FID+1]
      
      new.dat <- subset(new.dat, cellID %in% cellids.to.look.at)
      
      dat <- bind_rows(dat, new.dat)
      
    } else {
      paste0('ERROR!! SOMETHING WRONG WITH BULK NO. ', i)
      
      
      
      break
    }
  }
  
  print(paste0("DONE WITH BINDING ROWS FOR YEAR ", yyyy))
  
  cellids.wo.ndvi[[as.character(yyyy)]] <- cellids.to.look.at[!cellids.to.look.at %in% unique(dat$cellID)]
  
  assign(paste0('NDVIEVIdat.subsetted.', yyyy), dat)
  
  save(list = paste0('NDVIEVIdat.subsetted.', yyyy),
       file = paste0("F:/Offset_250m/NdviEvi/", "NdviEviSubsetted_", yyyy, ".RData"))
  #save the subsetted data
  
  rm(list = paste0('NDVIEVIdat.subsetted.', yyyy))
  
  dat$utc <- ymd(substr(dat$`system:index`,1,10), tz = "UTC")
  
  modis.qcs.in.this.dataset <- unique(dat$DetailedQA)
  modis.qcs.in.this.dataset.bit <- sapply(modis.qcs.in.this.dataset, 
                                          number2binary, noBits = 16) 
  #QCs turned into binary bits
  
  modis.qcs.list <- vector("list", length = length(modis.qcs.in.this.dataset))
  names(modis.qcs.list) <- as.character(modis.qcs.in.this.dataset)
  
  #vi quality: 0 - good, 1 & 2 - bad, 3 - not produced
  
  modis.viusefull.list <- vector("list", length = length(modis.qcs.in.this.dataset))
  names(modis.viusefull.list) <- as.character(modis.qcs.in.this.dataset)
  
  #usefulness: 0 > 1 > 2 > 4 > 8 > 9 > 10 > 12 > 13 (12, 13 flagged) > 14 > 15
  
  for (i in 1:length(modis.qcs.in.this.dataset)) {
    
    first.two.digits <- modis.qcs.in.this.dataset.bit[1:2,i]
    next.four.digits <- modis.qcs.in.this.dataset.bit[3:6,i]
    
    modis.qcs.list[i] <- first.two.digits[1]*1 + first.two.digits[1]*2
    modis.viusefull.list[i] <- next.four.digits[1]*1 + next.four.digits[2]*2 + next.four.digits[3]*4 + next.four.digits[4]*8
    
  }
  
  dat$vi.qc.label <- unname(unlist(modis.qcs.list[match(dat$DetailedQA, names(modis.qcs.list))]))
  dat$vi.useful.label <- unname(unlist(modis.viusefull.list[match(dat$DetailedQA, names(modis.viusefull.list))]))
  
  dat$vi.qc.flagged <- !dat$vi.qc.label==0
  dat$vi.useful.flagged <- dat$vi.useful.label >= 12
  
  dat$YEAR <- year(dat$utc)
  dat$MONTH <- month(dat$utc)
  dat$DAY <- day(dat$utc)
  
  dat %>% 
    dplyr::select(cellID, YEAR, MONTH, DAY, utc, NDVI, EVI,
           vi.qc.label, vi.qc.flagged, vi.useful.label, vi.useful.flagged,
           sur_refl_b01, sur_refl_b02, sur_refl_b03, sur_refl_b07) %>% 
    dplyr::rename(UTC = utc) -> dat.final
  
  dat.final <- data.table::as.data.table(dat.final)
  
  cellids <- unique(dat$cellID)
  
  dat.final$doy <- yday(dat.final$UTC)
  
  tic()
  
  ndvi.evi.interim.dat <- data.frame(matrix(ncol = ncol(ndvi.evi.summ.dat), nrow = 0))
  colnames(ndvi.evi.interim.dat) <- colnames(ndvi.evi.summ.dat)
  
  for (i in 1:length(cellids)) {
  
    cellid <- cellids[i]
    
    dat.to.look <- dat.final[cellID == cellid,]
    # dat.to.look$doy <- yday(dat.to.look$UTC)
    
    dat.to.look.flagged <- dat.to.look[!dat.to.look$vi.qc.flagged,]
    
    if (nrow(dat.to.look.flagged)==0) {
      
    } else {
      ndvi.evi.interim.dat <- rbind(ndvi.evi.interim.dat, ndvi.evi.summ.giver(dat.to.look.flagged, T))  
    }
    
    ndvi.evi.interim.dat <- rbind(ndvi.evi.interim.dat, ndvi.evi.summ.giver(dat.to.look, F))
    
    
    if (i %% 1000==0) {
      
      ndvi.evi.summ.dat <- as.data.frame(data.table::rbindlist(list(ndvi.evi.summ.dat, ndvi.evi.interim.dat)))
      rm(ndvi.evi.interim.dat)
      
      ndvi.evi.interim.dat <- data.frame(matrix(ncol = ncol(ndvi.evi.summ.dat), nrow = 0))
      colnames(ndvi.evi.interim.dat) <- colnames(ndvi.evi.summ.dat)
      
      print(paste0("DONE WITH PIXEL NO. ", i))
      toc()
      tic()
    }
    
    
  }
  
  rm(dat)
  rm(dat.final)
  
  assign(paste0("ndvi.evi.summ.dat.", yyyy), ndvi.evi.summ.dat)
  
  save(list= c(paste0("ndvi.evi.summ.dat.", yyyy), "cellids.wo.ndvi"),
       file = paste0("Step2_2_NDVI_SUMMARIZED_", yyyy, ".RData"))
  
  rm(ndvi.evi.summ.dat)
  rm(list = paste0("ndvi.evi.summ.dat.", yyyy))
  
  ndvi.evi.summ.dat <- data.frame(cellID = character(),
                                  YEAR = integer(),
                                  flagged = logical(),
                                  NDVI.max = numeric(),
                                  EVI.max = numeric(),
                                  NDVI.min = numeric(),
                                  EVI.min = numeric(),
                                  NDVI.max.doy = integer(),
                                  NDVI.min.date = integer(),
                                  EVI.max.date = integer(),
                                  EVI.min.date = integer(),
                                  ndvi.ns.coef0 = numeric(),
                                  ndvi.ns.coef1 = numeric(),
                                  ndvi.ns.coef2 = numeric(),
                                  ndvi.ns.coef3 = numeric(),
                                  evi.ns.coef0 = numeric(),
                                  evi.ns.coef1 = numeric(),
                                  evi.ns.coef2 = numeric(),
                                  evi.ns.coef3 = numeric())
  
  print(paste0("DONE WITH YEAR ", yyyy))
}

#sanity check

load("Step2_2_NDVI_SUMMARIZED_2012.RData")

ndvi.evi.summ.dat.2012$proj.from.crosswalk <- unique.cells.data[match(ndvi.evi.summ.dat.2012$cellID, unique.cells.data$cellID),]$proj.from.crosswalk

ndvi.evi.summ.dat.2012$lon <- unique.cells.data[match(ndvi.evi.summ.dat.2012$cellID, unique.cells.data$cellID),]$lon
ndvi.evi.summ.dat.2012$lat <- unique.cells.data[match(ndvi.evi.summ.dat.2012$cellID, unique.cells.data$cellID),]$lat

cafr0041.ndvi <- ndvi.evi.summ.dat.2012[ndvi.evi.summ.dat.2012$proj.from.crosswalk=="CAFR0041",]

cafr0041.ndvi %>% st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = ndvi.ns.coef1)) + geom_sf()
cafr0041.ndvi %>% filter(ndvi.ns.coef0 > 8500) %>%  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% ggplot(aes(color = ndvi.ns.coef0)) + geom_sf()


# Distance to roads =====

### ========= Checking road.vect data (do not run) ========= ###

road.vect <- terra::vect('E:/OFFSET_DATA/North_American_Roads.shp')
road.vect <- terra::project(road.vect, "epsg:4326")

road.vect.df <- as_tibble(road.vect)

road.vect.df %>% 
  group_by(CLASS) %>% 
  summarise(mean(LANES), mean(SPEEDLIM)) 
# see that class 1 and 2 have high speed limit & large lanes
# class 1: Interstate, class 2: other freeways & expressways, 3: other principal arterial
# 4: Minor arterial, 5: Major collector, 6: Minor collector

# road.vect.class1 <- road.vect[road.vect$CLASS==1,] #just interstate
#road.vect.class123 <- road.vect[road.vect$CLASS %in% c(1,2,3),] 
# principal roads (interstate, other freeways, other principal arterial)
#road.vect.class456 <- road.vect[road.vect$CLASS %in% c(4,5,6),] 
# minor arterial, major collector, minor collector

#rm(road.vect)

road.vect.agged.alltheway <- terra::aggregate(road.vect)
#road.vect.class123.agged.alltheway <- terra::aggregate(road.vect.class123)
#road.vect.class456.agged.alltheway <- terra::aggregate(road.vect.class456)
#all fused into one line vector

### ========= Checking road.vect data (end) ========= ###

load("Step2_1_Fixed characteristics_NDVI250m.RData")

cell.df.unique.shp <- terra::vect("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_ALL_NDVI250m.shp")

arboc.map <- terra::vect("ARBOC Map/OffsetsMapForestSimplified.shp")
# cafr0030 <- arboc.map[arboc.map$ARB_Projec=="CAFR0030",]
# cafr0030 <- terra::project(cafr0030, "epsg:3857")

#road.vect.conus <- terra::vect("D:/OFFSET_DATA/Roads_CONUS.shp")
road.vect.conus <- terra::vect("E:/OFFSET_DATA/Roads_CONUS.shp")
#EPSG:4326 for these

#road.vect.class123.conus <- terra::crop(road.vect.class123.agged.alltheway, county.conus.shp)
#road.vect.class456.conus <- terra::crop(road.vect.class456.agged.alltheway, county.conus.shp)

road.vect.conus.reproj <- terra::project(road.vect.conus, "epsg:3857")

cell.df.unique.shp.reproj <- terra::project(cell.df.unique.shp, "epsg:3857")
#take the WGS84, to address a bug with terra distance function
#terra gives distance in degrees, not meters, when one of the two vectors are line/polygon AND the CRS are in lon/lat

load("Step2_2_cellIDs_YEAR_CROSSWALKS_UPDATED.RData")
load("Step3_REDUNDANT_PROJECTS.RData")

cell.df <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv")

cell.df.split <- split(cell.df, cell.df$proj.id)

rm(cell.df)

# cafr0030.dots <- terra::vect(unique.cells.data[unique.cells.data$cellID %in% cell.df.split$CAFR0030$cellid,],
#                              geom = c("lon", "lat"), crs = "epsg:4326")
# cafr0030.dots <- terra::project(cafr0030.dots, "epsg:3857")
# 
# #cafr0030.nearest <- terra::nearest(cafr0030.dots,road.vect.conus.reproj)
# 
# cafr0030.dist <- terra::distance(cafr0030.dots, road.vect.class123.cafr0030)
# cafr0030.dist.allroad <- terra::distance(cafr0030.dots, road.vect.cafr0030)
# 
# cafr0030.df <- unique.cells.data[unique.cells.data$cellID %in% cell.df.split$CAFR0030$cellid,]
# 
# cafr0030.df$distance.class123.road <- cafr0030.dist
# cafr0030.df$distance.to.road <- cafr0030.dist.allroad
# 
# cafr0030.df %>% 
#   st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
#   ggplot(aes(color = distance.class123.road)) + scale_color_paletteer_c("grDevices::Fall", direction = 1) + geom_sf() + theme_bw()
# 
# cafr0030.df %>% 
#   st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
#   ggplot(aes(color = clm_DEM_ALL)) + scale_color_paletteer_c("grDevices::Fall", direction = 1) + geom_sf() + theme_bw()
# 
# cafr0030.df %>% 
#   st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
#   ggplot(aes(color = distance.to.road)) + scale_color_paletteer_c("grDevices::Fall", direction = 1) + geom_sf() + theme_bw()

nonred.proj <- names(cell.df.split)[!names(cell.df.split) %in% red.df$redundant.proj]

cell.df.split.nonred <- list()

for (p in nonred.proj) {
  cell.df.split.nonred[[p]] <- cell.df.split[[p]]
}

dist.df <- data.frame(cellID = integer(),
                      proj = character(),
                      distance = numeric())

for (p in 1:length(names(cell.df.split.nonred))) {

  this.p <- names(cell.df.split.nonred)[p]

  dots.for.this.p <- terra::vect(unique.cells.data[unique.cells.data$cellID %in% cell.df.split.nonred[[this.p]]$cellid,],
                                 geom = c("lon", "lat"), crs = "epsg:4326")
  dots.for.this.p <- terra::project(dots.for.this.p, "epsg:3857")

  # dots.for.this.p.rast <- terra::rast(dots.for.this.p, resolution = 0.0025, vals = 0)

  this.p.shp <- arboc.map[arboc.map$ARB_Projec==this.p,] # in EPSG:4326
  this.p.shp.reproj <- terra::project(this.p.shp, "epsg:3857") #CRS in meters

  this.p.shp.reproj.simplified <- terra::simplifyGeom(this.p.shp.reproj, 1000)
  # #simplify to ease buffers (not really major constraint, as the sole purpose is to obtain a buffer & crop the road)
  # #the buffer is 100km, so an error of 1km permitted seems fine
  
  this.p.shp.simplified <- terra::simplifyGeom(this.p.shp, 0.01)
  #Also simplified to 1km
  
  if (nrow(dots.for.this.p) >= 10^5) {
    #if this project area is massive
    this.p.shp.buffered <- terra::buffer(this.p.shp.simplified, 10^5*0.5) 
  } else {
    this.p.shp.buffered <- terra::buffer(this.p.shp.simplified, 10^5)
  }

  print(paste0("BUFFERING COMPLETE!!"))

  road.vect.for.this.p <- terra::crop(road.vect.conus, this.p.shp.buffered)
  #use EPSG4326 data (road.vect.conus and this.p.shp)
  #because lon-lat CRS has faster buffer

  print(paste0("DONE WITH BUFFER-CROP OF PROJ NO. ", p))

  road.vect.for.this.p <- terra::project(road.vect.for.this.p, "epsg:3857")

  dist.allroad <- terra::distance(dots.for.this.p, road.vect.for.this.p)

  # dist.allroad <- terra::distance(dots.for.this.p, road.vect.conus.reproj)
  #
  dist.df <- rbind(dist.df,
                   data.frame(cellID = dots.for.this.p$cellID,
                              proj = this.p,
                              distance = dist.allroad))
  
  dist.df <- dist.df[!duplicated(dist.df$cellID),]

  print(paste0("DONE WITH PROJ NO. ", p))
  
  if (p %% 10 == 0) {
    save(list = "dist.df",
         file = paste0('Step2_DistanceDFUntil', p, ".RData"))
  }

}

#unique.cells.data <- unique.cells.data[,-which(grepl("distance|nearest", colnames(unique.cells.data)))]

#unique.cells.data <- unique.cells.data[,-which(grepl("distance|nearest", colnames(unique.cells.data)))]

unique.cells.data$distance.to.road <- dist.df[match(unique.cells.data$cellID, dist.df$cellID),]$distance

unique.cells.data %>% 
  filter(proj.from.crosswalk == "CAFR0030") %>% 
  st_as_sf(.,coords=c("lon","lat"),crs=4326,remove=F) %>% 
  ggplot(aes(color = distance.to.road)) + scale_color_paletteer_c("grDevices::Fall", direction = 1) + geom_sf() + theme_bw()
#behaving quite well

save(list = "unique.cells.data", file = 'Step2_1_Fixed characteristics_NDVI250m.RData')

# Forest types and forest groups =====

cell.df.unique.buffered.shp.sf <- st_read("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m.shp")

tree.type.dat <- terra::rast("D:/OFFSET_DATA/Forest_Type/conus_foresttype.img")
terra::activeCat(tree.type.dat) <- is.factor("ForestType") # Activate "ForestType" for this column

cell.df.unique.buffered.shp.sf <- st_transform(cell.df.unique.buffered.shp.sf, terra::crs(tree.type.dat))

tree.type.extracted <- exact_extract(tree.type.dat, cell.df.unique.buffered.shp.sf[1,], fun = "mode")

nrow(cell.df.unique.buffered.shp.sf) # the number of FIDs

tic()

for (p in 1:(4157943 %/% 50000)) {
  
  tree.type.extracted.for.50000 <- exact_extract(tree.type.dat, 
                                               cell.df.unique.buffered.shp.sf[(50000*(p-1) +2):(50000*p+1) ,], 
                                               fun = "mode", max_cells_in_memory = 1e+08)
  
  tree.type.extracted <- c(tree.type.extracted,
                           tree.type.extracted.for.50000)
  
  print(paste0("DONE WITH ", p, "TH BUNCH!"))
  
  toc()
  tic()
  
  if (p %% 10 == 0) { # every 500,000th FID
    save(list = "tree.type.extracted",
         file = paste0("Step2_TreeType_UNTIL", p*50000, ".RData"))
  }
  
}

tree.type.extracted.rest <- exact_extract(tree.type.dat, 
                                        cell.df.unique.buffered.shp.sf[(50000*p +2):nrow(cell.df.unique.buffered.shp.sf),], 
                                        fun = "mode")

tree.type.extracted.all <- c(tree.type.extracted, tree.type.extracted.rest)

tree.type.df <- data.frame(cellID = cell.df.unique$cellid,
                           tree.type = tree.type.extracted.all)

tree.type.df$tree.group <- as.factor(paste0(substr(tree.type.df$tree.type,1,2), "0"))

save(list = "tree.type.df",
     file = "Step2_TreeType_ALL.RData")
