library(tidyverse);library(lubridate);library(exactextractr);library(sf)
library(tictoc);library(splines);library(paletteer)

cell.df <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv") # crosswalk file
cell.df.unique <- cell.df[!duplicated(cell.df$cellid),]
cell.df.unique.shp <- terra::vect("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_ALL_NDVI250m.shp")
cell.df.unique.buffered.shp.sf <- st_read("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m.shp")

# Outcome variable: biomass ======

# biomass.dat <- terra::rast("E:/OFFSET_DATA_ALL/Biomass/lt-stem_biomass_nbcd_v0.1_conus_median_2000_2017.tif")

biomass.2018 <-terra::rast("../BiomassData/composite_2018_median.tif")

cell.df.unique.buffered.shp.sf <- st_transform(cell.df.unique.buffered.shp.sf, terra::crs(biomass.2018))

biomass.extracted <- exact_extract(biomass.2018, cell.df.unique.buffered.shp.sf[1,], fun = "mean")

nrow(cell.df.unique) # the number of FIDs

tic()

for (p in 1:(4157943 %/% 100000)) {
  
  biomass.extracted.for.100000 <- exact_extract(biomass.2018, 
                                               cell.df.unique.buffered.shp.sf[(100000*(p-1) +2):(100000*p+1) ,], 
                                               fun = "mean", max_cells_in_memory = 1e+08)
  
  biomass.extracted <- c(biomass.extracted, biomass.extracted.for.100000)
  
  print(paste0("DONE WITH ", p, "TH BUNCH!"))
  
  toc()
  tic()
  
}

biomass.extracted.rest <- exact_extract(biomass.2018, 
                                        cell.df.unique.buffered.shp.sf[(100000*p +2):nrow(cell.df.unique.buffered.shp.sf),], 
                                        fun = "mean")

biomass.2018.all <- c(biomass.extracted, biomass.extracted.rest)

biomass.new.df <- data.frame(cellID = cell.df.unique$cellid)

biomass.new.df[["biomass_2018"]] <- biomass.2018.all

for (yyyy in 2019:2023) {
  biomass.this.yyyy <- terra::rast(paste0("../BiomassData/composite_", yyyy, "_median.tif"))
  
  biomass.extracted.this.yyyy <- exact_extract(biomass.this.yyyy, cell.df.unique.buffered.shp.sf[1,], fun = "mean")
  
  nrow(cell.df.unique) # the number of FIDs
  
  tic()
  
  for (p in 1:(4157943 %/% 300000)) {
    
    biomass.extracted.for.300000 <- exact_extract(biomass.this.yyyy, 
                                                  cell.df.unique.buffered.shp.sf[(300000*(p-1) +2):(300000*p+1) ,], 
                                                  fun = "mean", max_cells_in_memory = 5e+08)
    
    biomass.extracted.this.yyyy <- c(biomass.extracted.this.yyyy, biomass.extracted.for.300000)
    
    print(paste0("DONE WITH ", p, "TH BUNCH!"))
    
    toc()
    tic()
    
  }
  
  biomass.extracted.rest <- exact_extract(biomass.this.yyyy, 
                                          cell.df.unique.buffered.shp.sf[(300000*p +2):nrow(cell.df.unique.buffered.shp.sf),], 
                                          fun = "mean")
  
  biomass.all.this.yyyy <- c(biomass.extracted.this.yyyy, biomass.extracted.rest)
  
  biomass.new.df[[paste0("biomass_,", yyyy)]] <- biomass.all.this.yyyy
  
}

save(list = "biomass.new.df",
     file = paste0("Step2_BIOMASS_NDVI250m_2018to2023.RData"))
