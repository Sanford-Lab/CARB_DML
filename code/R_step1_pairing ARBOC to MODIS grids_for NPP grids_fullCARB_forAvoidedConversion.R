library(tidyverse); library(readxl)

cell.df <- read_csv("STEP1_ALL_PROJS_LONLATS_NDVI250m.csv") # crosswalk file

arboc.dat.class <- read_excel("OFFSET_ISSUANCE DATA_ProjectType.xlsx")
refor.and.ac.projs <- arboc.dat.class[arboc.dat.class$type!="IFM",]$proj.name

refor.and.ac.projs.in.dat <- refor.and.ac.projs[refor.and.ac.projs %in% cell.df$proj.id]
# reforestation projects

# see that the projects that are not in the dataset have never been given any credit

cells.in.refor.and.ac.projs <- cell.df[cell.df$proj.id %in% refor.and.ac.projs.in.dat,] 

rm(cell.df)

cells.in.refor.and.ac.projs.unique <- cells.in.refor.and.ac.projs[!duplicated(cells.in.refor.and.ac.projs$cellid),] # data with only cells that are unique
# identified for data retrieval purposes

cell.df.unique.ac.shp <- terra::vect(cells.in.refor.and.ac.projs.unique[, c("cellid", "lon", "lat")], geom = c("lon", "lat"))
terra::crs(cell.df.unique.ac.shp) <- "epsg:4326" 

terra::writeVector(cell.df.unique.ac.shp, "Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_ALL_NDVI250m_AC.shp", overwrite = T)

cell.df.unique.ac.buffered <- terra::simplifyGeom(terra::buffer(cell.df.unique.ac.shp, sqrt(250*250/3.14)),
                                                    tolerance = 0.0005)

terra::writeVector(cell.df.unique.ac.buffered,
                   paste0("Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m_AC.shp"),
                   overwrite = T)

nrow(cell.df.unique.ac.buffered)

for (i in 1:8) {
  if (i != 8) {
    terra::writeVector(cell.df.unique.ac.buffered[(50000*(i-1)+1):(50000*i), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/ForAC_POLYGONS/BUFFAC_", i, ".shp"))
    
    
  } else if (i==8) {
    terra::writeVector(cell.df.unique.ac.buffered[(50000*(i-1)+1):nrow(cell.df.unique.ac.buffered), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/ForAC_POLYGONS/BUFFAC_", i, ".shp"))
    
    
  }
  
  print(paste0("DONE WITH CELL GROUP ", i))
}
 
save(list=c("cells.in.refor.and.ac.projs.unique", "cells.in.refor.and.ac.projs"),
     file = "STEP1_AC_LONLATS_NDVI250m.RData")
