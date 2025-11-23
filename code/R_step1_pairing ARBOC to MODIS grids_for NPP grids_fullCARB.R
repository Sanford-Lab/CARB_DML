library(tidyverse);library(readxl)

## Difference of this code with the original "step 1" codes ##

# Subsetting ARBOC and forming boundaries =======

arboc.map <- terra::vect("ARBOC Map/OffsetsMapForest.shp")

arboc.map.agg <- terra::aggregate(arboc.map, by = "ARB_Projec")

rm(arboc.map)

arboc.dat <- read_excel("OFFSET_ISSUANCE DATA.xlsx", sheet = "ARB Offset Credit Issuance")

arboc.dat$project.name <-str_split_fixed(arboc.dat$`CARB Issuance ID`, "-",2)[,1] 
arboc.dat <- arboc.dat[arboc.dat$State !="AK",] #leave out Alaska (No climate dataset)

projs.to.assess <- unique(arboc.dat[arboc.dat$`Project Type`=="Forest", ]$project.name)

length(projs.to.assess) #145
length(which(projs.to.assess %in% arboc.map.agg$ARB_Projec)) #137

arboc.dat[arboc.dat$project.name %in% projs.to.assess[!projs.to.assess %in% arboc.map.agg$ARB_Projec],]
# missing from the shapefile: mostly projects that have been reversed

arboc.map.projs <- arboc.map.agg[arboc.map.agg$ARB_Projec %in% projs.to.assess,]
rm(arboc.map.agg)

arboc.map.projs <- terra::project(arboc.map.projs, "epsg:4326")

arboc.map.projs.simpl <- terra::simplifyGeom(arboc.map.projs, tolerance = 0.0025)
#make approx. 250m simplifications
  # required to reduce computational load for buffers
terra::writeVector(arboc.map.projs.simpl, filename = "ARBOC Map/OffsetsMapForestSimplified.shp")

terra::plot(arboc.map.projs[15,])
terra::plot(arboc.map.projs.simpl[15,], add = T, border = 'red')
#seems like a decent simplification 
  # note that the pixels right outside the boundary are left out, 
    # so extreme precision is not required

states <- unique(arboc.dat$State)

arboc.map.projs.b.250 <- terra::buffer(arboc.map.projs.simpl, 250)
arboc.map.projs.b.500 <- terra::buffer(arboc.map.projs.simpl, 500)
arboc.map.projs.b.1k <- terra::buffer(arboc.map.projs.simpl, 1000)
arboc.map.projs.b.2k <- terra::buffer(arboc.map.projs.simpl, 2000)
arboc.map.projs.b.10k <- terra::buffer(arboc.map.projs.simpl, 10000)

arboc.map.projs.b.250.fixed <- terra::vect(terra::geom(arboc.map.projs.b.250)[!is.nan(terra::geom(arboc.map.projs.b.250)[,3]),], 
                                           type = 'polygons', crs = terra::crs(arboc.map.projs.b.250),
                                           atts = terra::as.data.frame(arboc.map.projs.b.250))
# remove invalid points

arboc.map.projs.b.500.fixed <- terra::vect(terra::geom(arboc.map.projs.b.500)[!is.nan(terra::geom(arboc.map.projs.b.500)[,3]),], 
                                           type = 'polygons', crs = terra::crs(arboc.map.projs.b.250),
                                           atts = terra::as.data.frame(arboc.map.projs.b.500))

arboc.map.projs.b.1k.fixed <- terra::vect(terra::geom(arboc.map.projs.b.1k)[!is.nan(terra::geom(arboc.map.projs.b.1k)[,3]),], 
                                          type = 'polygons', crs = terra::crs(arboc.map.projs.b.250),
                                          atts = terra::as.data.frame(arboc.map.projs.b.1k))

arboc.map.projs.b.2k.fixed <- terra::vect(terra::geom(arboc.map.projs.b.2k)[!is.nan(terra::geom(arboc.map.projs.b.2k)[,3]),], 
                                          type = 'polygons', crs = terra::crs(arboc.map.projs.b.250),
                                          atts = terra::as.data.frame(arboc.map.projs.b.2k))

arboc.map.projs.b.10k.fixed <- terra::vect(terra::geom(arboc.map.projs.b.10k)[!is.nan(terra::geom(arboc.map.projs.b.10k)[,3]),], 
                                           type = 'polygons', crs = terra::crs(arboc.map.projs.b.250),
                                           atts = terra::as.data.frame(arboc.map.projs.b.10k))


terra::writeVector(arboc.map.projs.simpl, "Step1_SHAPEFILE_ALL_PROJS_POINTS.shp",
                   overwrite = T)
terra::writeVector(arboc.map.projs.b.250.fixed, 
                   "Step1_SHAPEFILE_ALL_PROJS_B250.shp",
                   overwrite = T)
terra::writeVector(arboc.map.projs.b.500.fixed, 
                   "Step1_SHAPEFILE_ALL_PROJS_B500.shp",
                   overwrite = T)
terra::writeVector(arboc.map.projs.b.1k.fixed, 
                   "Step1_SHAPEFILE_ALL_PROJS_B1k.shp",
                   overwrite = T)
terra::writeVector(arboc.map.projs.b.2k.fixed, 
                   "Step1_SHAPEFILE_ALL_PROJS_B2k.shp",
                   overwrite = T)
terra::writeVector(arboc.map.projs.b.10k.fixed, 
                   "Step1_SHAPEFILE_ALL_PROJS_B10k.shp",
                   overwrite = T)

list.of.polygons.with.invalid.geom <- list()

for (p in c("250", "500", "1k", "2k", "10k")) {
  vect.to.assess <- get(paste0("arboc.map.projs.b.", p))
  
  list.of.polygons.with.invalid.geom[[p]] <- unique(terra::geom(vect.to.assess)[is.nan(terra::geom(vect.to.assess)[,3]),1])
}


projs.to.assess.in.shapefile <- unique(arboc.map.projs$ARB_Projec)


# Bringing in MODIS data and identifying treatment/control points ======

arboc.map.projs.simpl <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_POINTS.shp")
arboc.map.projs.b.250.fixed <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_B250.shp")
arboc.map.projs.b.500.fixed <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_B500.shp")
arboc.map.projs.b.1k.fixed <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_B1k.shp")
arboc.map.projs.b.2k.fixed <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_B2k.shp")
arboc.map.projs.b.10k.fixed <- terra::vect("Step1_SHAPEFILE_ALL_PROJS_B10k.shp")


modis.sample <- terra::rast("MODISSample_NDVI250m_2016Jan1.nc")

arboc.selected.projs.extracted <- terra::extract(modis.sample, 
                                                 terra::project(arboc.map.projs, modis.sample),
                                                 cells = T, xy = T)

arboc.selected.projs.extracted.to <- terra::extract(modis.sample, 
                                                    terra::project(arboc.map.projs, modis.sample),
                                                    cells = T, xy = T, touches =T)
# I use the exact "arboc.map.projs" vector for identifying treatment cells

# I use the buffer vectors (from simplified vector of treatment cells) to identify control treatment
  # 10k: For controls, b2.5 -b2k: to leave out spillover effects

arboc.selected.projs.extracted.b2.5 <- terra::extract(modis.sample, 
                                                      terra::project(arboc.map.projs.b.250.fixed, modis.sample),
                                                      cells = T, xy = T, touches =F)

arboc.selected.projs.extracted.b5 <- terra::extract(modis.sample, 
                                                    terra::project(arboc.map.projs.b.500.fixed, modis.sample),
                                                    cells = T, xy = T, touches =F)

arboc.selected.projs.extracted.b1k <- terra::extract(modis.sample, 
                                                     terra::project(arboc.map.projs.b.1k.fixed, modis.sample),
                                                     cells = T, xy = T, touches =F)

arboc.selected.projs.extracted.b2k <- terra::extract(modis.sample, 
                                                     terra::project(arboc.map.projs.b.2k.fixed, modis.sample),
                                                     cells = T, xy = T, touches =F)

arboc.selected.projs.extracted.b10k <- terra::extract(modis.sample, 
                                                      terra::project(arboc.map.projs.b.10k.fixed, modis.sample),
                                                      cells = T, xy = T, touches =F)

cell.list <- list()

for (i in 1:length(projs.to.assess.in.shapefile)) {
  
  cell.list[[as.character(i)]] <- list()
  
  cell.exact <- arboc.selected.projs.extracted[arboc.selected.projs.extracted$ID == i, "cell"]
  #cells that are exactly inside the boundary (not even touching)
  
  cell.touch <- arboc.selected.projs.extracted.to[arboc.selected.projs.extracted.to$ID == i, "cell"]
  cell.b2.5 <- arboc.selected.projs.extracted.b2.5[arboc.selected.projs.extracted.b2.5$ID == i, "cell"]
  cell.b5 <- arboc.selected.projs.extracted.b5[arboc.selected.projs.extracted.b5$ID == i, "cell"]
  cell.b1k <- arboc.selected.projs.extracted.b1k[arboc.selected.projs.extracted.b1k$ID == i, "cell"]
  cell.b2k <- arboc.selected.projs.extracted.b2k[arboc.selected.projs.extracted.b2k$ID == i, "cell"]
  
  cell.b10k <- arboc.selected.projs.extracted.b10k[arboc.selected.projs.extracted.b10k$ID == i, "cell"]
  #the set of cells that we would collect "control" data from
  
  cell.list[[as.character(i)]][["cell.exact"]] <- cell.exact
  cell.list[[as.character(i)]][["cell.touch"]] <- setdiff(cell.b10k, cell.touch)
  #what we would get if we only leave out the cells that have touched the boundary with the treated forest
  
  cell.list[[as.character(i)]][["cell.b2.5"]] <- setdiff(cell.b10k, cell.b2.5)
  #leave out the 250m buffer zone areas
  cell.list[[as.character(i)]][["cell.b5"]] <- setdiff(cell.b10k, cell.b5)
  cell.list[[as.character(i)]][["cell.b1k"]] <- setdiff(cell.b10k, cell.b1k)
  cell.list[[as.character(i)]][["cell.b2k"]] <- setdiff(cell.b10k, cell.b2k)
  
  cell.list[[as.character(i)]][["cell.prevmethod"]] <- setdiff(cell.b2k, cell.exact)
  #previous method: only sample within 2km distance
  
}

cell.df <- data.frame(proj.id = character(0), # the CAFR project id
                      vect.id = character(0), # the id that the project was given in the vector (1, 2, ...)
                      type = character(0),
                      cellid = character(0))

for (i in 1:length(projs.to.assess.in.shapefile)) {
  
  cell.list.element.to.look.at <- cell.list[[as.character(i)]]
  
  arboc.projs.except.this.proj <- arboc.selected.projs.extracted[arboc.selected.projs.extracted$ID!=i,]
  arboc.projs.except.this.proj.b2.5 <- arboc.selected.projs.extracted.b2.5[arboc.selected.projs.extracted.b2.5$ID!=i,]
  
  # cell.df <- rbind(cell.df,
  #                  data.frame(proj.id = projs.to.assess[i],
  #                             vect.id = i,
  #                             type = "exact.treated.noboundary",  #no buffer for the rest of the projects
  #                             cellid = setdiff(cell.list.element.to.look.at[["cell.exact"]],
  #                                              arboc.projs.except.this.proj$cell)))
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "exact.treated", #buffer for the rest of the projects (default option)
                              cellid = cell.list.element.to.look.at[["cell.exact"]]))
  #cells that are treated, exactly inside 
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.touch", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.touch"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  #cells within the 10km buffer that did not touch with the boundary
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.b2.5", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.b2.5"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  #cells within the 10km buffer outside 250m buffer from the project boundary (deprecated for GPP 500m)
    # we only look at the within-10km buffer that is not included in 250m buffer of other projects
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.b5", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.b5"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  #cells within the 10km buffer outside 500m buffer from the project boundary
  
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.b1k", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.b1k"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  #cells within the 10km buffer outside 1k buffer from the project boundary
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.b2k", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.b2k"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  #cells within the 10km buffer outside 2k buffer from the project boundary
  
  cell.df <- rbind(cell.df,
                   data.frame(proj.id = projs.to.assess.in.shapefile[i],
                              vect.id = i,
                              type = "control.pmethod", 
                              cellid = setdiff(cell.list.element.to.look.at[["cell.prevmethod"]],
                                               arboc.projs.except.this.proj.b2.5$cell)))
  
}

#one last step to obtain the lon/lat coordinates of the cells
cell.df.xy <- terra::vect(terra::xyFromCell(modis.sample, cell.df$cellid))
terra::crs(cell.df.xy)<- terra::crs(modis.sample)
cell.df.lonlat <- terra::project(cell.df.xy, "epsg:4326")

cell.df$lon <- terra::geom(cell.df.lonlat)[,"x"]
cell.df$lat <- terra::geom(cell.df.lonlat)[,"y"]

write.csv(cell.df, "STEP1_ALL_PROJS_LONLATS_NDVI250m.csv") # will be used as a crosswalk
# "We have collected data on this cell: what was it for?"

cell.df.unique <- cell.df[!duplicated(cell.df$cellid),] # data with only cells that are unique
# identified for data retrieval purposes

cell.df.unique.shp <- terra::vect(cell.df.unique[, c("lon", "lat")])
terra::crs(cell.df.unique.shp) <- "epsg:4326" 

terra::writeVector(cell.df.unique.shp, "Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_ALL_NDVI250m.shp")

cell.df.unique.buffered <- terra::buffer(cell.df.unique.shp, sqrt(250*250/3.14))
#make the circle size as similar as possible to 250*250

terra::writeVector(cell.df.unique.buffered, "Step1_SHAPEFILE_FOR_CELLS_TO_STUDY_BUFFERED_ALL_NDVI250m.shp",
                   overwrite = T)

## divvy up the shapefiles, as they are just too big (not fit for GEE "ndvi" dataset)
nrow(cell.df.unique.buffered)

#41 cuts (for NDVI)

for (i in 1:41) {
  if (i != 41) {
    terra::writeVector(cell.df.unique.shp[(100000*(i-1)+1):(100000*i), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/POINTS/POINT_", i, ".shp"))
    
    # terra::writeVector(cell.df.unique.buffered[(100000*(i-1)+1):(100000*i), ],
    #                    paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/BUFFPOLYGON_", i, ".shp"),
    #                    overwrite = T)
  } else if (i==41) {
    terra::writeVector(cell.df.unique.shp[(100000*(i-1)+1):nrow(cell.df.unique.shp), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/POINTS/POINT_", i, ".shp"))
    
    # terra::writeVector(cell.df.unique.buffered[(100000*(i-1)+1):nrow(cell.df.unique.shp), ],
    #                    paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/BUFFPOLYGON_", i, ".shp"),
    #                    overwrite =T)
  }
  
  print(paste0("DONE WITH CELL GROUP ", i))
}

#4 cuts (for fixed)

cell.df.unique.shp$cellID <- cell.df.unique$cellid

for (i in 1:4) {
  if (i != 4) {
    terra::writeVector(cell.df.unique.shp[(1000000*(i-1)+1):(1000000*i), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/POINT_CUT4_", i, ".shp"))
    
  } else if (i==4) {
    terra::writeVector(cell.df.unique.shp[(1000000*(i-1)+1):nrow(cell.df.unique.shp), ],
                       paste0("Step1_SHAPEFILES_DIVIDED_FOR_GEE_TIMEVAR/POINT_CUT4_", i, ".shp"))
    
  }
  
  print(paste0("DONE WITH CELL GROUP ", i))
}


