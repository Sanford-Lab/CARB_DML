library(tidyverse);library(terra);library(sf);library(paletteer)

load("Step2.5_FOR_PSCORE_ESTIMATE_FullCARB_NLCDUPDATED_NDVI45.RData")

projs.in.dat <- unique(proj.dats.for.prop.b2k$proj)

quad.and.oct.list <- list()

for (i in 1:length(projs.in.dat)) {
  
  proj.to.assess <- projs.in.dat[i]
  
  proj.df <- proj.dats.for.prop.b1k[proj.dats.for.prop.b1k$proj==proj.to.assess,]
  
  proj.df$quadrants <- NA
  proj.df$octants <- NA
  
  st_centroid(proj.df %>% st_as_sf(.,coords=c("lon","lat"),crs=4269,remove=F)) %>% 
    st_bbox() -> bbox
  
  cent <- unname(c((bbox["xmin"]+bbox["xmax"])/2, (bbox["ymin"]+bbox["ymax"])/2))
  
  proj.octilelines <- terra::vect(matrix(c(bbox["xmin"], bbox["xmax"], bbox["xmin"], bbox["xmax"], 
                                           bbox["xmax"], bbox["xmin"], bbox["xmin"], bbox["xmax"], 
                                           cent[1], cent[1],
                                           bbox["ymin"], bbox["ymax"], bbox["ymax"], bbox["ymin"], 
                                           cent[2], cent[2], bbox["ymin"], bbox["ymin"],
                                           bbox["ymin"], bbox["ymax"]),
                                         ncol = 2),
                                  'lines', crs = 'epsg:4269')
  
  bbox.terra <- terra::vect(terra::ext(terra::vect(proj.df, geom = c("lon", "lat"))), crs = 'epsg:4269')
  
  bbox.spliced <- terra::split(bbox.terra, proj.octilelines)
  
  octants <- terra::extract(bbox.spliced, terra::vect(proj.df, geom = c("lon", "lat")))
  
  rows.for.each.quad <- list()
  
  rows.for.each.quad[["1"]] <- which(proj.df$lon < cent[1] & proj.df$lat < cent[2]) 
  rows.for.each.quad[["2"]] <- which(proj.df$lon < cent[1] & proj.df$lat >= cent[2])
  rows.for.each.quad[["3"]] <- which(proj.df$lon >= cent[1] & proj.df$lat < cent[2])
  rows.for.each.quad[["4"]] <- which(proj.df$lon >= cent[1] & proj.df$lat >= cent[2])
  
  for (j in 1:4){
    if (length(rows.for.each.quad[[as.character(j)]])!=0)  {
      proj.df$quadrants[rows.for.each.quad[[as.character(j)]]] <- j
      
    }
  }
  
  proj.df$octants <- octants[!duplicated(octants[,1]),"id.x"]
  
  proj.df %>% 
    select(proj, treat, FID, quadrants, octants) %>% 
    rename(treat.for.quad = treat)-> proj.df.selected-> proj.df.selected
  
  quad.and.oct.list[[proj.to.assess]] <- proj.df.selected
  
  print(paste0("DONE WITH NO. ", i))
  
  
}

#i = 30, 80, 88, 94 --> errors
projs.errors.quadrants <- projs.in.dat[c(30,80,88,94)]

quad.and.oct <- bind_rows(quad.and.oct.list)

quad.and.oct$proj.and.fid <- paste0(quad.and.oct$proj, "-", quad.and.oct$FID)

for (a in c('proj.dats.for.prop.b1k', 'proj.dats.for.prop.b2k', 
            'proj.dats.for.prop.b2k.qcT')) {
  
  dat <- get(a)
  
  dat$proj.and.fid <- paste0(dat$proj, "-", dat$FID)
  
  dat <- left_join(dat,
                   select(quad.and.oct, c("proj.and.fid", "quadrants")),
                   by = "proj.and.fid")
  
  dat <- dat[,-which(colnames(dat)=="proj.and.fid")]
  
  assign(a, dat)
  
}

save(list = c("quad.and.oct", "projs.errors.quadrants"),
     file = "Step2.5_QUADRANTS.RData")


for (a in c('proj.dats.for.prop.b1k', 'proj.dats.for.prop.b2k', 
            'proj.dats.for.prop.b2k.qcT')) {
  
  dat <- get(a)
  
  dat <- dat[,-which(colnames(dat)=="quadrants.x")]
  
  dat %>% 
    rename(quadrants = quadrants.y) -> dat
  
  assign(a, dat)
  
}

save(list=c('proj.dats.for.prop.b2k', 'proj.dats.for.prop.b2k.qcT', 
            'proj.dats.for.prop.pmethod', 'proj.dats.for.prop.b1k'),
     file = "Step2.5_FOR_PSCORE_ESTIMATE_FullCARB_NLCDUPDATED.RData")

# addition: horizontal/vertical binary & quad2 as centroid from treatment bbox =====

load("Step2.5_FOR_PSCORE_ESTIMATE_FullCARB_NLCDUPDATED_NDVI45.RData")

quad2.and.horiz.vert.list <- list()

for (i in 1:length(projs.in.dat)) {
  
  proj.to.assess <- projs.in.dat[i]
  
  proj.df <- proj.dats.for.prop.b1k[proj.dats.for.prop.b1k$proj==proj.to.assess,]
  
  proj.df$quadrants2 <- NA
  proj.df$binary.horiz <- NA
  proj.df$binary.vert <- NA
  
  st_centroid(proj.df %>% filter(treat==1) %>%  st_as_sf(.,coords=c("lon","lat"),crs=4269,remove=F)) %>% 
    st_bbox() -> bbox.treat
  
  cent <- unname(c((bbox.treat["xmin"]+bbox.treat["xmax"])/2, (bbox.treat["ymin"]+bbox.treat["ymax"])/2))
  
  rows.for.each.quad <- list()
  
  rows.for.each.quad[["1"]] <- which(proj.df$lon < cent[1] & proj.df$lat < cent[2]) 
  rows.for.each.quad[["2"]] <- which(proj.df$lon < cent[1] & proj.df$lat >= cent[2])
  rows.for.each.quad[["3"]] <- which(proj.df$lon >= cent[1] & proj.df$lat < cent[2])
  rows.for.each.quad[["4"]] <- which(proj.df$lon >= cent[1] & proj.df$lat >= cent[2])
  
  for (j in 1:4){
    if (length(rows.for.each.quad[[as.character(j)]])!=0)  {
      proj.df$quadrants2[rows.for.each.quad[[as.character(j)]]] <- j
      
    }
  }
  
  proj.df$binary.vert[proj.df$lat < cent[2]] <- 1
  proj.df$binary.vert[proj.df$lat >= cent[2]] <- 2
  
  proj.df$binary.horiz[proj.df$lon < cent[1]] <- 1
  proj.df$binary.horiz[proj.df$lon >= cent[1]] <- 2
  
  proj.df %>% 
    select(proj, treat, FID, quadrants2, binary.vert, binary.horiz) %>% 
    rename(treat.for.quad = treat)-> proj.df.selected
  #treat.for.quad = same as treat, but just to avoid confusion with original "treat"
  
  quad2.and.horiz.vert.list[[proj.to.assess]] <- proj.df.selected
  
  print(paste0("DONE WITH NO. ", i))
  
  
}

quad2.and.horiz.vert <- bind_rows(quad2.and.horiz.vert.list)

quad2.and.horiz.vert$proj.and.fid <- paste0(quad2.and.horiz.vert$proj, "-", quad2.and.horiz.vert$FID)

for (a in c('proj.dats.for.prop.b1k', 'proj.dats.for.prop.b2k', 
            'proj.dats.for.prop.b2k.qcT')) {
  
  dat <- get(a)
  
  dat$proj.and.fid <- paste0(dat$proj, "-", dat$FID)
  
  dat <- left_join(dat,
                   select(quad2.and.horiz.vert, c("proj.and.fid", "quadrants2", "binary.vert", 'binary.horiz')),
                   by = "proj.and.fid")
  
  dat <- dat[,-which(colnames(dat)=="proj.and.fid")]
  
  assign(a, dat)
  
}

save(list=c('proj.dats.for.prop.b2k', 'proj.dats.for.prop.b2k.qcT', 
            'proj.dats.for.prop.pmethod', 'proj.dats.for.prop.b1k'),
     file = "Step2.5_FOR_PSCORE_ESTIMATE_FullCARB_NLCDUPDATED_NDVI45.RData")

# addition: Quadrants of convex hulls & mean/median of lon and lat =====

cv.hull.list <- list()

for (i in 1:length(projs.in.dat)) {
  
  proj.to.assess <- projs.in.dat[i]
  
  proj.df <- proj.dats.for.prop.b1k[proj.dats.for.prop.b1k$proj==proj.to.assess,]
  
  proj.vect <- terra::vect(proj.df, geom = c("lon", "lat"), crs = "epsg:4269", keepgeom = T)
  
  proj.vect.entire.convhull <- terra::convHull(proj.vect)
  proj.vect.treat.convhull <- terra::convHull(proj.vect[proj.vect$treat == 1,])
  
  proj.df$quadrants.cvENTIRE <- NA
  proj.df$quadrants.cvTREAT <- NA
  proj.df$quadrants.treatBRUTE <- NA
  proj.df$quadrants.treatBRUTEMED <- NA
  
  cent.entire <- terra::geom(terra::centroids(proj.vect.entire.convhull))
  cent.treat <- terra::geom(terra::centroids(proj.vect.treat.convhull))
  cent.treat.brute <- c(mean(proj.df[which(proj.df$treat == 1),]$lon),
                        mean(proj.df[which(proj.df$treat == 1),]$lat))
  cent.treat.brute.med <- c(median(proj.df[which(proj.df$treat == 1),]$lon),
                        median(proj.df[which(proj.df$treat == 1),]$lat))
  
  rows.for.each.quad.entire <- list()
  
  rows.for.each.quad.entire[["1"]] <- which(proj.df$lon < cent.entire[, "x"] & proj.df$lat < cent.entire[, "y"]) 
  rows.for.each.quad.entire[["2"]] <- which(proj.df$lon < cent.entire[, "x"] & proj.df$lat >= cent.entire[, "y"])
  rows.for.each.quad.entire[["3"]] <- which(proj.df$lon >= cent.entire[, "x"] & proj.df$lat < cent.entire[, "y"])
  rows.for.each.quad.entire[["4"]] <- which(proj.df$lon >= cent.entire[, "x"] & proj.df$lat >= cent.entire[, "y"])
  
  for (j in 1:4){
    if (length(rows.for.each.quad.entire[[as.character(j)]])!=0)  {
      proj.df$quadrants.cvENTIRE[rows.for.each.quad.entire[[as.character(j)]]] <- j
      
    }
  }
  
  rows.for.each.quad.treat <- list()
  
  rows.for.each.quad.treat[["1"]] <- which(proj.df$lon < cent.treat[, "x"] & proj.df$lat < cent.treat[, "y"]) 
  rows.for.each.quad.treat[["2"]] <- which(proj.df$lon < cent.treat[, "x"] & proj.df$lat >= cent.treat[, "y"])
  rows.for.each.quad.treat[["3"]] <- which(proj.df$lon >= cent.treat[, "x"] & proj.df$lat < cent.treat[, "y"])
  rows.for.each.quad.treat[["4"]] <- which(proj.df$lon >= cent.treat[, "x"] & proj.df$lat >= cent.treat[, "y"])
  
  for (j in 1:4){
    if (length(rows.for.each.quad.treat[[as.character(j)]])!=0)  {
      proj.df$quadrants.cvTREAT[rows.for.each.quad.treat[[as.character(j)]]] <- j
      
    }
  }
  
  rows.for.each.quad.brute <- list()
  
  rows.for.each.quad.brute[["1"]] <- which(proj.df$lon < cent.treat.brute[1] & proj.df$lat < cent.treat.brute[2]) 
  rows.for.each.quad.brute[["2"]] <- which(proj.df$lon < cent.treat.brute[1] & proj.df$lat >= cent.treat.brute[2])
  rows.for.each.quad.brute[["3"]] <- which(proj.df$lon >= cent.treat.brute[1] & proj.df$lat < cent.treat.brute[2])
  rows.for.each.quad.brute[["4"]] <- which(proj.df$lon >= cent.treat.brute[1] & proj.df$lat >= cent.treat.brute[2])
  
  for (j in 1:4){
    if (length(rows.for.each.quad.brute[[as.character(j)]])!=0)  {
      proj.df$quadrants.treatBRUTE[rows.for.each.quad.brute[[as.character(j)]]] <- j
      
    }
  }
  
  rows.for.each.quad.bruteMED <- list()
  
  rows.for.each.quad.bruteMED[["1"]] <- which(proj.df$lon < cent.treat.brute.med[1] & proj.df$lat < cent.treat.brute.med[2]) 
  rows.for.each.quad.bruteMED[["2"]] <- which(proj.df$lon < cent.treat.brute.med[1] & proj.df$lat >= cent.treat.brute.med[2])
  rows.for.each.quad.bruteMED[["3"]] <- which(proj.df$lon >= cent.treat.brute.med[1] & proj.df$lat < cent.treat.brute.med[2])
  rows.for.each.quad.bruteMED[["4"]] <- which(proj.df$lon >= cent.treat.brute.med[1] & proj.df$lat >= cent.treat.brute.med[2])
  
  for (j in 1:4){
    if (length(rows.for.each.quad.bruteMED[[as.character(j)]])!=0)  {
      proj.df$quadrants.treatBRUTEMED[rows.for.each.quad.bruteMED[[as.character(j)]]] <- j
      
    }
  }
  
  proj.df %>% 
    select(proj, treat, FID, quadrants.cvENTIRE, quadrants.cvTREAT, quadrants.treatBRUTE,
           quadrants.treatBRUTEMED) %>% 
    rename(treat.for.quad = treat)-> proj.df.selected
  #treat.for.quad = same as treat, but just to avoid confusion with original "treat"
  
  cv.hull.list[[proj.to.assess]] <- proj.df.selected
  
  print(paste0("DONE WITH NO. ", i))
  
  
}

cv.hull <- bind_rows(cv.hull.list)

cv.hull$proj.and.fid <- paste0(cv.hull$proj, "-", cv.hull$FID)

for (a in c('proj.dats.for.prop.b1k', 'proj.dats.for.prop.b2k', 
            'proj.dats.for.prop.b2k.qcT')) {
  
  dat <- get(a)
  
  dat$proj.and.fid <- paste0(dat$proj, "-", dat$FID)
  
  dat <- left_join(dat,
                   select(cv.hull, c("proj.and.fid", "quadrants.cvENTIRE", "quadrants.cvTREAT", 'quadrants.treatBRUTE',
                                                  'quadrants.treatBRUTEMED')),
                   by = "proj.and.fid")
  
  dat <- dat[,-which(colnames(dat)=="proj.and.fid")]
  
  assign(a, dat)
  
}

save(list=c('proj.dats.for.prop.b2k', 'proj.dats.for.prop.b2k.qcT', 
            'proj.dats.for.prop.pmethod', 'proj.dats.for.prop.b1k'),
     file = "Step2.5_FOR_PSCORE_ESTIMATE_FullCARB_NLCDUPDATED_NDVI45.RData")

# Quadrants sanity check =====
quad.checker <- function(dat.to.assess, quad.col) {
  projs <- unique(dat.to.assess$proj)
  
  projs.with.no.treat <- c()
  
  for (p in projs) {
    dat.for.this.proj <- filter(dat.to.assess, proj == p)
    
    dat.table <- table(dat.for.this.proj[[quad.col]], dat.for.this.proj[["treat.for.quad"]])
    
    if (any(dat.table[,"1"] == 0)) {
      projs.with.no.treat <- c(projs.with.no.treat, p)
    }
    
  }
  
  return(projs.with.no.treat)
}

length(quad.checker(quad2.and.horiz.vert, "quadrants2")) # 24 projects
length(quad.checker(quad.and.oct, "quadrants")) # 22 projects

length(quad.checker(cv.hull, "quadrants.cvENTIRE")) #22 projects

length(quad.checker(cv.hull, "quadrants.cvTREAT")) # 20 projects

length(quad.checker(cv.hull, "quadrants.treatBRUTE")) # 12 projects

length(quad.checker(cv.hull, "quadrants.treatBRUTEMED")) # 5 projects

cvTREAT.notreat <- quad.checker(cv.hull, "quadrants.cvTREAT")
treatBRUTE.notreat <- quad.checker(cv.hull, "quadrants.treatBRUTE")
treatBRUTEMED.notreat <- quad.checker(cv.hull, "quadrants.treatBRUTEMED")
orig.notreat <- quad.checker(quad.and.oct, "quadrants")

setdiff(orig.notreat, treatBRUTE.notreat) # solved the problem with 12 projects

save(list = c("quad.and.oct", "quad2.and.horiz.vert", "cv.hull", 
              "cvTREAT.notreat", "treatBRUTE.notreat", "treatBRUTEMED.notreat", 'orig.notreat'),
     file = "Step2.5_QUADRANTS.RData")

