
## Fabio Alexander Castro - Llanos 
## Alliance Bioversity - CIAT 
## June 18th - 2024

## Get the index worldwide

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, ggspatial, ggrepel, rnaturalearthdata, rnaturalearth, climateStability, cptcity, spatialEco, parallelDist, tidyr, tidyverse, gridExtra, rgeos, glue, gtools, readxl, dtw, dtwclust)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------

## Directories (raster data)
root.tasm <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5'
root.prec <- '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/Chirps'

## Get a mask
mask <- dir_ls(root.prec) %>% .[1] %>% rast() %>% ifel(. < 0, 0, .)
mask <- mask * 0 + 1; names(mask) <- 'mask'

## Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50) %>% vect()

## To extract by mask 
mask <- terra::mask(mask, wrld)
mask <- mask * 0 + 1

## Years 
year <- 2015:2023

# To calculate index ------------------------------------------------------

## Number of dry days / Consecutive dry days

### Function to use =---------=
calc.ndd <- function(yr){
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  fls <- dir_ls(root.prec) %>% as.character() 
  fls <- grep(paste0('.', yr, '.'), fls, value = T)
  
  ## To apply by each month
  rst <- map(.x = 1:12, .f = function(m){
    
    ## Filtering the month and reclassified by negative values
    cat('To processing: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), m)
    trr <- rast(grep(paste0('.', yr, '.', mnt, '.'), fls, value = T))
    trr[trr < 0] <- 0
    trr <- terra::mask(trr, wrld)
    
    ## Number of dry days
    ndd <- terra::app(x = trr, fun = function(x){ ndd <- sum(x < 1, na.rm = T); return(ndd)})
    names(ndd) <- glue('ndd_{yr}-{mnt}')
    
    ## Number of consecutive dry days
    calc_cdd <- function(PREC, p_thresh = 1){
      runs <- rle(PREC < p_thresh)
      cons_days <- max(runs$lengths[runs$values==1], na.rm=TRUE)
      return(cons_days)
    }
    cdd <- terra::app(x = trr, fun = function(x){y = calc_cdd(PREC = x, p_thresh = 1); return(y)})
    cdd <- terra::mask(x = cdd, mask = wrld)
    names(cdd) <- glue('cdd_{yr}-{mnt}')
    
    ## To write the rasters 
    terra::writeRaster(x = ndd, filename = glue('./data/tif/index_world/ndd_{yr}-{mnt}.tif'), overwrite = TRUE)
    terra::writeRaster(x = cdd, filename = glue('./data/tif/index_world/cdd_{yr}-{mnt}.tif'), overwrite = TRUE)
    
    ## Return the two final rasters
    rm(trr); gc(reset = TRUE)
    cat('Done :)\n')
    
  })
  
}

### To apply the function =---------=
map(year, calc.ndd)

## Number of days with temperature above 35°C

### Function to use =---------=
calc.ntx <- function(yr){
  
  ## To list the files
  cat('>>> Process: ', yr, '\n')
  dir <- dir_ls(root.tasm) %>% as.character() %>% grep('temperature', ., value = T) %>% grep('maximum', ., value = T)
  fls <- as.character(dir_ls(dir))
  fls <- grep(paste0('_', yr), fls, value = T)
  
  map(.x = 1:12, .f = function(m){
    
    ## Filtering the month
    cat('To process: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    fle <- grep(paste0('_', yr, mnt), fls, value = T)
    
    ## To read the raster and convert the units
    rst <- rast(fle)
    rst <- terra::crop(rst, ext(mask))
    rst <- rst - 273.15
    thr <- 35
    
    ## To calculate the index and write the final 
    ntx <- terra::app(x = rst, fun = function(x){ntxv = sum(x >= thr, na.rm = T); return(ntxv)})
    ntx <- terra::mask(ntx, wrld)
    names(ntx) <- glue('ntx35_{yr}-{mnt}')
    
    ## To write the final raster and finish
    terra::writeRaster(x = ntx, filename = glue('./data/tif/index_world/ntx_{yr}-{mnt}.tif'), overwrite = TRUE)
    rm(rst, ntx); gc(reset = TRUE)
    cat('Done!\n')
    
  })
  
}

### To apply the function =-----------=
map(year, calc.ntx)


















