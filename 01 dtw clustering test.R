
## Fabio Alexander Castro - Llanos 
## Alliance Bioversity - CIAT 
## May 9th - 2024

## Get the DTW clustering

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, ggspatial, ggrepel, climateStability, cptcity, spatialEco, parallelDist, tidyr, tidyverse, gridExtra, rgeos, glue, gtools, readxl, dtw, dtwclust)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions ---------------------------------------------------------------
inv <- function(x){((x - max(x, na.rm = TRUE)) * -1) + min(x, na.rm = TRUE)}

# Load data ---------------------------------------------------------------

## Tabular data
tble <- read_csv('../V1/r/tbl/values/climate-vars_monthly.csv', show_col_types = FALSE)
tble <- mutate(tble, iso = str_sub(id, nchar(id) - 2, nchar(id)))
isos <- unique(tble$iso)
cntr <- 'PER'
dfrm <- filter(tble, iso == cntr)
crds <- dfrm %>% distinct(id, lon, lat) 

## Raster data 
root.tasm <- '//catalogue/WFP_ClimateRiskPr1/1.Data/AgERA5'
root.prec <- '//CATALOGUE/WFP_ClimateRiskPr1/1.Data/Chirps'

## Vector data 
lim0 <- geodata::gadm(country = cntr, level = 0, path = './tmpr')
lim1 <- geodata::gadm(country = cntr, level = 1, path = './tmpr')
lim1 <- st_as_sf(lim1)

# Extract by mask  --------------------------------------------------------

## Temperature ---------------------------------------- 
dirs.tasm <- as.character(grep('_temperature-', (dir_ls(root.tasm, type = 'directory')), value = T))
map(.x = dirs.tasm, .f = function(d){
  
  ## To list the files -----------
  cat('To process: ', basename(d), '\t')
  fls <- as.character(dir_ls(d))
  vrb <- basename(d)
  yrs <- 1993:2023
  fls <- grep(paste0(yrs, collapse = '|'), fls, value = T)  
  fls %>% basename() %>% str_split(., pattern = '_') %>% map_chr(., 4) %>% str_sub(., 1, 4) %>% unique() %>% sort()
  
  ## To extract by mask -----------
  rstr.mnth <- map(.x = yrs, .f = function(y){
    
    ### To extract
    cat('>>> Processing: ', y, '\n')
    rst <- rast(grep(paste0('_', y), fls, value = T))
    rst <- terra::crop(rst, lim0)
    rst <- terra::mask(rst, lim0)
    
    ### To average by each year
    rst.mnt <- map(.x = 1:12, .f = function(m){
      cat('Month: ', m, '\n')
      rs <- mean(rst[[grep(paste0(y, '-', ifelse(m < 10, paste0('0', m), as.character(m)), '-'), time(rst))]])
      names(rs) <- glue('{y}-{m}')
      return(rs)
    }) %>% reduce(., c)
    
    ## Finish and return
    return(rst.mnt)
    
  })
  
  ## To write the raster
  rstr.mnth <- reduce(rstr.mnth, c)
  terra::writeRaster(x = rstr.mnth, filename = glue('./data/tif/peru/{vrb}_{yrs[1]}-{yrs[length(yrs)]}.tif'), overwrite = TRUE)
  
})

## Precipitation --------------------------------------
fles <- as.character(dir_ls(root.prec))
fles <- grep(paste0(1993:2023, collapse = '|'), fles, value = T)
prec <- map(.x = 1993:2023, .f = function(y){
  cat('To process: ', y, '\n')
  fls <- grep(paste0('.', y, '.'), fles, value = T)
  ppt <- map(.x = 1:12, .f = function(m){
    cat('Month: ', month.abb[m], '\n')
    mnt <- ifelse(m < 10, paste0('0', m), as.character(m))
    rst <- rast(grep(paste0(y, '.', mnt), fls, value = T))
    rst <- terra::crop(rst, lim0)
    rst <- terra::mask(rst, lim0)
    sma <- sum(rst)
    names(sma) <- glue('prec_{y}-{m}')
    return(sma)    
  }) %>% 
    reduce(., c)
  return(ppt)
})
prec <- reduce(prec, c)
terra::writeRaster(x = prec, filename = glue('./data/tif/peru/prec_1993-2023.tif'), overwrite = TRUE)

# Read the results --------------------------------------------------------
prec <- terra::rast('./data/tif/peru/prec_1993-2023.tif')
tmin <- terra::rast('./data/tif/peru/2m_temperature-24_hour_minimum_1993-2023.tif')
tmax <- terra::rast('./data/tif/peru/2m_temperature-24_hour_maximum_1993-2023.tif')
map_dbl(list(prec, tmin, tmax), nlyr)

# To calculate the average by year  ---------------------------------------
prec.sum <- map(.x = 1993:2023, .f = function(i){r <- sum(prec[[grep(paste0('_', i, '-'), names(prec))]]); r[r < 0] <- 0; names(r) <- glue('prec_{i}'); return(r)}) %>% reduce(., c)
tmin.avg <- map(.x = 1993:2023, .f = function(i){r <- mean(tmin[[grep(paste0('', i, '-'), names(tmin))]]); names(r) <- glue('tmin_{i}'); return(r)}) %>% reduce(., c)
tmax.avg <- map(.x = 1993:2023, .f = function(i){r <- mean(tmax[[grep(paste0('', i, '-'), names(tmax))]]); names(r) <- glue('tmax_{i}'); return(r)}) %>% reduce(., c)
prec.sum <- terra::resample(prec.sum, tmin.avg, method = 'bilinear')
stck.avg <- c(prec.sum, tmin.avg, tmax.avg)

# To extract the values for the points ------------------------------------
extr.vles <- function(ide){
  crd <- filter(crds, id == ide)
  rsl <- as_tibble(cbind(crd, terra::extract(stck.avg, crd[,c('lon', 'lat')])))
  return(rsl)
}
vles <- map_dfr(.x = unique(crds$id), .f = extr.vles)
tble <- drop_na(terra::as.data.frame(stck.avg, cell = T, xy = T))
indx <- dplyr::select(tble, -x, -y)

# Check the cell IDs ------------------------------------------------------
ref <- rast(tble[,2:4], type = 'xyz') * 0 + 1
n  <- terra::cellFromXY(ref, as.matrix(crds[,2:3]))
cellID <- tble$cell

ppt <- cbind(cellID = cellID, dplyr::select(tble, contains('prec')))
tmx <- cbind(cellID = cellID, dplyr::select(tble, contains('tmax')))
tmn <- cbind(cellID = cellID, dplyr::select(tble, contains('tmin')))

# Function ----------------------------------------------------------------
DTW_mul <- function(pp, tx, tn, rf){
  
  # pp <- ppt
  # tx <- tmx
  # tn <- tmn
  # rf <- n[1]
  
  dist_mul <- lapply(1:nrow(tn), function(j){
    
    index_ref_1 <- pp[which(pp$cellID == rf), -1]
    colnames(index_ref_1) <- as.factor(1:31)
    
    index_ref_2 <- tx[which(tx$cellID == rf), -1]
    colnames(index_ref_2) <- as.factor(1:31)
    
    index_ref_3 <- tn[which(tn$cellID == rf),-1]
    colnames(index_ref_3) <- as.factor(1:31)
    
    index_ref <- rbind(index_ref_1,index_ref_2,index_ref_3)
    
    #### VS 
    
    index_1 <- pp[j,-1]
    colnames(index_1) <- as.factor(1:31)
    
    index_2 <- tx[j,-1]
    colnames(index_2) <- as.factor(1:31)
    
    index_3 <- tn[j,-1]
    colnames(index_3) <- as.factor(1:31)
    
    indexes  <- rbind(index_1,index_2,index_3)
    
    df <- rbind(index_ref,indexes)
    
    sample.matrix <- as.matrix(df)
    
    sample.matrix1 <- sample.matrix[1:3,]
    sample.matrix2 <- sample.matrix[4:6,]
    
    sample.List <- list(sample.matrix1, sample.matrix2)
    dist <- parDist(x = sample.List, method = "dtw")
    result <- data.frame(cellID= tmx$cellID[j], tdist= dist[1])
    return(result)
    
  })
  
  dtw <- do.call(rbind, dist_mul)
  dtw$tdist[which(dtw$tdist == Inf)] <- mean(dtw$tdist)
  return(dtw)
  
}

bm <- lapply(1:2, function(i){
  cat(paste0("Procesando ref ::: ",i, "\n"))
  DTW_mul(pp = ppt , tx = tmx, tn = tmn, rf= n[i])
})

# Results to raster  ------------------------------------------------------
r1 <- terra::rast(cbind(tble[,2:3], rsl_1 = bm[[1]][,2]))
plot(r1)
points(crds[1,]$lon, crds[1,]$lat, pch = 16, col = 'red')

r2 <- terra::rast(cbind(tble[,2:3], rsl_2 = bm[[2]][,2]))
plot(r2)
points(crds[2,]$lon, crds[2,]$lat, pch = 16, col = 'red')

# To invert the values ----------------------------------------------------
rr <- c(rescale0to1(r1), rescale0to1(r2))
rr <- terra::app(rr, function(x) 1 - x)
rr <- 1 - rr
names(rr) <- crds$id

# To make the map ---------------------------------------------------------
rr.tb <- terra::as.data.frame(rr, xy = T) %>% as_tibble() %>% gather(var, value, -c(x, y))

grsl <- ggplot() + 
  geom_tile(data = rr.tb, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var) +
  scale_fill_gradientn(colors = cpt(pal = 'imagej_gyr_centre', n = 10, rev = TRUE)) +
  geom_point(data = crds, aes(x = lon, y = lat), col = 'red', pch = 16) +
  geom_text_repel(data = crds, aes(x = lon, y = lat, label = id), bg.color = 'white', bg.r = 0.25) +
  labs(fill = 'DTW mult normalized') +
  geom_sf(data = lim1, fill = NA, col = 'grey30') +
  coord_sf() + 
  theme_void() + 
  theme(legend.position = 'bottom', 
        legend.key.width = unit(3, 'line'), 
        strip.text = element_text(face = 'bold')) +
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = grsl, filename = './jpg/map_1.jpg', units = 'in', width = 9, height = 7, dpi = 300)
