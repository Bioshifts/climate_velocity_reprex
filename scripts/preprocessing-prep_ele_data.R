# Jake's script for downloading global 1km x 1km rasters 
# from Sunday lab server and cropping them to target study areas. 


path.in <- "/Volumes/general-lab-share/bs-climate-data-processed"
list.files(path.in)


library(BioShiftR)
library(terra)
get_shifts()

get_shifts() %>%
    filter(article_id == "A100") %>%
    add_methods() %>%
    add_trends() %>%
    distinct(midpoint_firstperiod, midpoint_lastperiod, trend_temp_mean)
# ok so we'll need temperature data from 1973 to 1995
start_year <- 1973
end_year <- 1995

for(i in start_year:end_year){
    
    file <- paste0("bios_Ter_",i,".tif")
    
    r <- rast(file.path(path.in, "res-25km","terrestrial",
                       file))
    
    writeRaster(r, 
                file.path("data-raw",
                          "temperature-data",
                          "res-25km",
                          "terrestrial",
                          file),
                overwrite=T)
    
    print(i)
    
}





# get lapse rate ----------------------------------------------------------
# prepare lapse rate raster for elevational CV calculations


# here, we'll upload the lapse rate raster from Chen 2024 (0.5x0.5° grid)
# and resample it to a 1km*1km resolution (without interpolation). 


# libraries ---------------------------------------------------------------
library(terra)
library(dplyr)
library(sf)
library(ggplot2)

# data  -------------------------------------------------------------------

lr <- rast("../BS_climate_velocities/data_local/lapse_rate/cru_ts4.05.2011-2020.malr.tif")
terra::res(lr)
# upload one single temp layer from Sunday server
path <- "/Volumes/general-lab-share/BrunnoData/Land/cruts/bio_proj_1km_temperature"
list.files(path)
r <- rast(file.path(path, list.files(path)[1]))
terra::res(r)
setequal(crs(r), crs(lr))
ext(r)
ext(lr)


# resample lapse rate to small grid ---------------------------------------
# 1. crop lapse rate raster to MAT raster extent
lr_cropped <- crop(lr, r)

# 2. resample without interpolation (nearest neighbor)
lr_resampled <- resample(
    lr_cropped, r, method="near",
    filename="data-raw/lapse-rate/lapse-rate-1km-grid.tif",
    overwrite=TRUE
) # this file gets saved locally, but probably we'll delete it before pushing to github

# 3. checks
compareGeom(lr_resampled, r, stopOnError = FALSE)  # should say TRUE for alignment
all.equal(res(lr_resampled), res(r))
origin(lr_resampled); origin(r)
ext(lr_resampled); ext(r)

plot(lr_resampled)


library(BioShiftR)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyterra)
# for ele, let's use:
# A001_P1, A001_P2, and A005_P1
polys <- get_shifts(type = "ELE") %>%
    add_methods() %>%
    distinct(article_id, poly_id, midpoint_firstperiod, midpoint_lastperiod, duration) %>%
    slice(c(1,2,4)) %>%
    add_polygons() 
bbox <- get_shifts(type = "ELE") %>%
    add_methods() %>%
    distinct(article_id, poly_id, midpoint_firstperiod, midpoint_lastperiod, duration) %>%
    slice(c(1,2,4)) %>%
    add_polygons() %>%
    st_union() %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_as_sf() %>%
    st_buffer(100000) %>%
    st_bbox() %>%
    st_as_sfc() %>%
    st_as_sf()

ggplot() +
    geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf")) +
    geom_sf(data = bbox, color = "red", fill = "transparent") +
    geom_sf(data = polys, color = "purple", fill = "purple", alpha = .2) +
    coord_sf(xlim = c(0, 20),
             ylim = c(35,55))


lr_cropped_2 <- crop(lr_resampled, ext(bbox))


ggplot() +
    geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf")) +
    geom_spatraster(data = lr_cropped_2) +
    scale_fill_viridis_c(na.value = "transparent") +
    geom_sf(data = polys, color = "red", fill = "red", alpha = .7) +
    coord_sf(xlim = c(-5, 20),
             ylim = c(35,55)) 
    


plot(lr_cropped_2)
writeRaster(lr_cropped_2,
            "data-raw/lapse-rate/lapse-rate-1km-cropped.tif",
            overwrite=T)

min(polys$midpoint_firstperiod)
max(polys$midpoint_lastperiod)

# we need: temperature data from 1973 to 2002
temp_path_in <- "/Volumes/general-lab-share/BrunnoData/Land/cruts/bio_proj_1km_temperature"
bbox <- bbox %>% vect()
bbox_coords <- ext(bbox)
for(i in 1975:2002){
    
    
    in_file <- file.path(temp_path_in,
                            paste0("bios_Ter_",i,".tif"))
    
    out_file <-  file.path(
        "data-raw/temperature-data/res-1km",
        paste0("bios_Ter_",i,".tif"))
    
    r_1km <- rast(in_file)
    print(paste0(i," uploaded."))
    
    
    crop(r_1km, ext(bbox),
         filename = file.path(out_file),
         overwrite = T)
    print(paste0(i," saved."))
    
    rm(r_1km)
    gc()
}




