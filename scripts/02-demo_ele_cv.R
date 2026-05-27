# CV demo for elevation


# libraries ---------------------------------------------------------------
library(BioShiftR)
library(ggplot2)
library(terra)
library(sf)
library(tidyterra)
library(dplyr)
library(stringr)


# get cv functions --------------------------------------------------------
source("cv-functions/velocity_functions.R")


# settings ----------------------------------------------------------------
res <- "1km"
# only run 1km for elevation velocities

# get polygons ------------------------------------------------------------
# here, we'll use A001_P1, A001_P2, and A005_P1 to demonstrate, 
# since they are all geographically close (Alps) and over similar durations.
# due to file size, I have cropped 1km temperature data and lapse rate data
# to a bbox surrounding the alps. 

ele_polygons <- 
    # get all latitudinal range shifts
    get_shifts(type = "ELE") %>%
    # add methods (bc we need duration)
    add_methods() %>%
    # find distinct geographies and timeframes
    distinct(article_id, poly_id, eco, midpoint_firstperiod, midpoint_lastperiod) %>%
    # slice the example polygons
    slice(c(1,2,4)) %>%
    # add geometries
    add_polygons() 


# change start/end dates to match data availability -------------------------
# we only have cruts MAT data from 1902- 
# and ORAS SST data from 1959-
# change cv start dates to those years when study duration starts earlier
# (not applicable to this reprex)
ele_polygons <- ele_polygons %>%
    mutate(
        midpoint_firstperiod_cv = case_when(
            eco == "Ter" & midpoint_firstperiod >= 1902 ~ midpoint_firstperiod,
            eco == "Ter" & midpoint_firstperiod < 1902 ~ 1902,
            eco == "Mar" & midpoint_firstperiod >= 1959 ~ midpoint_firstperiod,
            eco == "Mar" & midpoint_firstperiod < 1959 ~ 1959
        )
    )


# create out list
out.list <- vector(mode = "list",
                   length = nrow(ele_polygons))


# upload lapse rate raster ------------------------------------------------
lapse_rate <- rast("data-raw/lapse-rate/lapse-rate-1km-cropped.tif")
temp_data_path <- file.path("data-raw/temperature-data/res-1km")

# loop --------------------------------------------------------------------
for(i in 1:nrow(ele_polygons)){
    
    print("start loop")
    
    poly <- vect(ele_polygons[i,])
    #plot(poly)
    
    start_year <- floor(poly$midpoint_firstperiod)
    end_year <- ceiling(poly$midpoint_lastperiod)
    eco <- poly$eco
    upload_years <- start_year:end_year
    if(eco == "Ter") {
        upload_files <- paste0("bios_Ter_",upload_years,".tif")
    }
    
    # upload and rasterize files
    r <- rast(file.path(temp_data_path,upload_files))
    # add year to layer names
    names(r) <- paste0("mat_", upload_years)
    r <- r * 1
    # plot to see a layer:
    #plot(r[[1]])
    

    # 1. crop to buffered polygon ---------------------------------------------
    # create polygon buffer around the study polygon to incorporate cells from approximately
    # 1 cell width outside the polygon into calculations 
    # (Climate velocity uses 9-cell neighborhoods, so minimizing NAs around focal cells helps)
    buffer_size_m <- as.numeric(str_remove(res,"km"))*1000 # buffer size is in meters >> *1000 converts resol from km to meters
    buffered_poly <- terra::buffer(poly, width = buffer_size_m)
    
    # crop raster
    terra::window(r) <- terra::ext(buffered_poly)
    r <- terra::mask(r, buffered_poly)
    terra::window(r) <- NULL
    
    # convert to real degrees (chelsa uses degrees*10)
    if(eco == "Ter") r <- r/10 
    gc() # clear memory after reducing raster size
    
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = r[[c(1,round(nlyr(r)/2),nlyr(r))]]) +
    #    facet_wrap(~lyr) +
    #    scale_fill_viridis_c(na.value = "transparent") +
    #    geom_spatvector(data = poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1) + 
    #    geom_spatvector(data = buffered_poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1,
    #                    linetype = "dashed")
    
    # crop lapse rate raster to match
    lr <- crop(lapse_rate, buffered_poly)
    
    lr <- mask(lr, buffered_poly)
    
    # un-scale by 100
    lr <- lr/100 # note that this is now in c/km
    
    # 2.1 find ttrend ----------------------------------------------
    ttrend = temp_trend(
        r,
        th = 0.25*nlyr(r)) ## set minimum N obs. to 1/4 time series length
    print("trend calculated")
    
    # 2.2 find avg r ----------------------------------------------
    # note that we don't need to find the spatial gradient here,
    # but we'll calculate average to get a baseline
    avg_r <- terra::app(
        r, 
        fun = mean, na.rm = TRUE)
    #plot(avg_r)
    
    # 3. find CV in the Y-direction (up elevation)
    cv_ele <- ttrend/lr
    names(cv_ele) <- "cv_ele"
    print("velocity calculated")
    
    
    # 4. find average inside polygon ------------------------------------------
    # calculate a few metrics of CV and variance, weighted by cell area and coverage fraction
    
    # get all values of raster cells
    vals_cv <- values(cv_ele)[,1] #  latitudinal CVs
    vals_ttrend <- values(ttrend)[,1] # temperature trends
    vals_baseline <- values(avg_r)[,1]
    
    # make weights vector: coverage fraction of cells * cell size
    weights <- values(
        exactextractr::coverage_fraction(cv_ele, st_as_sf(poly))[[1]] * 
            cellSize(cv_ele, mask = T)
    )[,1]
    
    
    # calculate all weighted values:
    wt.mean <- Hmisc::wtd.mean(vals_cv, weights)
    wt.stdev <- sqrt(Hmisc::wtd.var(vals_cv, weights))
    wt.quants <- Hmisc::wtd.quantile(vals_cv, weights, probs = c(.25,.5,.75)) %>% unname()
    wt.mean.ttrend <- Hmisc::wtd.mean(vals_ttrend, weights)
    wt.stdev.ttrend <- sqrt(Hmisc::wtd.var(vals_ttrend, weights))
    wt.mean.baseline <- Hmisc::wtd.mean(vals_baseline, weights)
    wt.stdev.baseline <- sqrt(Hmisc::wtd.var(vals_baseline, weights))
    print("means calculated")
    
    out.vals <- data.frame(
        "cv_ele.wtd.mean" = wt.mean,
        "cv_ele.wtd.stdev" = wt.stdev,
        "cv_ele.wtd.q25" = wt.quants[1],
        "cv_ele.wtd.q50" = wt.quants[2],
        "cv_ele.wtd.q75" = wt.quants[3],
        "ttrend.wtd.mean" = wt.mean.ttrend,
        "ttrend.wtd.stdev" = wt.stdev.ttrend,
        "baseline.wtd.mean" = wt.mean.baseline,
        "baseline.wtd.stdev" = wt.stdev.baseline
    )
    rm(wt.mean, wt.stdev, wt.quants, wt.mean.ttrend, wt.stdev.ttrend, wt.mean.baseline, wt.stdev.baseline)
    
    # convert cvs to m/year instead of km/year to match elevational range shifts
    out.vals <- out.vals %>%
        mutate(across(contains("cv_"), ~.*1000)) 
    
    
    # 5. save average values --------------------------------------------------
    out.row <- as.data.frame(poly) %>% cbind(out.vals)
    rm(out.vals)
    out.list[[i]] <- out.row
    
    
    # 7. clock ----------------------------------------------------------------
    print(paste0(i," of ",nrow(ele_polygons)))
    
    rm(ttrend, avg_r, r, buffered_poly, cv_ele,
       out.row, eco, end_year, start_year,
       i, upload_years, upload_files,
       weights, vals_cv, vals_ttrend, lr)
    
    gc()
    
}

# bind list
out_bind <- out.list %>%
    bind_rows()



# check saved values for these study areas/durations in bioshifts
get_shifts(type = "ELE") %>%
    add_methods() %>%
    distinct(article_id, poly_id, eco, type, method_id, midpoint_firstperiod, midpoint_lastperiod) %>%
    slice(c(1,2,4)) %>%
    add_trends(res = "1km") %>%
    add_cv(res = "1km") %>%
    add_baselines(res = "1km") %>%
    distinct(article_id, poly_id, midpoint_firstperiod, midpoint_lastperiod,
             across(contains("temp_mean"))) %>%
    glimpse()

# compare to values we just calculated
out_bind %>% 
    select(article_id, poly_id, midpoint_firstperiod,midpoint_lastperiod,
           ttrend.wtd.mean, 
           cv_ele.wtd.mean,
           baseline.wtd.mean) %>%
    glimpse()

rm(list = ls())





