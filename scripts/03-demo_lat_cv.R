# Script to Demo Latitudinal Climate Velocity Calculation

# here, we will use article_polygon A100_P1 as an example, 
# which studies range shifts in Finland over multiple decades.
# This example is pertinent for climate velocities because the study
# area warms in some decades and cools in others.

# In order to limit local file size, I have provided global temperature
# data at 25km x 25km resolution from 1973 to 1995. Thus, we will use 
# the 1974-1984 period (cooling period) and 1984-1994 period (warming period)


# libraries ---------------------------------------------------------------
library(BioShiftR)
library(terra)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyterra)
library(sf)
library(fields)
library(exactextractr)
theme_set(theme_bw())


# source climate velocity functions ---------------------------------------
source("cv-functions/velocity_functions.R")
source("cv-functions/JL_plotting_checks.R")



# make polygon dataframe -------------------------------------------------------
# find distinct geographies and timeframes in BioShifts
# in the "real" project repository, this would list through all study polygons,
# but here, we use only the 2 aforementioned examples. 

lat_polygons <- 
    # get all latitudinal range shifts
    get_shifts(type = "LAT") %>%
    # add methods (bc we need duration)
    add_methods() %>%
    # find distinct geographies and timeframes
    distinct(article_id, poly_id, eco, midpoint_firstperiod, midpoint_lastperiod) %>%
    add_polygons() %>%
    
    # filter to our target polygon/timeframes
    filter(article_id == "A100",
           midpoint_firstperiod %in% c(1974.5, 1984.5),
           midpoint_lastperiod %in% c(1984.5, 1994.5)) 
    

# settings ----------------------------------------------------------------
res <- "25km" # only 25km data provided in this reprex
temp_data_path <- file.path("data-raw/temperature-data")
# set inbound temperature data path (determined by res)


# change start/end dates to match data availability -------------------------
# we only have cruts MAT data from 1902-2016
# and ORAS SST data from 1959-
# change cv start dates to those years when study duration starts earlier
# (not applicable to this reprex)
lat_polygons <- lat_polygons %>%
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
                   length = nrow(lat_polygons))



# loop through polygons and calculate exposure ----------------------------
for(i in 1:nrow(lat_polygons)){
    
    # isolate polygon
    poly <- vect(lat_polygons[i,])
    
    # get focal years
    start_year <- floor(poly$midpoint_firstperiod) # round down for start year
    end_year <- ceiling(poly$midpoint_lastperiod) # round up for end year
    eco <- poly$eco # Mar or Ter; designates which temperature dataset to use
    
    # edit to beginning of temperature series if range shift detection is too long
    if(eco == "Ter" & start_year < 1902) start_year <- 1902
    if(eco == "Mar" & start_year < 1959) start_year <- 1959
    
    # list temperature years to upload
    upload_years <- start_year:end_year
    # make vector of file names to upload
    if(eco == "Ter") {
        upload_files <- paste0("bios_Ter_",upload_years,".tif")
    }
    if(eco == "Mar"){
        upload_files <-paste0("bios_Mar_",upload_years,".tif")
    }
    
    # upload raster files
    # upload rasters 
    path <- file.path(temp_data_path, 
                      paste0("res-",res),
                      ifelse(eco == "Ter","terrestrial","marine"))
    
    # upload and rasterize files
    r <- rast(file.path(path,upload_files))
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
    
    # 2.1 find temperature trend ----------------------------------------------
    # `temp_grad` applies a custom function (`temp_gradFun`) to every raster cell,
    # resulting in a raster layer of the warming trend of each cell (the slope of the linear 
    # model of temperature in each cell over years).
    ttrend = temp_trend(
        r,
        th = 0.25*nlyr(r)) ## set minimum N obs. to 1/4 time series length
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = ttrend) +
    #    scale_fill_gradient2(low = "cornflowerblue", high = "tomato2", na.value = "transparent") +
    #    geom_spatvector(data = poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1) + 
    #    geom_spatvector(data = buffered_poly, 
    #                    fill = "transparent", 
    #                    color = "black",
    #                    linewidth = 1,
    #                    linetype = "dashed")
    
    # 2.2 find spat grad -------------------------------------------------
    # find average temperature in each cell:
    avg_r <- terra::app(r, fun = mean, na.rm = TRUE)
    
    # plot to view:
    #ggplot() +
    #    geom_spatraster(data = avg_r) +
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
    
    # calculate spatial gradient of mean temperature 
    # (the angle and magnitude of average warming)
    spgrad = spatial_grad_JL(avg_r)
    
    # plotting function to view outputs (returns a list of 2 cowplot plots)
    #cv_checks <- jl_plot2_spatgrad(avg_rast = avg_r, 
    #                  spatgrad = spgrad,
    #                  ttrend = ttrend,
    #                  SA = poly,
    #                  gradfact = 20, # multiplier for length of spokes (visual tool)
    #                  spoke_fact = 5 # proportion of spokes to plot (5 = for every 5th cell; visual tool to avoid clutter)
    #                  )
    # view first set:
    # - avg_temp,
    # - WE gradient
    # - NS gradient
    # - Directional spatial gradient of mean temp
    #cv_checks[[1]]
    # view second set:
    # - temperature trend (warming/cooling)
    # - spatial gradient of mean temp (arrows point down)
    # note that from these plots, 
    # we can infer the direction that isotherms should move: 
    # if spatial gradient is cooler in the north, isotherms will
    # be stratified across latitudes, with decreasing values. 
    # thus, a warming trend will make isotherms shift north,
    # and a cooling trend will make isotherms shift south. 
    #cv_checks[[2]]

    
    
    
    # 3.1 find directional climate velocity -----------------------------------
    if(ncell(spgrad) > 100){ truncate = T} else {truncate = F}
    CV <- gVelocity_JL(grad = spgrad, slope = ttrend, truncate = truncate)
    
    # 3.2 Find latitudinal CV --------------------------------------------------
    # extract the northward vector component:
    CV_lat <- CV$Vel * cos(deg_to_rad(CV$Ang))
    # crop back to cells within the focal polygon (removing the buffer cells):
    CV_lat <- mask(CV_lat, poly)
    CV <- mask(CV, poly)
    
    # 3.3 Flip signpoly# 3.3 Flip sign for southern hemispherre -----------------------------------
    # if polygon is in southern hemisphere, flip the sign of lat velocity
    # so that it becomes "poleward" instead of "northward" velocity
    lat <- terra::init(CV_lat, "y") |> terra::mask(CV_lat)
    min_lat <- as.numeric(global(lat,min,na.rm=TRUE))
    if(min_lat < 0){
        CV_lat <- ifel(lat < 0, CV_lat * -1, CV_lat)
    }
    rm(lat, min_lat)
    
    # plot CV and Latitudinal CV
    #cv_checks_2 <- jl_plot3_gvelLat(rast_gvel = CV, 
    #                 rast_gvelLat = CV_lat, 
    #                 avg_rast = avg_r,
    #                 SA = poly,
    #                 spoke_fact = 4,
    #                 gradfact = 1
    #                 )
    #cv_checks_2
    # here, black arrows represent the "directional" cv (at some angle), 
    # and red arrows represent the vector component of the angle facing 0°N
    
    # 4. find average inside polygon ------------------------------------------
    # calculate a few metrics of CV and variance, weighted by cell area and coverage fraction
    
    # get all values of raster cells
    vals_cv <- values(CV_lat)[,1] #  latitudinal CVs
    vals_ttrend <- values(ttrend)[,1] # temperature trends
    vals_baseline <- values(avg_r)[,1]
    # make weights vector: coverage fraction of cells * cell size
    weights <- values(
        exactextractr::coverage_fraction(CV_lat, st_as_sf(poly))[[1]] * 
            cellSize(CV_lat, mask = T)
    )[,1]
    
    # calculate all weighted values:
    wt.mean <- Hmisc::wtd.mean(vals_cv, weights)
    wt.stdev <- sqrt(Hmisc::wtd.var(vals_cv, weights))
    wt.quants <- Hmisc::wtd.quantile(vals_cv, weights, probs = c(.25,.5,.75)) %>% unname()
    wt.mean.ttrend <- Hmisc::wtd.mean(vals_ttrend, weights)
    wt.stdev.ttrend <- sqrt(Hmisc::wtd.var(vals_ttrend, weights))
    wt.mean.baseline <- Hmisc::wtd.mean(vals_baseline, weights)
    wt.stdev.baseline <- sqrt(Hmisc::wtd.var(vals_baseline, weights))
    
    out.vals <- data.frame(
        "cv_lat.wtd.mean" = wt.mean,
        "cv_lat.wtd.stdev" = wt.stdev,
        "cv_lat.wtd.q25" = wt.quants[1],
        "cv_lat.wtd.q50" = wt.quants[2],
        "cv_lat.wtd.q75" = wt.quants[3],
        "ttrend.wtd.mean" = wt.mean.ttrend,
        "ttrend.wtd.stdev" = wt.stdev.ttrend,
        "baseline.wtd.mean" = wt.mean.baseline,
        "baseline.wtd.stdev" = wt.stdev.baseline
    )
    rm(wt.mean, wt.stdev, wt.quants, wt.mean.ttrend, wt.stdev.ttrend, wt.mean.baseline, wt.stdev.baseline)
    
    
    # 5. save average values --------------------------------------------------
    out.row <- as.data.frame(poly) %>% cbind(out.vals)
    rm(out.vals)
    out.list[[i]] <- out.row
    
    # clock for loop
    print(i)
    
}


# bind list
out_bind <- out.list %>%
    bind_rows()

# check saved values for these study areas/durations in bioshifts
get_shifts() %>%
    filter(article_id == "A100") %>%
    add_methods() %>%
    filter((midpoint_firstperiod == 1974.5 & midpoint_lastperiod == 1984.5)|
               (midpoint_firstperiod == 1984.5 & midpoint_lastperiod == 1994.5)) %>%
    add_trends(res = "25km") %>%
    add_cv(res = "25km") %>%
    add_baselines(res = "25km") %>%
    distinct(article_id, poly_id, midpoint_firstperiod, midpoint_lastperiod,
             across(contains("temp_mean"))) %>%
    glimpse()

# compare to values we just calculated
out_bind %>% 
    select(article_id, poly_id, midpoint_firstperiod,midpoint_lastperiod,
           ttrend.wtd.mean, 
           cv_lat.wtd.mean,
           baseline.wtd.mean) %>%
    glimpse()
