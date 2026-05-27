# Script for plotting fake raster data and calculating CV

## create fake rasters with known shift directions as proof of concept for 
## climate velocity functions

library(terra)
library(tidyterra)
library(ggplot2)
library(dplyr)
library(sf)
library(circular)

source("cv-functions/velocity_functions.R")
source("cv-functions/JL_plotting_checks.R")


# function to make rasters with custom trends and gradients ---------------
make_raster <- function(cold_direction = "north",
                        trend = "both",
                        not_smooth = F,
                        hemisphere = "N",
                        n_lyr = 5,
                        grad_mag = "med", 
                        warming_mag = "med",
                        sd = 0.4){
    
    stopifnot(cold_direction %in% c("north","south","northeast","east"))
    
    # make raster half cooling half warming -----------------------------------
    r <- switch(
        hemisphere,
        "N" = rast(nrows = 18, ncols = 18, xmin = 0, xmax = 18, ymin = 0, ymax = 18),
        "S" = rast(nrows = 18, ncols = 18, xmin = 0, xmax = 18, ymin = -18, ymax = 0)
    )
    
    
    low_temp_val <- 1
    high_temp_val <- switch(grad_mag,
                            "high" = 18,
                            "med" = 9,
                            "low" = 3)
    
    warming_val <- switch(
        warming_mag,
        "low" = .1,
        "med" = .5,
        "high" = 1
    )

    # create default temp gradient --------------------------------------------
    
    if(cold_direction == "south"){
        
        # Generate elevation values with a gradient increasing from south to north
        # Assign higher values to rows closer to the north (ymax)
        gradient_values <- matrix(rep(seq(high_temp_val, low_temp_val, length.out = 18), each = 18), nrow = 18, ncol = 18, byrow = TRUE)
        
    }
    
    if(cold_direction == "north"){
        
        # Generate elevation values with a gradient increasing from south to north
        # Assign higher values to rows closer to the north (ymax)
        gradient_values <- matrix(rep(seq(low_temp_val, high_temp_val, length.out = 18), each = 18), nrow = 18, ncol = 18, byrow = TRUE)
        
    }
    
    if(cold_direction == "northeast"){
        
        # Generate a gradient with the highest elevation at the top-right corner
        # Create horizontal and vertical gradients
        horizontal_gradient <- seq(high_temp_val, low_temp_val, length.out = 18) # Increasing left to right
        vertical_gradient <- seq(low_temp_val, high_temp_val, length.out = 18)  # Increasing bottom to top
        
        # Combine the gradients to form the top-right corner as the highest point
        gradient_values <- outer(vertical_gradient/2, horizontal_gradient/2, "+")
        
        
    }
    
    if(cold_direction == "east"){
        
        # generate a gradient with the coldest values on the right
        gradient_values <- matrix(rep(seq(high_temp_val, low_temp_val, length.out = 18), times = 18), nrow = 18, ncol = 18, byrow = TRUE)
        
    }
    
    # Assign the gradient values to the raster
    values(r) <- gradient_values
    
    # Set raster layer name to "avg_temp"
    names(r) <- "avg_temp"
    
    # Visualize the raster
    # plot(r, main = paste0("Cold Direction: ", cold_direction))
    
    
    # add trend ---------------------------------------------------------------
    # add some randomness
    sd <- sd
    
    if(trend == "warming") {
        
        r_stack <- c(rep(r, times = n_lyr))
        names(r_stack) <- paste0("yr",c(1:n_lyr))
        
        
        for(i in seq_along(1:n_lyr)){
            
            values(r_stack[[i]]) <-  values(r) + rnorm(ncell(r), warming_val*(i-mean(1:n_lyr)), sd)
            print(i)
        }
        
        #return(r_stack)
    }
    
    
    if(trend == "cooling") {
        
        r_stack <- c(rep(r, times = n_lyr))
        names(r_stack) <- paste0("yr",c(1:n_lyr))
        
        
        for(i in seq_along(1:n_lyr)){
            
            values(r_stack[[i]]) <-  values(r) - rnorm(ncell(r), warming_val*(i-mean(1:n_lyr)), sd)
            print(i)
        }
        
    }
    
    
    # if trend is both, make the left half cool and the right half warm
    if(trend == "both") {
        # get matrix of values
        extent <-  switch(
            hemisphere,
            "N" = c(0, 18, 0, 18),
            "S" = c(0, 18, -18, 0)
        )
        
        
        m <- t(matrix(values(r), nrow = 18, ncol = 18) )
        
        list <- vector(mode = "list", length = n_lyr)
        for(i in seq_along(1:n_lyr)){
            
            list[[i]] <- m
            list[[i]][,1:9] <- m[,1:9] + rnorm(ncell(r)/2, warming_val*(i-mean(1:n_lyr)), sd)
            list[[i]][,10:18] <- m[,10:18] - rnorm(ncell(r)/2, warming_val*(i-mean(1:n_lyr)),sd)
            print(i)
        }
        
        list <- purrr::map(
            .x = list,
            .f = ~rast(.x, 
                       extent = extent
            )
        )
        
        r_stack <- rast(list)
        names(r_stack) <- paste0("yr",c(1:n_lyr))
        crs(r_stack) <- crs(r)
        
        # # yr 1
        # m1 <- m
        # m1[,1:9] <- m[,1:9] - rnorm(162, 2, sd)
        # m1[,10:18] <- m[,10:18] + rnorm(162, 2, sd)
        # #  plot(rast(m1))
        # # yr 2
        # m2 <- m
        # m2[,1:9] <- m[,1:9] - rnorm(162, 1, sd)
        # m2[,10:18] <- m[,10:18] + rnorm(162, 1, sd)
        # #  plot(rast(m2))
        # # yr 3
        # m3 <- m
        # m3[,] <- m[,] + rnorm(162*2, 0, sd)
        # #  plot(rast(m3))
        # # yr 4
        # m4 <- m
        # m4[,1:9] <- m[,1:9] + rnorm(162, 1, sd)
        # m4[,10:18] <- m[,10:18] - rnorm(162, 1, sd)
        # #  plot(rast(m4))
        # # yr 4
        # m5 <- m
        # m5[,1:9] <- m[,1:9] + rnorm(162, 2, sd)
        # m5[,10:18] <- m[,10:18] - rnorm(162, 2, sd)
        # #plot(rast(m5))
        # 
        
        
        #  r_stack <- c(rast(m1),
        #               rast(m2),
        #               rast(m3),
        #               rast(m4),
        #               rast(m5))
        #  names(r_stack) <- paste0("yr",1:5)
        #  crs(r_stack) <- crs(r)
        #  
    } # end make half/half trends
    
    # these are currently outfitted so that the mean layer is a smooth gradient,
    # but things will get get weirder if its' not. displace the gradient here
    if(not_smooth){
        
        # displace all the values in the top half of the plot by 2
        add_offset <-  matrix(
            rep(c(2,0), each = ncell(r)/2), nrow = 18, ncol = 18, byrow = T
        )
        
        left_side_cells <- r_stack %>% 
            as.data.frame(cells = T, xy=T) %>% 
            glimpse() %>%
            filter(x <= 9) %>%
            pull(cell)
        
        r_stack[left_side_cells] <- r_stack[left_side_cells] + 10
        plot(r_stack)
        
        
    }
    
    
    return(r_stack)
    
}


make_poly = function(hemisphere = "N"){
    
    # Define points for a simple triangle
    lon <- c(2.5, 3, 14, 13)
    lat <- c(3, 14, 12, 2)
    if(hemisphere == "S")
        lat <- lat*-1
    points_matrix <- cbind(lon, lat)
    SA <- vect(points_matrix, type="polygons", crs="+proj=longlat +datum=WGS84")
    return(SA)
    
}



# make fake raster and SA -------------------------------------------------
r <- make_raster(trend = "cooling",
                 cold_direction = "north",
                 warming_mag = "med")
SA <- make_poly(hemisphere = "N")

res(r)

ggplot() +
    geom_spatraster(data = r) +
    facet_wrap(~lyr) +
    geom_spatvector(dat = SA,
                    color = "darkred", fill = "transparent", linewidth = 1) +
    scale_fill_viridis_c()

# track one isotherm:
points_df <- as.data.frame(r, xy=T) %>%
    tidyr::pivot_longer(cols = contains("yr"),
                        names_to = "lyr", values_to = "temp") 
points_15 <- points_df %>% filter(round(temp) == 5)
ggplot() +
    geom_spatraster(data = r) +
    geom_smooth(data = points_15,
                aes(x=x,y=y), color = "red") +
    geom_point(data = points_15, 
               aes(x=x, y=y), color = "red", alpha = .2) +
    facet_wrap(~lyr) +
    geom_spatvector(dat = SA,
                    color = "darkred", fill = "transparent", linewidth = 1) +
    scale_fill_viridis_c()


# calculate average movement of this isotherm using a loess model
loess_model_yr1 <- loess(data = points_15 %>% filter(lyr == "yr1"),
                             y ~ x)
loess_model_yr5 <- loess(data = points_15 %>% filter(lyr == "yr5"),
                         y ~ x)
# track the velocity of the target isotherm at every x value
# note the mean, because we'll use that to compare to gradient-based CV
x_coords <- unique(points_15$x)
purrr::map(
    .x = x_coords,
    .f = ~{
        
        pos_1 <- predict(loess_model_yr1, list(x = .x))
        pos_2 <- predict(loess_model_yr5, list(x = .x))
        vel <- (pos_2 - pos_1) / (nlyr(r)-1) * 111
        return(vel)
    }
) %>%
    bind_rows() %>%
    summarize(mean = mean(x, na.rm=T))
# note that tracking loess-curve fits like this works well enough for
# comparing to gradient-based climate velocity when the isotherms are 
# only moving directly up or down. If isotherms are diagonal 
# (cold_direction = "northeast", for example), tracking isotherms in 
# every x-cell is going to give inflated values, since the "shortest point" of isotherm
# tracking will be the diagonal, from which we would extract the diagonal
# distance, then the latitudinal component of that distance. 




# 1. calculate trend

# 2.1 find ttrend ----------------------------------------------
ttrend = temp_trend(
    r,
    th = 0.25*nlyr(r)) ## set minimum N obs. to 1/4 time series length
plot(ttrend)
hist(values(ttrend))
mean(values(ttrend))

# 2.2 find spat grad ----------------------------------------------
avg_r <- terra::app(
    r, 
    fun = mean, na.rm = TRUE)
plot(avg_r)

spgrad <- spatial_grad_JL(avg_r)

jl_plot2_spatgrad(avg_r, spgrad, ttrend, SA = SA)


# 3.1 Find directional CV --------------------------------------------------
CV <- gVelocity_JL(grad = spgrad, slope = ttrend, truncate = F)
hist(values(CV$Ang), breaks = 30)

# 3.2 Find latitudinal CV --------------------------------------------------
CV_lat <- CV$Vel * cos(deg_to_rad(CV$Ang))
#CV_lat <- mask(CV_lat, SA)
plot(CV_lat)
mean(values(CV_lat))
mean(values(CV$Vel))
hist(values(CV_lat$Vel))
# visualize
jl_plot3_gvelLat(CV, CV_lat, avg_r, SA, gradfact = .1)

# 3.3 Flip sign for southern hemispherre -----------------------------------
# if polygon is in southern hemisphere, flip the sign of lat velocity
# change sign of gVelLat if in the south hemisphere to reflect a velocity away of the tropics
lat <- terra::init(CV_lat, "y") |> terra::mask(CV_lat)
min_lat <- as.numeric(global(lat,min,na.rm=TRUE))
if(min_lat < 0){
    CV_lat <- ifel(lat < 0, CV_lat * -1, CV_lat)
}
rm(lat, min_lat)


# get all values of raster cells
vals_cv <- values(CV_lat)[,1] #  latitudinal CVs
vals_ttrend <- values(ttrend)[,1] # temperature trends
# make weights vector: coverage fraction of cells * cell size
weights <- values(
    exactextractr::coverage_fraction(CV_lat, st_as_sf(SA))[[1]] * 
        cellSize(CV_lat, mask = T)
)[,1]

# calculate all weighted values:
wt.mean <- Hmisc::wtd.mean(vals_cv, weights)
wt.stdev <- sqrt(Hmisc::wtd.var(vals_cv, weights))
wt.quants <- Hmisc::wtd.quantile(vals_cv, weights, probs = c(.25,.5,.75)) %>% unname()
wt.mean.ttrend <- Hmisc::wtd.mean(vals_ttrend, weights)
wt.stdev.ttrend <- sqrt(Hmisc::wtd.var(vals_ttrend, weights))
wt.mean.directional <- Hmisc::wtd.mean(values(CV$Vel), weights)




out.vals <- data.frame(
    "cv_lat.wtd.mean" = wt.mean,
    "cv_lat.wtd.stdev" = wt.stdev,
    "cv_directional.mean" = wt.mean.directional,
    "cv_directional.meanAng" = wt.mean.directional.angle,
    # "cv_lat.wtd.q25" = wt.quants[1],
    # "cv_lat.wtd.q50" = wt.quants[2],
    #  "cv_lat.wtd.q75" = wt.quants[3],
    "ttrend.wtd.mean" = wt.mean.ttrend,
    "ttrend.wtd.stdev" = wt.stdev.ttrend
)

out.vals$cv_lat.wtd.mean


rm(wt.mean, wt.stdev, wt.quants, wt.mean.ttrend, wt.stdev.ttrend)




