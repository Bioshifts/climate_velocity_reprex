# JL Plotting Functions

# functions to plot outputs of climate velocity steps


# plot normal resoliution and downscaled resolution side by side
# (irrelevant here since we only have .25 res)
JL_plot1_both_res <- function(original, downscaled, SA){
    
    # test both
    globalmin <- min(c(values(original[[1]], na.rm=T), values(downscaled[[1]], na.rm=T)))
    globalmax <- max(c(values(original[[1]], na.rm=T), values(downscaled[[1]], na.rm=T)))
    plotlims <- ext(terra::buffer(SA, width = 150000))
    
    test_p1 <- ggplot() +
        geom_spatraster(data = original[[1]]) +
        scale_fill_viridis_c(na.value = "transparent",
                             limits = c(globalmin, globalmax)) +
        geom_spatvector(data = SA,
                        color = "darkred", 
                        fill = "transparent",
                        linewidth = 1) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F)
    
    test_p2 <- ggplot() +
        geom_spatraster(data = downscaled[[1]]) +
        scale_fill_viridis_c(na.value = "transparent",
                             limits = c(globalmin, globalmax))+
        geom_spatvector(data = SA,
                        color = "darkred", 
                        fill = "transparent",
                        linewidth = 1) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F)
    
    p <- cowplot::plot_grid(test_p1, test_p2)
    return(p)
}



# Plot spatial gradient and ttrend
jl_plot2_spatgrad <- function(avg_rast, spatgrad, ttrend, SA, gradfact = 6, prefix = NULL, spoke_fact = 1 ){
    
    if(!is.null(prefix)) prefix <- paste0(prefix,": ")
    
    plotlims <- ext(terra::buffer(SA, width = 100000))
    
    linewidth <- ifelse(ncell(avg_rast) > 50000, .03, .3)
    
    
    test_p1 <- ggplot() +
        geom_spatraster(data = avg_rast) +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        scale_fill_viridis_c(na.value = "transparent") +
        geom_spatvector(data = SA,
                        fill = "transparent",
                        color = "darkred",
                        linewidth = 1) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F) +
        labs(title =paste0(prefix, "Avg Temp"),
             fill = "Mean") +
        theme(legend.position = "bottom")
    
    we_ang <- spatgrad[[1]] %>%
        as.data.frame(xy=T) %>%
        # geom_spoke uses "math" radians, e.g., 
        # radians calculated from "math" degrees, not "compass" degrees
        # (e.g., 0° is right and goes counterclockwise)
        # so convert compass angles to those here. 
        mutate(
            # positive values for "WEgrad" mean that east is warmer, 
            # so the default angle is 90 (compass degrees east)
            angle = 90,
            # convert to math degrees
            mathDegrees = (450 - angle) %% 360,
            # convert to math radians
            mathRadians =  DescTools::DegToRad(mathDegrees)) 
    
    test_p2_we <- ggplot() +
        
        geom_spatraster(data = spatgrad[[1]]) +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        
        scale_fill_gradient2(high = "purple",
                             low = "turquoise",
                             midpoint = 0,
                             na.value = "transparent") +
        geom_spatvector(data = SA,
                        fill = "transparent",
                        color = "darkred",
                        linewidth = 1) +
        geom_spoke(data = we_ang %>%
                       filter(row_number() %% spoke_fact == 0),
                   aes(x=x,
                       y=y,
                       angle = mathRadians,
                       radius = WE *gradfact),
                   arrow = arrow(length = unit(2,"mm"))) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F) +
        labs(title = paste0(prefix, "WE gradient"),
             fill = "Grad") +
        theme(legend.position = "bottom",
              axis.title = element_blank())
    
    
    
    # again, convert NSgrad to "math" degrees then "math" radians
    ns_ang <- spatgrad[[2]] %>%
        as.data.frame(xy=T) %>%
        mutate(
            # positive values in NSgrad mean north is warmer,
            # so default angle is 0 (compass bearing north)
            angle = 0,
            # convert to math degrees
            mathDegrees = (450 - angle) %% 360,
            # convert to math radians
            mathRadians =  DescTools::DegToRad(mathDegrees)) 
    
    test_p3_ns <- ggplot() +
        geom_spatraster(data = spatgrad[[2]]) +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        
        scale_fill_gradient2(high = "tomato2",
                             low = "cornflowerblue",
                             midpoint = 0,
                             na.value = "transparent") +
        geom_spatvector(data = SA,
                        fill = "transparent",
                        color = "darkred",
                        linewidth = 1) +
        geom_spoke(data = ns_ang %>%
                       filter(row_number() %% spoke_fact == 0),
                   aes(x=x,
                       y=y,
                       angle = mathRadians,
                       radius = NS *gradfact),
                   arrow = arrow(length = unit(2,"mm"))) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F) +
        labs(title = paste0(prefix, "NS gradient"),
             fill = "Grad") +
        theme(legend.position = "bottom",
              axis.title = element_blank())
    
    # plot angle
    angle_df <- spatgrad[[3:4]] %>%
        as.data.frame(xy=T)  %>%
        mutate(
            # convert compass degrees to math degrees
            mathDegrees = (450 - angle) %% 360 ,
            # convert math degrees to radians
            mathRadians = DescTools::DegToRad(mathDegrees))
    
    test_p4_tempspat <- angle_df %>%
        ggplot() +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        geom_tile(data = avg_rast %>% as.data.frame(xy=T),
                  aes(x=x,y=y, fill = mean),
                  alpha = .5,
                  color = "transparent") +
        scale_fill_viridis_c(na.value = "transparent") +
        geom_spoke(data = . %>%
                       filter(row_number() %% spoke_fact == 0), 
                   aes(x=x,
                       y=y,
                       radius = Grad*gradfact,
                       angle = mathRadians),
                   arrow = arrow(length = unit(2,"mm"))) +
        coord_equal() +
        geom_spatvector(data = SA,
                        fill = "transparent",
                        color = "darkred",
                        linewidth = .5) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F) +
        labs(title = paste0(prefix,"Temp SpatGrad"),
             subtitle = "Arrows Point towards Warm",
             fill = "Mean") +
        theme(legend.position = "bottom",
              axis.title = element_blank())
    
    
    
    p <- cowplot::plot_grid(test_p1, test_p2_we, test_p3_ns, test_p4_tempspat, nrow = 1)
    
    
    test_p5_ttrend <- 
        ggplot() +
        geom_spatraster(data = ttrend) +
        scale_fill_gradient2(low = "cornflowerblue",
                             high = "tomato2",
                             na.value = "transparent") +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        geom_spatvector(data = SA,
                        fill = "transparent",
                        color = "darkred",
                        linewidth = .5) +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F)  +
        labs(title = paste0(prefix,"Warming Trend"),
             subtitle = "Blue = Cooling; Red = Warming") +
        theme(legend.position = "bottom")
    
    p2 <- cowplot::plot_grid(test_p5_ttrend, test_p4_tempspat)
    
    return(list(p, p2))
    
}


# plot directional climate velocity and poleward climate velocity vectors:
jl_plot3_gvelLat <- function(rast_gvel, rast_gvelLat, avg_rast,SA, spoke_fact = 1, gradfact = 6, prefix = NULL ){
    
    if(!is.null(prefix)) prefix <- paste0(prefix,": ")
    
    plotlims <- ext(terra::buffer(SA, width = 100000))
    
    # set linewidth to be small / nonexistent when many cells
    linewidth <- ifelse(ncell(rast_gvel) > 50000, .03, .3)
    
    gVel_df <- rast_gvel %>%
        as.data.frame(xy=T) %>%
        mutate(
            # convert compass degrees to math degrees
            mathDegrees = (450 - Ang) %% 360 ,
            # convert math degrees to radians
            mathRadians = DescTools::DegToRad(mathDegrees)) 
    
    
    gVelLat_down_df <- rast_gvelLat %>%
        as.data.frame(xy=T) %>%
        mutate(
            Ang = 0,
            # convert compass degrees to math degrees
            mathDegrees = (450 - Ang) %% 360 ,
            # convert math degrees to radians
            mathRadians = DescTools::DegToRad(mathDegrees)) 
    
    
    test_p8 <- 
        ggplot() +
        geom_spatvector(data = SA,
                        fill = "grey80",
                        color = "transparent",
                        linewidth = .5) +
        geom_spatvector(data = as.polygons(avg_rast, aggregate = F),
                        color = "black",
                        fill = "transparent",
                        linewidth = linewidth) +
        
        geom_spoke(data = gVelLat_down_df %>%
                       filter(row_number() %% spoke_fact == 0),
                   aes(x=x,
                       y=y,
                       radius = Vel/20 * gradfact,
                       angle = mathRadians),
                   arrow = arrow(length = unit(2,"mm")),
                   color = "red") +
        geom_spoke(data = gVel_df %>% 
                       filter(row_number() %% spoke_fact == 0),
                   aes(x=x,
                       y=y,
                       radius = Vel/20 * gradfact,
                       angle = mathRadians),
                   arrow = arrow(length = unit(2,"mm")),
                   color = "black") +
        coord_equal() +
        coord_sf(xlim = c(plotlims$xmin, plotlims$xmax),
                 ylim = c(plotlims$ymin, plotlims$ymax),
                 expand = F) +
        labs(title = paste0(prefix,"Climate Velocity"),
             subtitle = paste0("Black: Full (mean = ",round(mean(gVel_df$Vel, na.rm=T),1),") ;\nRed: Latitudinal (mean = ",round(mean(gVelLat_down_df$Vel, na.rm=T)),")"),
             fill = "Mean") +
        theme(legend.position = "bottom",
              axis.title = element_blank())
    
    
    #egg::ggarrange(test_p5_ttrend, test_p4_tempspat, test_p8, nrow = 1)
    return(test_p8)
    
    
}


