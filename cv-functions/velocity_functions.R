# Authors: Brunno Oliveira, edited by Jake Lawlor
# Last update : Jan 2026
#--------
# This workflow is based on the functions from the climetrics package (https://github.com/shirintaheri/climetrics)
# and the gVelocity workflow (functions `tempTrend`, `spatGrad`, `gVoCC`) from the (depreciated) VoCC package (https://github.com/JorGarMol/VoCC)




# temp trend internal helpers ---------------------------------------------
# extract slope coefficient of linear model
temp_gradFun <- function(x, th) {
    # keep only non-NA cells 
    # (inherited; there should never really be any random NAs anyway)
    x <- x[!is.na(x)]
    # regress temp values over length
    # and extract coefficient 
    if (length(x) > th) {
        s <- lm(x~c(1:length(x)))
        s$coefficients[2] # extract slope
    } else NA
}

# temp trend function -----------------------------------------------------
# function to apply temp_gradFun to every raster cell
temp_trend <- function(x, th, ncores=NULL, file_name="",overwrite=TRUE) {
    if(is.null(ncores)){
        tmp <- terra::app(x,temp_gradFun,th=th,filename=file_name,overwrite=overwrite)
    } else {
        tmp <- terra::app(x,temp_gradFun,th=th,cores=ncores,filename=file_name,overwrite=overwrite)
    }
    names(tmp) <- "Trend" # name the layer
    return(tmp)
}



# spatial gradient internal helpers ---------------------------------------

#----
# Utils
.is_package_installed <- function(n) {
    names(n) <- n
    sapply(n, function(x) length(unlist(lapply(.libPaths(), function(lib) find.package(x, lib, quiet=TRUE, verbose=FALSE)))) > 0)
}
#----
# get projection type
.getProj <- function(x) {
    if (inherits(x,'Raster')) {
        if (!is.na(projection(x))) strsplit(strsplit(projection(x),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    } else {
        if (!is.na(crs(x))) strsplit(strsplit(crs(x,proj=TRUE),'\\+proj=')[[1]][2],' ')[[1]][1]
        else {
            if (all(extent(x)[1:2] >= -180 & extent(x)[1:2] <= 180 & extent(x)[3:4] >= -90 & extent(x)[3:4] <= 90)) 'longlat'
            else 'projected'
        }
    }
}
#-----
# weighted mean of 6 pairwise comparisons in a 9-cell neighborhood
# where pairs containing the focal cell (middle) are counted twice
.mnwm <- function(d1, d2, d3, d4, d5, d6){
    X <- sum(c(d1, d2*2, d3, d4, d5*2, d6), na.rm = T)
    w <- sum(c(1,2,1,1,2,1) * is.finite(c(d1, d2, d3, d4, d5, d6)))
    return(X/w)
}
#-----
# find angle of warmer temperature
# takes dx (change in temperature as you move right/east)
# and dy (change in temperature as you move up/north)
# and returns the compass angle of temperature increase
.ang <- function(dx, dy){
    ifelse(dy < 0, 180 + rad_to_deg(atan(dx/dy)),
           ifelse(dx < 0, 360 + rad_to_deg(atan(dx /dy )), rad_to_deg(atan(dx/dy))))
}
# JL NOTE: I think we could simplify this to:
# (atan2(dx, dy) * 180/pi) %% 360
# currently, the ifelse is "quadrant correcting", but whatever. 

# sanity check:
#.ang(dx = 1, dy = 0) # temperature warms to the right (90)
#.ang(dx = -1, dy = -1) # temperature warms to SW (225)
#.ang (dx = 1, dy = 1) # temperature warms NE (45)


#---
# change from degrees to radians
deg_to_rad <- function (degree) {
    (degree * pi) / 180
}
#---
# change from radians to degrees
rad_to_deg <-  function (radian) {
    (radian * 180) / pi
}
#---
d2km <- function (d, base.latitude = 1) 
{
    if (!requireNamespace("fields")) 
        stop("Required fields package is missing.")
    onerad_to_degree.dist <- fields::rdist.earth(matrix(c(0, base.latitude), ncol = 2), 
                                                 matrix(c(1, base.latitude), ncol = 2), 
                                                 miles = FALSE)[,1]
    out <- d * onerad_to_degree.dist
    return(out)
}



#--------------
# spatial gradient
# function to find the spatial gradient of mean temperature across space.
spatial_grad_JL <- function(rx, 
                            y_diff = 1, 
                            unit_out = "km") {
    
    if(!(unit_out == "km" | unit_out == "m")){
        stop("unit_out should be 'km' or 'm'")
    }
    
    # if there's more than one layer (i.e., input data weren't averaged), average them.
    if(nlyr(rx) > 1){ rx <- mean(rx,na.rm = TRUE) }
    
    # find y-distance between cells
    if (.getProj(rx) == 'longlat') {
        if(unit_out=="km"){
            y_dist <- d2km(res(rx)) # from degrees to km
        }
        if(unit_out=="m"){
            y_dist <- d2km(res(rx))*1000 # from degrees to m
        }
        
    } else {
        if(unit_out=="km"){
            y_dist <- res(rx) / 1000 # from meters to km
        }
        if(unit_out=="m"){
            y_dist <- res(rx)
        }
        y_diff <- NA
    }
    
    if (!.is_package_installed("dplyr") || !.is_package_installed('tidyr')) stop('The packages dplyr and tidyr are needed for this metric; Please make sure they are installed!')
    
    # find cells that touch each other
    y <- data.frame(adjacent(rx, cells=1:ncell(rx), directions=8,pairs=TRUE))
    y <- y[order(y$from, y$to),] # sort 
    y <- na.omit(y) # remove NAs (but there shouldn't be any)
   
    
    # add mean temp from each "to" cell
    y$temp <- rx[y$to][,1]
    # find the difference in rows between adjacent cells
    # -1 = focal cell is above; 0 = same row; 1 = focal cell is below
    y$sy <- rowFromCell(rx, y$from)-rowFromCell(rx, y$to)
    # find column difference between adjacent cells
    # -1 = focal cell is right; 0 = same column, 1 = focal cell is left
    y$sx <- colFromCell(rx, y$to)-colFromCell(rx, y$from)
    y$sx[y$sx > 1] <- -1
    y$sx[y$sx < -1] <- 1
    # NOTE: these last 2 lines are inherited from climetrics, but 
    # shouldn't do anything, as from - to should never be >1 or <-1 
    
    # paste the coordinates of the col and row difference
    y$code <- paste(y$sx, y$sy)
    
    # translate coordinates to plain English
    y$code1 <- eval(parse(text='dplyr::recode(y$code,
                           `1 0` = "tempE",
                           `-1 0` = "tempW",
                           `-1 1` = "tempNW",
                           `-1 -1` = "tempSW",
                           `1 1` = "tempNE",
                           `1 -1` = "tempSE",
                           `0 1` = "tempN",
                           `0 -1` = "tempS")'),envir =environment())
    
    # remove +1, 0, -1 directional columns, keep only plain english. 
    y3b <- eval(parse(text="dplyr::select(y,from, code1, temp)"),envir =environment())
    # spread so each focal cell is one row
    y3b <- eval(parse(text="tidyr::spread(y3b,code1, temp)"),envir =environment())
    # add temperature of the "from" cells. 
    y3b$tempFocal <- rx[y3b$from][,1]
    # add latitude (used to correct for differences in cell area across latitudes.)
    y3b$LAT <- yFromCell(rx, y3b$from)
    
    if(!is.na(y_diff)) {
        # add a latitudinal corrector for cell area for the 
        # focal cell, the cell above, and the cell below. 
        # basically, this adds a correction coefficient using cos of latitude
        # where cos(deg_to_rad(0)) = 1, and scales in an arch to cos(deg_to_rad(90)) = 0
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = cos(deg_to_rad(LAT + y_diff)),
                         latneg = cos(deg_to_rad(LAT - y_diff)),
                         latfocal = cos(deg_to_rad(LAT)))"),envir =environment())
    } else {
        
        y3b <- eval(parse(text="dplyr::mutate(y3b,
                         latpos = 1,
                         latneg = 1,
                         latfocal = 1)"),envir =environment())
    }
    
    
    
    # calculate the difference in temperature between 9 neighboring cells.
    # here, we multiply the longitude demonimators (the distance between cells)
    # by the radian position of latitude, (e.g., at the pole, cells are 0m 
    # apart, at the equator, they are y_dist apart)
    # note that:
    # NS gradients are (n - s), so a positive value means north is warmer.
    # WE gradients are (e - w), so positive values mean east is warmer. 
    y3c <- "dplyr::mutate(y3b,
                       gradWE1 = (tempN-tempNW)/ (latpos *  y_dist[1]),
                       gradWE2 = (tempFocal - tempW)/(latfocal * y_dist[1]),
                       gradWE3 = (tempS-tempSW)/(latneg * y_dist[1]),
                       gradWE4 = (tempNE-tempN)/(latpos * y_dist[1]),
                       gradWE5 = (tempE-tempFocal)/(latfocal * y_dist[1]),
                       gradWE6 = (tempSE-tempS)/(latneg * y_dist[1]),
                       gradNS1 = (tempNW-tempW)/y_dist[2],
                       gradNS2 = (tempN-tempFocal)/y_dist[2],
                       gradNS3 = (tempNE-tempE)/y_dist[2],
                       gradNS4 = (tempW-tempSW)/y_dist[2],
                       gradNS5 = (tempFocal-tempS)/y_dist[2],
                       gradNS6 = (tempE-tempSE)/y_dist[2])" 
    
    y3c <- eval(parse(text=y3c),envir=environment())
    
    y3c <- eval(parse(text="dplyr::rowwise(y3c)"),envir=environment())
    
    # basically takes a weighted average of gradient values in the cells,
    # with pairwise gradients including the focal cell counted twice. 
    y3c <- eval(parse(text="dplyr::mutate(y3c,
      WEgrad = .mnwm(gradWE1, gradWE2, gradWE3, gradWE4, gradWE5, gradWE6),
      NSgrad = .mnwm(gradNS1, gradNS2, gradNS3, gradNS4, gradNS5, gradNS6),
      angle = .ang(WEgrad, NSgrad))"),envir=environment())
    
    # ok, so in the weighted averages, 
    # positive NSgrad values mean that the North is warmer; negative values mean south is warmer.
    # e.g., NSgrad is the value that temperature changes as you move North!
    # positive WEgrad values mean that the East is warmer; negative values mean West is warmer. 
    # e.g., WEgrad is the value that temperature changes as you move East!  
    y3c <- eval(parse(text="dplyr::select(y3c,icell = from, WE = WEgrad, NS = NSgrad, angle = angle)"),envir=environment())
    
    # JL edit: replace grad and angle with NA if input raster had no value 
    y3c$angle[is.na(values(rx))] <- NA
    y3c$NS[is.na(values(rx))] <- NA
    y3c$WE[is.na(values(rx))]  <- NA
    
    NS <- y3c$NS # isolate NS gradient 
    WE <- y3c$WE # isolate EW gradient
    NS[is.na(NS)] <- 0 # change NAs to 0
    WE[is.na(WE)] <- 0 # change NAs to 0
    NAsort <- ifelse((abs(NS)+abs(WE)) == 0, NA, 1)
    # calculate spatial gradient (hypotneuse of NS and EW gradient)
    y3c$Grad <- NAsort * sqrt((WE^2) + (NS^2))
    # note GRAD is always going to be positive because both values are squared. 
    # which means that we will need to only interpret gradient along with angle
    # in the following steps, we have to consider angle when interpreting gradient
    
    # make output raster
    rx <- c(rx, rx, rx, rx)
    names(rx) <- names(y3c)[-1]
    rx[[1]][y3c$icell] <- y3c$WE
    rx[[2]][y3c$icell] <- y3c$NS
    rx[[3]][y3c$icell] <- y3c$angle
    rx[[4]][y3c$icell] <- y3c$Grad
    # note that in the finished raster, 
    # WE: positive values mean cooler in the west, negative values mean cooler in the east
    # NS: positive values mean cooler in the south, negative values mean cooler in the north
    return(rx)
}



#----------
# truncate is for bounding max and min values to upper (95%) and lower (5%) quantiles, respectively
gVelocity_JL <- function(grad, slope, grad_col = "Grad", truncate=FALSE) {
    
    v <- slope # this is temp trend 
    g <- grad[grad_col] # this is gradient == always positive
    v_ang <- grad["angle"] 
    
    # velocity angles have opposite direction to the spatial climatic gradient if warming and same direction (cold to warm) if cooling
    v_ang <- terra::ifel(slope > 0, v_ang + 180, v_ang) %% 360
    # NOTE: this line is inherited from VoCC/climetrics, but we make some changes below. 
    # this is true: gVel angles should be opposite to the spatial gradient if the trend is warming
    # and same angle as spatial gradient if cooling, 
    # but this strategy omits that the numerator of CV (ttrend) can be negative, 
    # and thus influences the interpretation of the ending vector. 
    # For example:
    # a cell with a downward gradient (angle = 180; points towards warmer), and a warming trend of 1
    # will result in a vector with magnitude 1 and angle 0 (towards North). That is CORRECT. 
    # However, cell with a downward gradient (angle = 180) and a cooling trend of (-1) **should**
    # result in a vector of magnitude 1 and angle 180 (down), but the current code results in a 
    # magnitude -1 and angle 180, which is essentially the same as a northward velocity (mag 1 at angle 0). 
    
    # POTENTIAL FIXES:
    # 1. keep current code (flipping angle when trend is positive), 
    # but then use abs(v)/g to calculate the vector magnitude
    # this way, velocities are always positive, and angles always point in the direction of isotherm movement. 
    # 2. standardize all spgrad to direction of cooling by flipping ALL angles, not just when ttrend is positive
    # since we already flipped half the angles, we'll implement Fix 1 here:
    
    # calculate velocity
    # note that this differs from VoCC/climetrics
    v <- abs(v) / g
    
    # truncate extreme values if truncate=T (if there are lots of cells in the raster)
    if(truncate){
        # find anomalously large velocities 
        # (usually the result of small spatial gradients in the denominator)
        .o <- as.matrix(global(v,fun=quantile,probs=c(0.05,0.95),na.rm=TRUE))[1,]
        # change 0-5% velocities equal to the 5%
        v[v < .o[1]] <- .o[1]
        # change 95-100% velocities to be equal to 95%
        v[v > .o[2]] <- .o[2] 
    }
    # note that in "weird" climate velocity values often occur when cells in the 
    # 9-cell neighborhood used to calculate spatial gradient are NA. Thus, our
    # addition of a buffer around the study polygon helps to reduce "weird" values
    # within the focal polygon by keeping edge cells, and making those cells "sacrificial"
    # (e.g., they won't, or will barely be, averaged in to the study area mean value),
    # instead of having a bunch of weird values driven by neighboring NAs inside our focal polygon. 
    
    # export velocity magnitude and velocity angle layers
    output <- c(v,v_ang)
    names(output) <- c("Vel", "Ang")
    return(output)
}
