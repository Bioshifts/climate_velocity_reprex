# 1. get BioShifts polygons from BioShiftR R package

devtools::install_github("bioshifts/BioShiftR", force=T)
library(BioShiftR)
library(dplyr)
library(ggplot2)
library(sf)

# download study polygons from external source (OSF) via BioShiftR
download_polygons()

# Identify a few polygons for reproducing:
# We'll find 2 latitudinal and 2 elevational polygons with 
# different trends (warming vs cooling), 
# and with short enough duration that we don't have to upload
# all the temperature data to this GitHub repo (< ~20 years I guess. Arbitrary.), 
# and small enough size to run locally on a laptop. 

get_shifts(type = "ELE") %>%
    add_trends() %>%
    add_methods() %>% filter(duration < 25) %>%
    filter(trend_temp_mean < 0) %>%
    distinct(article_id, poly_id, duration, midpoint_firstperiod, midpoint_lastperiod) %>%
    add_poly_info() %>%
    filter(area_km2 < 10000) %>%
    #select(midpoint_firstperiod, midpoint_lastperiod, duration) %>%
    slice(c(1,7)) %>%
    

    #ggplot() +  geom_histogram(aes(x = area_km2))
    add_polygons() %>%
    ggplot() + 
    geom_sf() + 
    geom_sf_text(aes(label = study_area))

get_shifts() %>%
    filter(article_id == "A041") %>%
    add_methods() %>%
    distinct(midpoint_firstperiod, midpoint_lastperiod, duration)

get_shifts(type = "ELE") %>%
    add_methods() %>%
    filter(duration < 25) %>%
    distinct(article_id, poly_id, eco, param, type, method_id, duration, midpoint_firstperiod, midpoint_lastperiod) %>%
    slice(1:2) %>%
    add_polygons() %>%
    add_trends() %>%
    pull(trend_temp_mean)
    add_poly_info() %>% pull(study_area)
    ggplot() +
    geom_sf() 

all_ele_shifts <- get_shifts() %>%
    filter(type == "ELE") %>%
    # add methods so we can get duration
    add_methods() %>%
    # find distinct articles, polygons, and durations
    distinct(article_id, poly_id, type, method_id, eco, duration)

# find 

all_ele_trends <- all_ele_shifts %>% 
    add_trends() 

all_ele_trends %>%
    filter(trend_temp_mean > 0,
           duration < 25) %>%
    add_polygons() %>%
    ggplot() +
    geom_sf()



# latitudinal examples ----------------------------------------------------

get_shifts(type = "LAT") %>%
    distinct(article_id, poly_id) %>%
    slice(3:5) %>%
    add_polygons() %>% 
    ggplot() + 
    geom_sf() +
    geom_sf_label(aes(label = paste0(article_id,"_", poly_id)))

get_shifts(type = "LAT") %>%
    filter(article_id %in% c("A006","A008")) %>%
    add_trends() %>%
    add_methods() %>%
    distinct(article_id, trend_temp_mean, duration, midpoint_firstperiod, midpoint_lastperiod)



get_shifts() %>%
    add_methods() %>% 
    filter(duration < 25) %>%
    distinct(article_id, poly_id) %>%
    add_poly_info() %>%
    filter(stringr::str_detect(study_area,"Sweden|Finland")) %>%
    pull(article_id) -> test

get_shifts() %>%
    filter(article_id %in% test) %>%
    add_methods() %>%
    filter(duration < 25) %>%
    add_trends() %>%
    distinct(article_id, poly_id, trend_temp_mean, midpoint_firstperiod, midpoint_lastperiod)

get_shifts() %>%
    filter(article_id %in% c("A006","A098","A117","A190")) %>%
    add_methods() %>%
    add_trends() %>%
    distinct(article_id, midpoint_firstperiod, midpoint_lastperiod, duration, trend_temp_mean)


get_shifts(type = "ELE") %>%
    add_trends() %>%
    add_methods() %>%
    distinct(article_id, poly_id, duration, midpoint_firstperiod, midpoint_lastperiod) %>%
    mutate(art_poly_id)
    mutate(article)

    
get_shifts(type = "LAT") %>%
    filter(article_id == "A100") %>%
    distinct(article_id, poly_id) %>%
    add_polygons() %>%
    ggplot() + geom_sf()
    

get_shifts(type = "LAT") %>%
    filter(article_id == 'A100') %>%
    add_trends() %>%
    add_methods()  %>%
    distinct(midpoint_firstperiod, midpoint_lastperiod, 
             duration, trend_temp_mean)

