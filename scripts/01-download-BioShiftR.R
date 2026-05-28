# 1. get BioShifts polygons from BioShiftR R package

remotes::install_github("bioshifts/BioShiftR", force=T)
library(BioShiftR)
library(dplyr)
library(ggplot2)
library(sf)

# download study polygons from external source (OSF) via BioShiftR
download_polygons()


# In this reprex, we are going to reproduce 2 examples of latitudinal 
# climate velocities: one that is cooling, and one that is warming, 
# both in Sweden:
# - A100 P1 with duration 1974.5-1984.5
# - A100 P1 with duration 1984.5-1994.5

# And also, 3 examples of elevational climate velocities, chosen because
# the study areas are in similar locations and time durations (and here,
# we could not upload the full global temperature data at 1km resulution 
# due to file size constraints on GitHub):
# - A001 P1 with duration 1979-2002
# - A001 P2 with duration 1974-2000
# - A005 P1 with duration 1985-1998

# proceed to scripts 02-demo_ele_cv.R and 03-demo_lat_cv.R