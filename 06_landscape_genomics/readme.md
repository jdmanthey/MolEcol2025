## Landscape Genetics 1 - 16 October 2025

Today we will be working with some geographic data to create conductance maps. Remember, conductance maps are the opposite
of resistance maps. The steps we will be taking today (all steps outlined below):

1. Install necessary packages and download data
2. Plot sampling points on landscape data
3. Look up information about this bird (_Sitta carolinensis_) to formulate an expert opinion
4. Reclassify a landscape feature into a conductance map based on your new knowledge
5. Plot the conductance map and modify or re-make as desired
6. Check with Dr. Manthey about your hypothesis before you leave

On Thursday, we will be creating least-cost path distances across your conductance maps and seeing which hypothesis best
fits the genetic differentiation between populations.

Dr. Manthey will perform two analyses based on:
1. geographic distance between sampling sites (i.e., isolation by distance)
2. overall climatic differences between sampling sites (i.e., isolation by environment)

Data sources:
1. Genetic data: 2015_Manthey_et_al_ME.pdf (included in download directory)
2. Tree cover data: http://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.6.html
3. Land cover data: https://swregap.org/data/landcover/
4. Elevation data: http://www.worldclim.org/v2/#/version1
5. NDVI data: https://search.earthdata.nasa.gov/
6. Political boundaries: https://gadm.org/data.html

### 1. Install Packages

We will need four new packages for the work today: raster, rgdal, ecodist, and gdistance.

    install.packages("raster")
    install.packages("gdistance")
    install.packages("ecodist")
    install.packages("sf")
    
Load the new libraries to make sure the installs worked. Red text is fine, just check there are no errors. 

    library(raster)
    library(gdistance)
    library(ecodist)
    library(sf)

Download the .zip file from this GitHub directory, unzip it, and set that as your new working directory for today. We will be
using the same working directory on Friday.

### 2. Plot sampling points on landscape data

First, we will read in the raster data. Rasters are a way of visualizing continuous data on a landscape by breaking up data
into a grid. Here, the grid size is 30 arc seconds (https://en.wikipedia.org/wiki/Minute_and_second_of_arc).

    elevation <- raster("elevation_resample.grd")
    landcover <- raster("landcover_resample.grd")
    treecover <- raster("tree_resample.grd")
    ndvi <- raster("ndvi_resample.grd")
    
These layers are (1) elevation, (2) land cover classification into multiple categories, (3) % tree cover, and (4) NDVI 
(https://en.wikipedia.org/wiki/Normalized_difference_vegetation_index).

Next, we will read in a shapefile of USA political boundaries. A shapefile is a way to represent points and polygons in 
geographic space. For those of you that have taken a geographic information systems (GIS) class, you know about coordinate
systems. All the files here are pre-processed to have the same geographic coordinate system 
(https://desktop.arcgis.com/en/arcmap/latest/map/projections/about-geographic-coordinate-systems.htm). 

    usa <- read_sf("gadm36_USA_shp", "gadm36_USA_1")
    usa <- st_transform(usa, crs=4326)
    arizona <- usa[usa$NAME_1 %in% "Arizona",]

Now we will read in the genetic sampling locations from the study:

    # read in sampling points and re-structure format
    sampling_points <- read.csv("sitta_sampling.csv", header=T, stringsAsFactors=F)
    sampling_points2 <- as.matrix(sampling_points[,3:2])

And finally, plot each of the rasters, with the sampling points overlaid:

Elevation:

    plot(elevation, main="elevation")
    plot(st_geometry(arizona), add=T)
    points(sampling_points2)

Land Cover. Here we define specific color classes for the different classifications (lightgray = barren, darkgreen = forest, 
lightgreen = scrub, yellow = herbaceous, blue = water, magenta = developed, red = disturbed)
    
    colors <- c("lightgray", "darkgreen", "lightgreen", "yellow", "blue", "magenta", "red")
    plot(landcover, col=colors, main="landcover")
    plot(st_geometry(arizona), add=T)
    points(sampling_points2)

Tree Cover %:

    plot(treecover, main="tree cover")
    plot(st_geometry(arizona), add=T)
    points(sampling_points2)
    
NDVI:

    plot(ndvi, main="ndvi")
    plot(st_geometry(arizona), add=T)
    points(sampling_points2)

### 3. Look up information about this bird (_Sitta carolinensis_) to formulate an expert opinion about factors shaping connectivity:

General info: https://en.wikipedia.org/wiki/White-breasted_nuthatch

More: https://www.allaboutbirds.org/guide/White-breasted_Nuthatch/overview

This requires an account to view, but is an interesting resource if you want to check it out:
Distribution in breeding season in Arizona over past few years (zoom in to see specific points): 
https://ebird.org/map/whbnut?neg=true&env.minX=-119.43028641658736&env.minY=30.52623221873654&env.maxX=-106.09288407283736&env.maxY=36.502502711573584&zh=true&gp=true&ev=Z&excludeExX=false&excludeExAll=false&mr=on&bmo=6&emo=6&yr=range&byr=2023&eyr=2025

Sampling locality names: Check out the map in the pdf distributed here: 2015_Manthey_et_al_ME.pdf

Hint: Compare with satellite imagery in Google Maps to help get ideas of what shapes the breeding distribution of this bird.

### 4. Reclassify a landscape feature into a conductance map based on your new knowledge

Next, you will need to reclassify a raster of your choice (choose 1) into a conductance map. All conductance values should be
between 0 and 1. If you want to get more complicated, you can combine multiple rasters into a single conductance map.

Examples of how to do this (but don't use these values specifically, they are made-up). These examples use a ???? to define
the test object. You should change that to one of the names of the raster layers (e.g., treecover, ndvi, elevation):

    # Define the test raster
    test <- ????
    # change all values lower than 1500 to a conductance of 0.2
    values(test)[values(test) < 1500] <- 0.2
    # change all values greater than or equal to 1500 and less than or equal to 2500 to a conductance of 0.8
    # if you want another range of values, just copy the line below and add another range of values
    values(test)[values(test) >= 1500 & values(test) <= 2500] <- 0.8
    # change all values greater than 2500 to a conductance of 0.5
    values(test)[values(test) > 2500] <- 0.5
    # plot your new test conductance map
    plot(test)

Examples for landcover reclassification. Make sure to redefine all values otherwise the map will not make sense.

    # initial values: 1000 = barren, 2000 = forest, 3000 = scrub, 4000 = herbaceous, 
    #5000 = water, 6000 = developed (i.e., urban), 7000 = disturbed (e.g., agriculture, invasive species, destroyed areas)
    test <- landcover
    values(test)[values(test) == 1000] <- 0.1
    values(test)[values(test) == 2000] <- 0.5
    values(test)[values(test) == 3000] <- 0.7
    values(test)[values(test) == 4000] <- 0.2
    values(test)[values(test) == 5000] <- 0.4
    values(test)[values(test) == 6000] <- 0.3
    values(test)[values(test) == 7000] <- 0.1
    plot(test)

Now, create a conductance map that you think is reasonable.

### 5. Plot your conductance map. Remake if necessary.

    plot(test)

Redo the classification of the layer if you don't like the results.

### 6. Check with Dr. Manthey about your scheme. 

If the test object looks good and reasonable, we can write it to a new file:

    writeRaster(test, file="test")

If approved, save a pdf file of the plot of your conductance raster, and email to Dr. Manthey to be shown next week.
   


