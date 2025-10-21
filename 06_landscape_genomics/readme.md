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
   


## Landscape Genetics 2 - 21 October 2025

Today we will be using the conductance maps we made last week to estimate least cost path (LCP) distances between 
populations, and see which model best fits the data.

### 1. Answer question #5 on your worksheet.

### 2. Change the directory to the one you used last week. Load packages.

Load the packages we need for today:

      library(raster)
      library(gdistance)
      library(ecodist)
      library(sf)

### 3. Estimate distances between populations for your conductance layer.

1. Read your conductance raster. If you named it something different, you will need to change the code.

        x <- raster("test.grd")

2. Load the political boundaries and make the Arizona shapefile again: 

        usa <- read_sf("gadm36_USA_shp", "gadm36_USA_1")
        usa <- st_transform(usa, crs=4326)
        arizona <- usa[usa$NAME_1 %in% "Arizona",]

3. Read in sampling points and modify:

        # read in sampling points and re-structure format
        sampling_points <- read.csv("sitta_sampling.csv", header=T, stringsAsFactors=F)
        sampling_points2 <- as.matrix(sampling_points[,3:2])

4. Make sure they loaded properly by plotting. 

        plot(x)
        plot(st_geometry(arizona), add=T)
        points(sampling_points2)
    
5. Create a cost distance object (functions are part of the gdistance package).

        # first step is to make all NA values in the raster to 0 conductance
        x[is.na(x)] <- 0
        
        # next we create a transition layer that includes costs of moving between cells in all directions (not plottable)
        transition_layer <- transition(x, mean, directions=8)
        
        # corrects the transition layer with information based on how close cells are to one another and the coordinate system
        transition_layer2 <- geoCorrection(transition_layer, "c", scl=T)
        
        # create the cost-distance object, calculating LCP for each pairwise comparison of sampling locations
        cost_distance <- costDistance(transition_layer2, sampling_points2)
        
        # look at the object that was created and pairwise LCP distances
        cost_distance
        
### 4. Statistics of the result using multiple regression of distance matrices (MRM).

    # read in genetic differentiation data
    fst <- read.table("distance_all.csv", header=T, sep=",", stringsAsFactors=F)
        
    # convert the cost distance matrix to a vector
    lcp <- as.vector(cost_distance)
        
    # combine your cost_distance object with the fst object
    fst <- cbind(fst, lcp)
    
    # look at the fst table (now includes your LCP data)
    fst
        
    # do the stats!
    mrm_results <- MRM(fst$fst ~ fst$lcp, nperm=100)
        
    # look at the results:
    mrm_results$r.squared
        
### 5. Write your results on the white-board for comparison with everyone else.

### 6. Make a plot of LCP between populations

If we want to visualize the LCP between two populations, we can do so as follows:

    # see what numbers each of the populations are:
    sampling_points
    
    # get the LCP between population 2 and population 5
    plot(x)
    plot(st_geometry(arizona), add=T)
    lines(shortestPath(transition_layer2, sampling_points2[2, ], sampling_points2[5, ], output="SpatialLines"))
   
If you want to plot all the LCP between population 1 and another population, run the following. If you want to change which
population you are starting with, choose a different number assigned to the variable "population_to_use"

    plot(x)
    plot(st_geometry(arizona), add=T)
    population_to_use <- 1
    for(a in 1:9) {
	    lines(shortestPath(transition_layer2, sampling_points2[population_to_use,], sampling_points2[a,], output="SpatialLines"), lwd=2)
    }
    points(sampling_points[,3], sampling_points[,2], cex=1.5, lwd=2)
    points(sampling_points[population_to_use,3], sampling_points[population_to_use,2], pch=19, cex=2)








