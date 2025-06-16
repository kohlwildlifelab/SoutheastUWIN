##########################
#
# This code is used in Lavendar et al., " Species Richness Response to Urbanization across mid-size cities in the  Southeastern United-States"
#
# Landcover and landscapemetrics for a city (repeat code for each city)
#
# Lavendar Harris and Michel T. Kohl 
#
#########################

# Load packages 
library(landscapemetrics)
library(raster)
library(glue)
library(dplyr)
library(tidyverse)
library(sf)
library(tidycensus)
library(mapview)
library(tigris)
library(crsuggest)
library(reshape2)
library(rgdal)
library(ggplot2)
library(viridis)
library(terra)

#### Load data for Athens

# Load OSM layer for Athens
athnlcd <- raster::raster("./OSM/Athens_OSM_5.tif")

# Plot
plot(athnlcd)

# Check that raster can be used in landscapemetrics package

check_landscape(athnlcd)

# Load camera locations

sites <- read.csv("./Athens_Camera_Sites.csv")
head(sites)

# Create a vector of locations
s <- vect(sites, c("Long", "Lat"), crs="+proj=longlat")

# Project camera locations using raster crs
ath.p <- project(s, crs(athnlcd))
crs(ath.p) <- "+proj=utm +zone=17 +datum=WGS84"

# Plot locations over raster
plot(athnlcd);points(ath.p)


#### Plot for manuscript 

# Convert raster to dataframe
landuse_df <- as.data.frame(athnlcd, xy = TRUE)

# Convert camera location vector to sf
gps_sf <- st_as_sf(sites, coords = c("Long", "Lat"), crs = 4326)

# Transform camera locations to match raster CRS
gps_sf <- st_transform(gps_sf, crs(athnlcd))

# Make raster a factor

landuse_df$Athens_OSM_5 <- as.factor(landuse_df$Athens_OSM_5)

# Pull out classified layers

levels(landuse_df$Athens_OSM_5) <- c("Urban development", "Unclassified development", "Natural areas", "Water", "Physical infrastructure")

# Plot using ggplot2

Athens_map <- ggplot() +
  geom_raster(data = landuse_df, aes(x = x, y = y, fill = Athens_OSM_5)) + # Adjust fill column
  scale_fill_viridis(name = "Land Cover Type", discrete = TRUE)+
  geom_sf(data = gps_sf, shape = 21, fill = "red", color = "black", size = 2)+ # Overlay GPS points
  labs(title = "Athens, Georgia")+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5, size=16, face='bold'),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#ggsave("Athens_map.jpeg", height = 5, width = 7, units = c("in"))

#### Landscape metrics analysis 

# Check metrics

metric.names<-lsm_abbreviations_names #store names
tail(metric.names) #show last 6 observations
view(metric.names) 

# Establish buffers
buffer_sizes=c(500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500)

# Buffer around camera locations

buffer_ath = buffer_sizes %>%
  set_names()%>%
  map_dfr(~sample_lsm(list(athnlcd),
                      what=c("lsm_p_area", "lsm_l_ed", "lsm_l_shei",
                             "lsm_c_area_sd", "lsm_c_te", "lsm_l_te"),
                      y=sites_sf,
                      size=.),
          .id="buffer")

# Look at current order of buffer

levels(factor(buffer_ath$buffer))

# Rename the classes based off of OSM NLCD class legend/description
unique(buffer_ath$class)

buffer_ath$class_name<- plyr::revalue(factor(buffer_ath$class),
                                      c("1"="Landuse", "2"="Developed", "3"="Natural areas",
                                        "4"="Water", "5"="Infrastructure"))

## Filter different metrics from results

# Across landscape

landscape_results<-filter(buffer_ath,level=="landscape")

# Across classes

class_results<-filter(buffer_ath,level=="class")

# Across patches

patch_results<-filter(buffer_ath,level=="patch")

# Create a dataframe with appropriate class buffer results 
ath_class <- dcast(class_results, buffer+plot_id~metric, value.var="value", fun.aggregate=mean)

# Rename covariates to differentiate layers 

ath_class<-ath_class %>%
  rename(c_area_sd=area_sd,
         c_te=te)

# Create a dataframe with appropriate landscape buffer results 

ath_landscape <- dcast(landscape_results, buffer+plot_id~metric, value.var="value", fun.aggregate=mean)

# Rename covariates

ath_landscape<-ath_landscape%>%
  rename(l_ed=ed,
         l_shei=shei,
         l_te=te)

# Create a dataframe with appropriate patch buffer results 

ath_patch <- dcast(patch_results, buffer+plot_id~metric, value.var="value", fun.aggregate=mean)

# Rename covariates

ath_patch<-ath_patch %>%
  rename(p_area=area)

# Join the covariates together 
dat<-right_join(ath_landscape,ath_class, by=c("buffer","plot_id"))
dat2<-right_join(dat,ath_patch, by=c("buffer","plot_id"))

#write.csv(x=dat2,file="Final_athnlcd_analysis.csv")

#
#
#
#### Load in Shannon Diversity Index data and merge with landscape metrics data
sdi <- read.csv("species_index.csv")
head(sdi)
out <- read.csv("Final_athnlcd_analysis.csv")
head(out)

output <- merge(out, sdi, by = "plot_id")
#write.csv(x=output, file="athens_frag_data.csv")

