
##########################
#
# This code implements the anaylsis used in Lavendar et al., " Species Richness Response to Urbanization across mid-size cities in the  Southeastern United-States"
#
# Formatting Data and Quantifying Scales for Generalized Linear Mixed Effect Model
#
# Lavendar Harris and Michel T. Kohl 
#
# 
#########################

## Load packages

library(tidyverse)
library(dplyr)
library(usdm)
library(lme4)
library(MuMIn)
library(ggplot2)
library(patchwork)


## Format data for analysis
# Load city scale data

city_scale <- read.csv("./data/covariates/city_scale_data.csv")

# Load landscape metrics data

city_landscape <- read.csv("./data/covariates/city_landscape_data.csv")

# Load city point values data

# Load city point values data and replace NAs with 0 

city_pointvalues <- read.csv("./data/covariates/city_point_values_data.csv")
city_pointvalues[is.na(city_pointvalues)] <- 0

## Replace the city scale data NA's with point data values

currentdf <- city_scale
pointdf <- city_pointvalues

# Grab the data columns for all covariates

columnlist <- c(5:16)

# Blank dataframe to put new data on

city_df <- data.frame()

## Run loop to replace NA's

for(i in columnlist){
  # Grab one column at a time
  small.df <- currentdf[,c(2:4,17:20, i)]
  small.point.df <- pointdf[,c(2, i-2)] 
  
  head(small.df)
  head(small.point.df)
  
  # Grab name of column need to replace na's with
  name <- colnames(small.df[8])
  name
  
  corrected.df <- small.df %>% inner_join(small.point.df, by= "Site")
  
  corrected.df$NewValue <- ifelse(is.na(corrected.df[,8]), corrected.df[,9], corrected.df[,8])
  
  corrected.df <- corrected.df %>%
       dplyr::select(c(Site, Lat, Long, Scale, RICHNESS, SHANNON, StudyArea, NewValue))
  
  colnames(corrected.df)[ncol(corrected.df)] <- name
  
  if(nrow(city_df)<1){
    city_df <-corrected.df
  } else {city_df <- cbind(city_df, corrected.df[,ncol(corrected.df)])
  colnames(city_df)[ncol(city_df)] <- name}
}

city_df

## Combine the city scale and city landscape data to create one large dataframe
# Need to change the scale names from each dataset so they match

unique(city_df$Scale)
unique(city_landscape$buffer)

# Change scale names 

city_df <- city_df %>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 500 _m", "500", as.character(Scale))) %>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 600 _m", "600", as.character(Scale))) %>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 700 _m", "700", as.character(Scale))) %>% 
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 800 _m", "800", as.character(Scale))) %>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 900 _m", "900", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1000 _m", "1000", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1100 _m", "1100", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1200 _m", "1200", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1300 _m", "1300", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1400 _m", "1400", as.character(Scale)))%>%
  mutate(Scale = ifelse(as.character(Scale) == "Buffer_ 1500 _m", "1500", as.character(Scale)))


# Change from "scale" to "buffer" to merge datasets

colnames(city_df)[4]<- "buffer"

# Make buffer a number

city_df$buffer  <- as.integer(city_df$buffer)

# Combine our two datasets to create a final dataframe

final_data <-merge(city_df, city_landscape, by=c("Site","buffer", "RICHNESS", "SHANNON", "StudyArea"))

# Check dataframe

head(final_data)

# Standardize covariates for model

standat <- final_data %>% mutate_at(c('imperv.mean', 'tree.mean', 'income_numE',
                               'density.km2', 'housing.km2', 'p_white', 'p_black',
                               'l_ed','l_te', 'c_area_sd', 'c_te', 'p_area'), ~(scale(.) %>% as.vector))

# Check
standat

### AIC for habitat and socioeconomic spatial scales for all cities within model

# Pull out scales 

scales <- unique(final_data$buffer)[c(1:11)]
scales

# AIC for Habitat scale

aic_tabs <- data.frame()
for(i in scales){
  scale <- subset(standat, buffer == i)
  mod <- lmer(SHANNON ~ imperv.mean + tree.mean + l_ed + (1 |StudyArea), data=scale) # habitat covariates
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs <- rbind(table, aic_tabs)
}

habitat <- aic_tabs
habitat <- habitat[order(habitat$AIC),]
habitat

# library(writexl)
# write_xlsx(habitat, "habitat.xlsx")

# AIC for Socioeconomic scale

aic_tabs1 <- data.frame()
for(i in scales){
  scale <- subset(standat, buffer == i)
  mod <- lmer(SHANNON ~ income_numE + density.km2 + p_white + (1 |StudyArea), data = scale) # change this to the socioeconomic variables
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs1 <- rbind(table, aic_tabs1)
}

socio <- aic_tabs1
socio <- socio[order(socio$AIC),]
socio

#write_xlsx(socio, "socio.xlsx")

## Separate and combine scales (determined by AIC) into a datframe to use for GLM model

analysis.df <- standat[which(standat$buffer == "1000"), c(1:2, 4:5, 18:19,21,23,26)]
analysis.df2 <- standat[which(standat$buffer == "500"), c(1:2, 10, 12, 14)]

analysis.df3 <- merge(analysis.df, analysis.df2, by = "Site")

######### MODEl ANALYSIS

# Full GLMM model

dat.fullmodel <- lmer(SHANNON ~ imperv.mean + tree.mean + l_ed + p_area + income_numE + density.km2 + p_white + (1 |StudyArea), data = analysis.df3, REML = F)

options(na.action = "na.fail") # Required for dredge to run

# Dredge model

fullmodel.dredge <- dredge(dat.fullmodel, beta = F, evaluate = T, rank = AICc)

head(fullmodel.dredge)

# Model averaging for full GLMM model

summary(model.avg(fullmodel.dredge, subset=delta<=2))

#
#
#
######### Analysis for appropriate spatial scale for each individual city
# AIC for Habitat and socioeconomic covariates by each city

# Athens habitat AIC

at <- subset(standat, StudyArea  == "Athens")
aic_tabs4 <- data.frame()
for(i in scales){
  scale <- subset(at, buffer == i)
  mod <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed, data=scale)
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs4 <- rbind(table, aic_tabs4)
}

ath_habitat <- aic_tabs4
ath_habitat <- ath_habitat[order(ath_habitat$AIC),]
ath_habitat

# Athens socioeconomic AIC

at1 <- subset(standat, StudyArea  == "Athens")
aic_tabs5 <- data.frame()
for(i in scales){
  scale <- subset(at1, buffer == i)
  mod <- glm(SHANNON ~ income_numE + density.km2 + p_white, data = scale) ### change this to the socioeconomic variables
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs5 <- rbind(table, aic_tabs5)
}

ath_socio <- aic_tabs5
ath_socio <- ath_socio[order(ath_socio$AIC),]
ath_socio

# Select and combine scales determined by AIC

adf <- at[which(at$buffer == "1200"), c(1:2, 4:5, 18:19,21,23,26)]
adf2 <- at[which(at$buffer == "600"), c(1:2, 10, 12, 14)]

adf3 <- merge(adf, adf2, by = "Site")

# Jackson habitat AIC

jk <- subset(standat, StudyArea  == "Jackson")
aic_tabs7 <- data.frame()
for(i in scales){
  scale <- subset(jk, buffer == i)
  mod <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed, data=scale)
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs7 <- rbind(table, aic_tabs7)
}

jk_habitat <- aic_tabs7
jk_habitat <- jk_habitat[order(jk_habitat$AIC),]
jk_habitat

# Jackson socioeconomic AIC

jk1 <- subset(standat, StudyArea  == "Jackson")
aic_tabs8 <- data.frame()
for(i in scales){
  scale <- subset(jk1, buffer == i)
  mod <- glm(SHANNON ~ income_numE + density.km2 + p_white, data = scale) ### change this to the socioeconomic variables
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs8 <- rbind(table, aic_tabs8)
}

jk_socio <- aic_tabs8
jk_socio <- jk_socio[order(jk_socio$AIC),]
jk_socio

# Separate and combine scales based off AIC
''
jdf <- jk[which(jk$buffer == "1100"), c(1:2, 4:5, 18:19,21,23,26)]
jdf2 <- jk[which(jk$buffer == "500"), c(1:2, 10, 12, 14)]

jdf3 <- merge(jdf, jdf2, by = "Site")

# Little Rock habitat AIC

lr <- subset(standat, StudyArea  == "Little Rock")
aic_tabs3 <- data.frame()
for(i in scales){
  scale <- subset(lr, buffer == i)
  mod <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed, data=scale)
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs3 <- rbind(table, aic_tabs3)
}

lr_habitat <- aic_tabs3
lr_habitat <- lr_habitat[order(lr_habitat$AIC),]
lr_habitat

# Little Rock socioeconomic AIC

lr1 <- subset(standat, StudyArea  == "Little Rock")
aic_tabs6 <- data.frame()
for(i in scales){
  scale <- subset(lr1, buffer == i)
  mod <- glm(SHANNON ~ income_numE + density.km2 + p_white, data = scale) ### change this to the socioeconomic variables
  table <- data.frame(buffer = i,
                      loglik = logLik(mod),
                      AIC = AIC(mod))
  aic_tabs6 <- rbind(table, aic_tabs6)
}

lr_socio <- aic_tabs6
lr_socio <- lr_socio[order(lr_socio$AIC),]
lr_socio

# Separate and combine scales based off AIC

ldf <- lr[which(lr$buffer == "500"), c(1:2, 4:5, 18:19,21,23,26)]
ldf2 <- lr[which(lr$buffer == "500"), c(1:2, 10, 12, 14)]

ldf3 <- merge(ldf, ldf2, by = "Site")


######### Model analysis for individual cities
#
#
# Athens GLM model
model.phys.athens <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed + income_numE + density.km2 + p_white , data = adf3)

# Summary
summary(model.phys.athens)

# Confidence intervals
c_athens <- confint.default(model.phys.athens)

# Jackson GLM model
model.phys.jackson <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed + income_numE + density.km2 + p_white , data = jdf3)

# Summary
summary(model.phys.jackson)

# Confidence intervals
c_jack<- confint.default(model.phys.athens)

# Little Rock GLM model
model.phys.lr <- glm(SHANNON ~ imperv.mean + tree.mean + l_ed +income_numE + density.km2 + p_white, data = ldf3)

# Summary
summary(model.phys.lr)

# Confidene intervals
c_lrock<- confint.default(model.phys.athens)


######## FINAL PLOTS 
#
#
# Plot covariates based off full model average

# % Impervious surface

imperv <-ggplot(data = analysis.df3, aes(x = imperv.mean,y = SHANNON, colour=StudyArea)) +
  geom_smooth(method=lm)+
  scale_color_manual(values=c("Athens"="#920000", "Jackson"="#009e73", "Little Rock"="#006ddb"), drop = F)+
  theme_classic() +
  labs(x='Mean Impervious Surface', y='Shannon Index', title='Response to Impervious Surface') +
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold'))

imperv

# ggsave("Cities_imperv.jpeg", height = 5, width = 7, units = c("in"))

# Edge density

edge <-ggplot(data = analysis.df3, aes(x = l_ed,y = SHANNON, colour=StudyArea)) +
  geom_smooth(method=lm)+
  scale_color_manual(values=c("Athens"="#920000", "Jackson"="#009e73", "Little Rock"="#006ddb"), drop = F)+
  theme_classic() +
  labs(x='Edge Density', y='Shannon Index', title='Response to Edge Density') +
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold'))

edge

# ggsave("Cities_edge.jpeg", height = 5, width = 7, units = c("in"))

# Tree canopy cover

tree <-ggplot(data = analysis.df3, aes(x = tree.mean,y = SHANNON, colour=StudyArea)) +
  geom_smooth(method=lm)+
  scale_color_manual(values=c("Athens"="#920000", "Jackson"="#009e73", "Little Rock"="#006ddb"), drop = F)+
  theme_classic() +
  labs(x='Tree Canopy Cover', y='Shannon Index', title='Response to Tree Canopy Cover') +
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold'))

tree

# ggsave("Cities_tree.jpeg", height = 5, width = 7, units = c("in"))

# Human population density
den<-ggplot(data = analysis.df3, aes(x = density.km2,y = SHANNON, colour=StudyArea)) +
  geom_smooth(method=lm)+
  scale_color_manual(values=c("Athens"="#920000", "Jackson"="#009e73", "Little Rock"="#006ddb"), drop = F)+
  theme_classic() +
  labs(x='Human Density', y='Shannon Index', title='Response to Human Density') +
  theme(axis.title.x = element_text(size = 10, face = "bold"))+
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = "bold"))+
  theme(plot.title = element_text(hjust=0.5, size=15, face='bold'))

den

# ggsave("Cities_den.jpeg", height = 5, width = 7, units = c("in"))

# Put plots together
patchwork.plot <- ((imperv.line | edge.line) /
                     (tree.line| den)) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect', axes = 'collect') &
  theme(legend.position ="bottom")

# ggsave("full_plot.jpeg", dpi = 400, width = 12, height = 7, units = "in")


