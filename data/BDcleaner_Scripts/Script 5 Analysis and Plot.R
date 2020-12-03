#####################################################################################
## Script 5 Analysis and Plot
## Author: Jing Jin, Jun Yang*
## Manuscript: BDcleaner: a workflow for cleaning taxonomic and geographic errors
##             in occurrence data archived in biodiversity databases
######################################################################################

# load packages
library(sp)
library(rgdal)
library(maptools)
library(biogeo)
library(data.table)
library(raster)
library(ggplot2)
library(rworldxtra)
library(rworldmap)
library(data.table)
library(magrittr)

load("YourPath/Cleandataset.RData")
worldmap <- data(wrld_simpl)

# Number of occurrence records in different countries
tb.countrycode <- table(df.bio$countryCode)
tb2 <- sort(tb.countrycode,decreasing = T)
vec.num <- rep(NA,nrow(wrld_simpl@data))
for(i in 1:length(tb.countrycode))
{
  x.row <- which(wrld_simpl@data$ISO2 == names(tb.countrycode)[i])
  vec.num[x.row] <- tb.countrycode[i]
}
wrld_simpl@data$NumOCC <- vec.num
writeOGR(wrld_simpl,'YourPath','YourFileName',driver="ESRI Shapefile", overwrite_layer=TRUE)


df.bio$Overland %>% table()
df.bio$overUrbanType %>% table()

# Number of species and number of occurrence records in each database 
load("YourPath/EightDatabase.RData")

Name.ala <- df.ALA$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.ala  <- get.speNo(Name.ala,use.spe.list)
nrow(df.ALA)

Name.idig <- df.iDig$scientificname %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.idig  <- get.speNo(Name.idig,use.spe.list)
nrow(df.iDig)

Name.rainbio <- df.rainbio$species %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.rainbio  <- get.speNo(Name.rainbio,use.spe.list)
nrow(df.rainbio)

Name.biotime <- df.BioTime$spe %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.biotime  <- get.speNo(Name.biotime,use.spe.list)
nrow(df.BioTime)

Name.spelink <- df.SpeLink$scientificname %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.spelink  <- get.speNo(Name.spelink,use.spe.list)
nrow(df.SpeLink)

Name.bison <- df.Bison$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.bison  <- get.speNo(Name.bison,use.spe.list)
nrow(df.Bison)

Name.gbif <- df.GBIF$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.gbif  <- get.speNo(Name.gbif,use.spe.list)
nrow(df.GBIF)

Name.bien <- df.BIEN$verbatim_scientific_name %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.bien  <- get.speNo(Name.bien,use.spe.list)
nrow(df.BIEN)

Name.tot  <- df.bio$verbatim_scientific_name %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num.tot   <- get.speNo(Name.tot ,use.spe.list)
nrow(df.bio)


###################################
## Taxonomy dimension analysis
###################################

## Percentage of different names' type
table(df.bio$TaxonType)
table(df.bio$TaxonType)/nrow(df.bio)

## Three levels of cleaning
vec1 <- rep(NA,3)
vec2 <- rep(NA,3)

# Match name type
tax.mat <- c("Exact Match","LD Match_Type1","LD Match_Type2","Processing Exact")

# High level
subuse1 <- subset(df.bio, (TaxonStatus %in% c('Accepted','Synonym')) 
                  &(TaxonConf == 3)  & (TaxonType %in% tax.mat) )

Name1 <- subuse1$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num1  <- get.speNo(Name1,use.spe.list)
num1
vec1[1] <- nrow(subuse1); vec2[1] <- nrow(subuse1)/nrow(df.bio)
rm(subuse1)

# Medium level
subuse2 <- subset(df.bio, (TaxonStatus %in% c('Accepted','Synonym')) 
                  &(TaxonConf == 2)  & (TaxonType %in% tax.mat) )
Name2 <- subuse2$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num2  <- get.speNo(Name2,use.spe.list)
num2
vec1[2] <- nrow(subuse2); vec2[2] <- nrow(subuse2)/nrow(df.bio)
rm(subuse2)

# Low level
subuse3 <- subset(df.bio, (TaxonStatus %in% c('Accepted','Synonym')) 
                  &(TaxonConf == 1)  & (TaxonType %in% tax.mat) )
Name3 <- subuse3$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num3  <- get.speNo(Name3,use.spe.list)
num3
vec1[3] <- nrow(subuse3); vec2[3] <- nrow(subuse3)/nrow(df.bio)
rm(subuse3)

# output relative results
df.tax <- data.frame(Num = vec1, Per = vec2, Spe = c(num1,num2,num3))

#####################################
## Geographical dimension
#####################################

## Percentage of different errors
# Zero values in Lat/Long
table(df.bio$GeoNotZero)
table(df.bio$GeoNotZero)/nrow(df.bio)

# Located in oceans
table(df.bio$Overland)
table(df.bio$Overland)/nrow(df.bio)

# Located outside of the national boundary 
# (Note: exclude the records located in ocean area)
table(df.bio$CountrySame)
num.boundary <- length(which(df.bio$CountrySame == 0 & df.bio$Overland != 0))

# Centroids of administrative areas
table(df.bio$FlagCent)
table(df.bio$FlagCent)/nrow(df.bio)

# Coordinates of institutions
table(df.bio$FlagIns)
table(df.bio$GeoNotZero)/nrow(df.bio)

# Coordinates of botanic gardens
table(df.bio$FlagBGCI)
table(df.bio$FlagBGCI)/nrow(df.bio)

## Three levels of cleaning
vec1 <- rep(NA,3)
vec2 <- rep(NA,3)

# High level
subuse1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=3) & (GeoNotZero == 1) 
                  & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) 
                  & (is.na(FlagIns)) & (is.na(FlagCent)))
Name1 <- subuse1$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num1  <- get.speNo(Name1,use.spe.list)
num1
vec1[1] <- nrow(subuse1); vec2[1] <- nrow(subuse1)/nrow(df.bio)
rm(subuse1)

# Medium level
subuse2 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==2) & (GeoNotZero == 1) 
                  & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) 
                  & (is.na(FlagIns)) & (is.na(FlagCent)))
vec1[2] <- nrow(subuse2); vec2[2] <- nrow(subuse2)/nrow(df.bio)
Name2 <- subuse2$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num2  <- get.speNo(Name2,use.spe.list)
num2
rm(subuse2)

# Low level
subuse3 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==1) & (GeoNotZero == 1) 
                  & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) 
                  & (is.na(FlagIns)) & (is.na(FlagCent)))
vec1[3] <- nrow(subuse3); vec2[3] <- nrow(subuse3)/nrow(df.bio)

Name3 <- subuse3$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num3  <- get.speNo(Name3,use.spe.list)
num3
rm(subuse3)

# output relative results
df.geo <- data.frame(Num = vec1, Per = vec2, Spe = c(num1,num2,num3))



############################################
## Taxonomy and Space
##########################################
vec1 <- rep(NA,3)
vec2 <- rep(NA,3)

# High level
subuse1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=3) & (TaxonStatus %in% c('Accepted','Synonym')) 
                  &(TaxonConf == 3)  & ((TaxonType %in% tax.mat))
                  & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))
vec1[1] <- nrow(subuse1); vec2[1] <- nrow(subuse1)/nrow(df.bio)
Name1 <- subuse1$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num1  <- get.speNo(Name1,use.spe.list)
num1
rm(subuse1)

# Medium level
subuse2.1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==2) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf >= 2)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

subuse2.2 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=3) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf == 2)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

subuse2 <- rbind(subuse2.1,subuse2.2)
rm(subuse2.1);rm(subuse2.2)

vec1[2] <- nrow(subuse2); vec2[2] <- nrow(subuse2)/nrow(df.bio)
Name2 <- subuse2$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num2  <- get.speNo(Name2,use.spe.list)
num2
rm(subuse2)


# Low level
subuse3.1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==1) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf >= 1)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

subuse3.2 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=2) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf == 1)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))
subuse3 <- rbind(subuse3.1,subuse3.2)
rm(subuse3.1);rm(subuse3.2)
vec1[3] <- nrow(subuse3); vec2[3] <- nrow(subuse3)/nrow(df.bio)
rm(subuse3)
Name3 <- subuse3$scientificName %>% as.character %>% unique %>% tolower %>% attain.latin() %>% unique()
num3  <- get.speNo(Name3,use.spe.list)
num3
rm(subuse3)

df.2dim <- data.frame(Num = vec1, Per = vec2, Spe = c(num1,num2,num3))


#######################################
## code of plots
#####################################
worldmap2 <- readOGR('YourPath/shapefile/worldmap2.shp') # Natural Earth world map
ww2       <- fortify(worldmap2)

# Transform df.bio to spatial points
pro.shp<- SpatialPointsDataFrame(df.bio[which(df.bio$ExGeo == 1),c('decimalLongitude','decimalLatitude')],
                                 data = df.bio[which(df.bio$ExGeo == 1),c("scientificName",'FlagCountry','countryCode')],
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))

# Font and ggplot2 set-up
windowsFonts(JP1 = 'cambria')
mytheme <- theme_bw()+theme(plot.title = element_text(size = rel(2),hjust = 0.5),
                            axis.title.x = element_text(size=rel(2)),
                            axis.title.y = element_text(size=rel(2)),
                            plot.margin = margin(0.6,0.6,0.6,0.6, "cm"),
                            axis.ticks.length = unit(0.3,'cm') ,
                            legend.text = element_text(size=rel(1.5),lineheight = 0.8),legend.key.height = unit(1.5,"cm"),
                            legend.title = element_text(size=rel(1.5)),legend.background = element_rect(fill = 'white',colour=NA),
                            axis.text.x = element_text(size=rel(2),margin = margin(0.3,0,0.3,0,"cm"),face = 'bold'),
                            axis.text.y = element_text(size=rel(2),margin = margin(0,0.3,0,0.3,"cm"),face = 'bold'),legend.position='bottom'
                            # ,panel.grid.major = element_blank()???
                            ,panel.grid.minor = element_blank()
                            ,panel.grid.major = element_blank()
                            ,panel.border = element_rect(size=1))

###########################################################################
## Plot1: The distribution of occurrence records in integrated data set
###########################################################################

# Create 1 degree grids
r1 <- raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
vals <- 1:ncell(r1)
r2   <- setValues(r1, vals)

# Calculate the number of records in each grid
use.ex.r2 <- raster::extract(r2,pro.shp)
tb.ex     <- table(use.ex.r2)
num.value <- rep(0,360*180)
use.num   <- names(tb.ex) %>% as.integer()
for(i in 1:length(tb.ex)) {num.value[use.num[i]] <- tb.ex[i]}
r3 <- setValues(r2,num.value)


# Out put the result as spatial pixel data frame
test_spdf <- as(r3,"SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c('value','x','y')
test_df$value[which(test_df$value == 0)] <- NA
test_df$value2 <- log10(test_df$value)
test_df$Number_of_Occurrence <- test_df$value2

# Plot use ggplot2 package

g1 <- ggplot()+geom_polygon(data = ww2,aes(x = long, y = lat, group = group),fill='grey50',color = 'black',size=0.25)
g2 <- g1+geom_raster(data = test_df,aes(x=x,y=y,fill = Number_of_Occurrence),alpha=0.8,na.rm = T)+
  scale_fill_gradientn(colours = c('mediumpurple2','cyan','lawngreen','greenyellow','yellow','orange','red'),
                       na.value = 'white',limits = c(0,6),breaks = seq(0,6,1),
                       labels =format(10^c(0:6),scientific = T))
g3 <- g2+mytheme+theme(legend.key.size=unit(0.8,'cm'))+theme(legend.key.width=unit(3,'cm'))+scale_x_continuous(limit=c(-180,180),breaks=seq(-180,180,60),expand = c(0.01,0.01))+scale_y_continuous(limit=c(-90,90),breaks=seq(-90,90,45))
g4 <- g3+theme(axis.title = element_blank())
ggsave('YourPath/distributin_all.png',dpi = 500,width = 16, height = 12)

########################################################################################
## Plot: distribution of records after different levels of taxonomy and space cleaning
#######################################################################################
tax.mat <- c("Exact Match","LD Match_Type1","LD Match_Type2","Processing Exact")

## High
subuse1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=3) & (TaxonStatus %in% c('Accepted','Synonym')) 
                  &(TaxonConf == 3)  & ((TaxonType %in% tax.mat))
                  & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

## Medium
subuse2.1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==2) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf >= 2)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

subuse2.2 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=3) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf == 2)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))
subuse2 <- rbind(subuse2.1,subuse2.2)
rm(subuse2.1);rm(subuse2.2)

## Low
subuse3.1 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision==1) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf >= 1)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))

subuse3.2 <- subset(df.bio, (ExGeo == 1) & (GeoPrecision>=2) & (TaxonStatus %in% c('Accepted','Synonym'))
                    &(TaxonConf == 1)  & ((TaxonType %in% tax.mat))
                    & (GeoNotZero == 1) & !is.na(Overland) & !is.na(CountrySame) & (is.na(FlagBGCI)) & (is.na(FlagIns)) & (is.na(FlagCent)))
subuse3 <- rbind(subuse3.1,subuse3.2)
rm(subuse3.1);rm(subuse3.2)

## Transform to spatial points
# High
pro.shp.l1 <- SpatialPointsDataFrame(subuse1[which(subuse1$ExGeo == 1),c('decimalLongitude','decimalLatitude')],
                                     data = subuse1[which(subuse1$ExGeo == 1),c("scientificName",'FlagCountry','countryCode')],
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
# Medium
pro.shp.l2 <- SpatialPointsDataFrame(subuse2[which(subuse2$ExGeo == 1),c('decimalLongitude','decimalLatitude')],
                                     data = subuse2[which(subuse2$ExGeo == 1),c("scientificName",'FlagCountry','countryCode')],
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
# Low
pro.shp.l3 <- SpatialPointsDataFrame(subuse3[which(subuse3$ExGeo == 1),c('decimalLongitude','decimalLatitude')],
                                     data = subuse3[which(subuse3$ExGeo == 1),c("scientificName",'FlagCountry','countryCode')],
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
# Extraction results
use.l1  <- raster::extract(r2,pro.shp.l1)
use.l2  <- raster::extract(r2,pro.shp.l2)
use.l3  <- raster::extract(r2,pro.shp.l3)
uselist <- list(use.l1,use.l2,use.l3)

# Plots
for(j in 1:3)
{
  use.raster <- uselist[[j]]
  num.value.l1 <- rep(0,360*180)
  tb.l1        <- table(use.raster)
  use.num.l1   <- names(tb.l1) %>% as.integer()
  for(i in 1:length(tb.l1)) {num.value.l1[use.num.l1[i]] <- tb.l1[i]}
  
  r5 <- setValues(r2,num.value.l1)
  level1.spdf <- as(r5,"SpatialPixelsDataFrame")
  test_level1 <- as.data.frame(level1.spdf)
  colnames(test_level1) <- c('value','x','y')
  
  test_level1$value[which(test_level1$value == 0)] <- NA
  test_level1$`Number_of_occurrence` <- log10(test_level1$value)
 
  p1 <- ggplot()+geom_polygon(data = ww2,aes(x = long, y = lat, group = group),fill='grey50',color = 'black',size=0.25)
  p2 <- p1+geom_raster(data = test_level1,aes(x=x,y=y,fill = `Number_of_occurrence`),alpha=0.8,na.rm = T)+
    scale_fill_gradientn(colours = c('mediumpurple2','cyan','lawngreen','greenyellow','yellow','orange','red'),
                         na.value = 'white',limits = c(0,6),breaks = seq(0,6,1),
                         labels =format(10^c(0:6),scientific = T))
  
  p3 <- p2+mytheme+theme(legend.key.size=unit(0.8,'cm'))+theme(legend.key.width=unit(3,'cm'))+scale_x_continuous(limit=c(-180,180),breaks=seq(-180,180,60),expand = c(0.01,0.01))+scale_y_continuous(limit=c(-90,90),breaks=seq(-90,90,45))
  p4 <- p3+theme(axis.title = element_blank())
  
  ggsave(paste0('YourPath/subclean_',j,'.png'),dpi = 300,width = 16, height = 12)
}


#######################################################
## The distribution of records that can be corrected
########################################################

# Extract the records
no.sub <- which(df.bio$Adj.Coord == 1)
subuse <- df.bio[no.sub,]

## location.adj: 0, not used; 1, adj used
location.adj <- rep(0,nrow(df.bio))
location.adj[no.sub] <- 1
coord.adj <- data.frame(array(NA,dim = c(length(no.sub),2)))
names(coord.adj) <- c('lon','lat')
num.adj <- which(df.bio$Adj.Coord == 1)

# Transform to data.frame
df.adj <- data.frame(num.row = num.adj, Type = df.bio$Adj.TransTypr[num.adj], 
                     Lon = df.bio$decimalLongitude[num.adj], 
                     Lat = df.bio$decimalLatitude[num.adj],
                     AdjLon = NA, AdjLat = NA)


# Spatial points: before correction
pro.before <- SpatialPointsDataFrame(df.adj[,c('Lon','Lat')],
                                     data = df.adj[,c("num.row",'Type')],
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))
# Spatial points: after correction
pro.after  <- SpatialPointsDataFrame(df.adj[,c('AdjLon','AdjLat')],
                                     data = df.adj[,c("num.row",'Type')],
                                     proj4string = CRS("+proj=longlat +datum=WGS84"))

# Before
p1 <- ggplot()+geom_polygon(data = ww2,aes(x = long, y = lat, group = group),fill='grey95',color = 'black',size=0.25)
p2 <- p1+geom_point(data = df.adj,aes(x = Lon, y = Lat),color = cb_palette[1],alpha = 0.9,size = 1.1)
p3 <- p2+mytheme+theme(legend.key.size=unit(0.8,'cm'))+theme(legend.key.width=unit(3,'cm'))+scale_x_continuous(limit=c(-180,180),breaks=seq(-180,180,60),expand = c(0.01,0.01))+scale_y_continuous(limit=c(-90,90),breaks=seq(-90,90,45))
p4 <- p3+theme(axis.title = element_blank())
ggsave('YourPath/beforpts.png',dpi = 500,width = 16, height = 10)


# After
cb_palette <- c(rainbow(8))
p1 <- ggplot()+geom_polygon(data = ww2,aes(x = long, y = lat, group = group),fill='grey95',color = 'black',size=0.25)
p2 <- p1+geom_point(data = df.adj,aes(x = AdjLon, y = AdjLat, color = Type))
p3 <- p2+mytheme+theme(legend.key.size=unit(0.8,'cm'))+theme(legend.key.width=unit(3,'cm'))+scale_x_continuous(limit=c(-180,180),breaks=seq(-180,180,60),expand = c(0.01,0.01))+scale_y_continuous(limit=c(-90,90),breaks=seq(-90,90,45))
p4 <- p3+theme(axis.title = element_blank())
p5 <- p4+scale_colour_manual(values=sample(cb_palette[-1]))+theme(legend.position = 'none')
ggsave('YourPath/afterpts.png',dpi = 500,width = 16, height = 10)

p6 <- p5+theme(legend.position = 'bottom')
ggsave('YourPath/afterpts_legend.png',dpi = 500,width = 16, height = 10)

