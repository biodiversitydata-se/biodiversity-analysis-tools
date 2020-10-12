####################################################
## Script 3. Data cleaning of the merged dataset
## Author: Jing Jin, Jun Yang*
## Manuscript: BDcleaner: a workflow for cleaning taxonomic and geographic errors
##             in occurrence data archived in biodiversity databases
####################################################


################################
# load merged data set: df.bio
################################
load("YourPath/dfbio.RData")

# load packages
library(tidyverse)
library(data.table)
library(magrittr)

# data format trasformation (the coordinates' format of some data base are "string")
df.bio$decimalLatitude  <- df.bio$decimalLatitude %>% as.numeric
df.bio$decimalLongitude <- df.bio$decimalLongitude %>% as.numeric


##########################################################################################
## Part 1. Taxon name cleaning
## Dimension: Taxonomy
## Important input variables:
##   (1) taxonbind (data.frame): the TPL name list downloaded fron the Plant List home page
## Important output variables:
##   (1) judge.taxon (data.frame): the TPL match result of unique scientif names in dataset
##   (2) TaxonType (vector): the TPL match names' type (based on edit distence)
##   (3) TaxonStatus (vector): the TPL match names' status (accepted,synonym,misapplied,unresolved)
##   (4) TaxonConf (vector): the TPL match names' confidence level (1-low,2-medium,3-high)
## Name status description in manuscript and scripts below: 
##   (1) Exact Match: Exact match
##   (2) Match after removing punctuations: Processing Exact
##   (3) Match after removing spelling errors:LD Match_Type1,LD Match_Type2
##   (4) No match: LD DisMatch_Type1, LD DisMatch_Type2, LD DisMatch_Type3
##   (5) NA value: no genus name match result
## If you only need consider binomial nomenclature, you can turn on the simplematch code
##########################################################################################

# load TPL and other taxon names information
load("YourPath/DataUse.RData") # From script 1. 

# Variable preparation
name.vec <- c('TPL_Match','TPL_Status','TPL_AccID','TPL_Conf','TPL_Line')
name.vec.j <- c("bindname","Taxonomic status in TPL","Accepted ID","conf","Line")
LD.name <- c('LD_Spe','LD_Aut','RLD_Spe','RLD_Aut')

# unique taxon names in dataframe
Name.occ <- unique(df.bio$scientificName) %>% as.character
Name.occ2 <- Name.occ[which(!str_detect(Name.occ,' [×xX] '))] # remove hybrids species


# Create a data frame to save status of different taxon names
judge.taxon <- data.frame(array(NA,dim = c(length(Name.occ2),17)))
names(judge.taxon) <- c('SciName_GBIF','Sciname',"Genus.name","Spe.name","Author.name","SpeType",
                        'MatchType','LD_Spe','LD_Aut','RLD_Spe','RLD_Aut',"ExAutInf",
                        'TPL_Match','TPL_Status','TPL_AccID','TPL_Conf','TPL_Line')
judge.taxon$SciName_GBIF <- Name.occ2

# taxon name processing
judge.use <- ProBlank(attain.sciname(judge.taxon$SciName_GBIF))
judge.taxon[,c("Genus.name","Spe.name","Author.name","SpeType")] <- attain.spname(judge.use)
judge.taxon$Sciname  <- judge.use
# Author information available?
judge.taxon$ExAutInf <- 1
judge.taxon$ExAutInf[which(judge.taxon$Author.name == "")] <- 0
judge.taxon$MatchSim <- NA

# Match taxon name with TPL
flag.next <- NULL
NA.spe    <- NULL
for(k in 1:nrow(judge.taxon))
{
  tpl.judge   <- subset(taxonbind,Bindname == judge.use[k])
  nrow.tj     <- nrow(tpl.judge)
  
  if(nrow.tj == 1) # Match the names in TPL exactly
  {
    judge.taxon[k,name.vec]  <- tpl.judge[,c("Bindname","Taxonomic status in TPL","Accepted ID","conf","Line")]
    judge.taxon$MatchType[k] <- 'Exact Match'
  }
  
  if(nrow.tj > 1) 
  {
    ru <- which(tpl.judge$conf == max(tpl.judge$conf))[1]
    judge.taxon[k,name.vec] <- tpl.judge[ru,c("Bindname","Taxonomic status in TPL","Accepted ID","conf","Line")]
    judge.taxon$MatchType[k] <- 'Exact Match'
  }
  
  if(nrow.tj == 0) # Other situation
  {
    # Genus name
    genus.use <- judge.taxon$Genus.name[k]
    # Specific epithet
    spe.use   <- tolower(judge.taxon$Spe.name[k])
    # Author information
    aut.use   <- tolower(judge.taxon$Author.name[k])
    
    if(is.na(spe.use)) 
    {
      NA.spe <- c(NA.spe,k)
      next 
    }
    
    # Ensure the name of Genus match exactly
    sub.pro1  <- subset(taxonbind,Genus == genus.use)
    if(nrow(sub.pro1) == 0)
    {
      flag.next <- c(flag.next,k)
      next
    }
    
    # Edit distance calculation of specific epithet
    LD.ch1    <- sub.pro1$bind_speepi
    LD.re1.ort <- sort(adist(spe.use, LD.ch1),index.return = T)
    LD.re1.rel <- LD.re1.ort[[1]]/str_length(spe.use)
    
    ## Best Matching
    judge.0 <- which(LD.re1.ort[[1]] == 0)
    if(length(judge.0) != 0)
    {
      # Edit distance calculation of author information
      sub.pro2 <- sub.pro1[LD.re1.ort[[2]][judge.0],]
      LD.ch2    <- sub.pro2$AuthorInf
      LD.re2.ort <- sort(adist(aut.use, LD.ch2),index.return = T)
      LD.re2.rel <- LD.re2.ort[[1]]/str_length(aut.use)
      judge.aut.zero <- which(LD.re2.ort[[1]] == 0)
      
      if(length(judge.aut.zero) != 0)
      {
        subpro2.line <- LD.re2.ort[[2]][judge.aut.zero[1]]
        use.data <- sub.pro2[subpro2.line,]
        judge.taxon[k,LD.name]   <- rep(0,4)
        judge.taxon[k,name.vec]  <- use.data[,c("Bindname","Taxonomic status in TPL","Accepted ID","conf","Line")]
        ## This classification is the scientif name can match the TPL name after removing special punctations 
        judge.taxon$MatchType[k] <- 'Processing Exact'
      } else
      {
        subpro2.line <- LD.re2.ort[[2]][1]
        use.data     <- sub.pro2[subpro2.line,]
        judge.taxon[k,LD.name]   <- c(0,LD.re2.ort[[1]][1],0,LD.re2.rel[1])
        judge.taxon[k,name.vec]  <- use.data[,c("Bindname","Taxonomic status in TPL","Accepted ID","conf","Line")]
        if(LD.re2.ort[[1]][1] < 4 & LD.re2.rel[1] < 0.3)
        {
          ## This classification is the scientif name can match the TPL name after removing spelling errors
          judge.taxon$MatchType[k] <- 'LD Match_Type1'
        } else
        {
          ## This classification is the scientif name cannot match the TPL name
          judge.taxon$MatchType[k] <- 'LD DisMatch_Type1'
          
          # if(aut.use == "") judge.taxon$MatchSim[k] <- 'SimpleNameMatch'
        }
      }
    } else
    {
      
      #Step1: judge the spe information distance
      spe.min    <- which(LD.re1.ort[[1]] == min(LD.re1.ort[[1]]))
      sub.pro2   <- sub.pro1[LD.re1.ort[[2]][spe.min],]
      LD.ch2     <- sub.pro2$AuthorInf
      LD.re2.ort <- sort(adist(aut.use, LD.ch2),index.return = T)
      LD.re2.rel <- LD.re2.ort[[1]]/str_length(aut.use)
      subpro2.line <- LD.re2.ort[[2]][1]
      use.data     <- sub.pro2[subpro2.line,]
      judge.taxon[k,LD.name]   <- c(LD.re1.ort[[1]][1],LD.re2.ort[[1]][1],LD.re1.rel[1],LD.re2.rel[1])
      judge.taxon[k,name.vec]  <- use.data[,c("Bindname","Taxonomic status in TPL","Accepted ID","conf","Line")]
      if(LD.re1.ort[[1]][1] < 4 & LD.re1.rel[1] < 0.3)
      {
        
        # if(aut.use == "") judge.taxon$MatchSim[k] <- 'SimpleNameMatch_2'
        
        if(LD.re2.ort[[1]][1] < 4 & (LD.re2.rel[1] < 0.3 | (aut.use == "" & use.data$AuthorInf == "")))
        {
          ## This classification is the scientif name can match the TPL name after removing spelling errors
          judge.taxon$MatchType[k] <- 'LD Match_Type2'
        } else
        {
          ## This classification is the scientif name cannot match the TPL name
          judge.taxon$MatchType[k] <- 'LD DisMatch_Type2'
        }
      } else
      {
        ## This classification is the scientif name cannot match the TPL name
        judge.taxon$MatchType[k] <- 'LD DisMatch_Type3'
      }
    }
  }
}

# save.image('Your path and file name') 


# Match the names status in judge.taxon with the df.bio dataframe
Name.occ.full <- df.bio$scientificName %>% as.character
TaxonStatus  <- rep(NA,nrow(df.bio)) # Status in TPL
TaxonConf    <- rep(NA,nrow(df.bio)) # Confidence level in TPL
TaxonType    <- rep(NA,nrow(df.bio)) # Match type in TPL
# TaxonType2   <- rep(NA,nrow(df.bio)) # Not consider author information

int.ra <- data.frame(sp.mt(1,nrow(df.bio),3000000))

# Match the name status with occurrence data in df.bio
for(j in 1:nrow(int.ra))
{
  Name.occ.use    <- Name.occ.full[int.ra$start[j]:int.ra$end[j]]
  TaxonStatus.sub <- rep(NA,length(Name.occ.use))
  TaxonType.sub   <- rep(NA,length(Name.occ.use))
  # TaxonType2.sub  <- rep(NA,length(Name.occ.use))
  TaxonConf.sub   <- rep(NA,length(Name.occ.use))
  judge.taxon.sub <- judge.taxon[which(judge.taxon$SciName_GBIF %in% Name.occ.use),]
  
  for(i in 1:nrow(judge.taxon.sub))
  {
    row.occ <- which(Name.occ.use == judge.taxon.sub$SciName_GBIF[i])
    TaxonStatus.sub[row.occ] <- judge.taxon.sub$TPL_Status[i]
    TaxonType.sub[row.occ]   <- judge.taxon.sub$MatchType[i]
    # TaxonType2.sub[row.occ]  <- judge.taxon.sub$MatchSim[i]
    TaxonConf.sub[row.occ]   <- judge.taxon.sub$TPL_Conf[i]
  }
  
  TaxonStatus[int.ra$start[j]:int.ra$end[j]] <- TaxonStatus.sub
  TaxonType[int.ra$start[j]:int.ra$end[j]]   <- TaxonType.sub
  # TaxonType2[int.ra$start[j]:int.ra$end[j]]  <- TaxonType2.sub
  TaxonConf[int.ra$start[j]:int.ra$end[j]]   <- TaxonConf.sub
}

# Mark correlated variables
df.bio$TaxonType   <- TaxonType
# df.bio$TaxonType2  <- TaxonType2
df.bio$TaxonStatus <- TaxonStatus
df.bio$TaxonConf   <- TaxonConf

# save.image('YourFile')
rm(list = setdiff(ls(),c('df.bio','add.var','taxonbind','judge.taxon','flag.next')))


###########################################################################
## Part 2. Geographical cleaning
## Dimension: Space
## Important output variables
##  (1)  ExGeo: the existence of geo-coordinates
##  (2)  GeoPrecision: the number of decimal digits of coordinates
##  (3)  TaxonType: match type of taxon name (with TPL)
##  (4)  TaxonStatus: name status of taxon names (TPL)
##  (5)  TxonComf: confidence level of taxon names (TPL)
##  (6)  Adj.Coord: whether the error coordinate can be corrected
##  (7)  Adj.TransTypr: type of coordinate transformation (correct)
##  (8)  CountrySame: whether the coordinates match the political boundary?
##  (9)  Overland: whether the coordinated located in terrsital area
##  (10) GeoNotZero: Lat/Lon == 0?
##  (11) FlagCent: centroids of administrative area
##  (12) FlagIns: location of biodiversity institution
##  (13) FlagBGCI: location of botanical garden
##  (14) over.urban: whether located in urban area
##  (15) overUrbanType: urban cluster of urban core or rural area
############################################################################

# Add new variable to df.bio dataframe
add.var <- c("ExGeo","GeoPrecision")
df.bio[,add.var] <- rep(NA,nrow(df.bio))

## ExGeo: existence of geo-coordinates
judge.geo         <- is.na(df.bio$decimalLongitude)+is.na(df.bio$decimalLatitude) # True: 0; False: 1
no.geo            <- which(judge.geo==0)
le.geo            <- length(no.geo)
df.bio$ExGeo[no.geo]  <- 1
df.bio$ExGeo[-no.geo] <- 0

# GeoPrecision Use: number of decimal digits of coordinates
clat <- df.bio$decimalLatitude[no.geo] %>% as.character %>% str_split(pattern = '[.]')
clon <- df.bio$decimalLongitude[no.geo] %>% as.character %>% str_split(pattern = '[.]')

geop.use <- rep(0,length = length(no.geo))
for(p in 1:length(no.geo))
{
  if(is.na(clat[[p]][2]) | is.na(clon[[p]][2]))
  {next} else
  {
    geop.use[p] <- max(str_length(clat[[p]][2]),str_length(clon[[p]][2]))
  }
}
df.bio$GeoPrecision[no.geo] <- geop.use

# remove useless variables
rm(list = setdiff(ls(),c('df.bio','add.var')))


##################################
## Geographical error detection
####################################
library(sp)
library(rgdal)
library(rworldxtra)
library(countrycode)
library(rworldmap)
library(biogeo)
library(maptools)
library(geonames)
library(data.table)

# Read world map data

data(wrld_simpl) # world map in maptools R package (coarse)
worldmap2 <- readOGR('YourPath/worldmap2.shp') # Natural earth wold map (medium)
worldmap  <- readOGR('YourPath/gadm28_adm0.shp') # GADM world map (fine)

# remove useless variables
rm(list = setdiff(ls(),c('df.bio','add.var','taxonbind',
                         'judge.taxon','flag.next','worldmap','worldmap2')))

# Add filter variable names
add.name <- c('GeoNotZero','OverLand','CountrySame','Adj.Overland','Adj.Type',
              'Adj.Correct','AdjLat','AdjLon')

############################
# GeoNotZero: lat/lon = 0?
############################
GeoNotZero  <- Coord.Notzero(df.bio$decimalLatitude,df.bio$decimalLongitude)
df.bio$GeoNotZero <- GeoNotZero
rm(GeoNotZero)

country.inf <- wrld_simpl@data
country.inf$ISO2 <- as.character(country.inf$ISO2)

# Some special lower case
df.bio$countryCode <- toupper(df.bio$countryCode %>% as.character)

# "" not NA if countryCode is NA

###################################
## GeoFilter 2: On the Land?
## GeoFilter 3: Country boundary
## We used three world maps
###################################
Line.use <- which(df.bio$ExGeo == 1)
df.use   <- df.bio[Line.use,]
L.line   <- length(Line.use)

pro.shp <- SpatialPointsDataFrame(df.use[,c('decimalLongitude','decimalLatitude')],
                                  data = df.use[,c("scientificNamDe",'FlagCountry','countryCode')],
                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
rm(df.use)

overlay.land <- rep(0,L.line)
country.same <- rep('Simple',L.line)

# Simple world map (in maptools package)
pro.shp@proj4string <- wrld_simpl@proj4string

# split data to ensure enough computer memory
use.ci <- sp.mt(1,nrow(pro.shp),3000000)

# Initializing the loop process
o1 <- over(pro.shp[1,],wrld_simpl)
over.shp1 <- data.frame(array(NA,dim = c(nrow(pro.shp),ncol(o1))))
names(over.shp1) <- names(o1)

## Overlay analysis
for(i in 1:nrow(use.ci))
{
  over.pro <- pro.shp[use.ci$start[i]:use.ci$end[i],]
  overuse  <- over(over.pro,wrld_simpl)
  
  overuse$FIPS      <- as.character(overuse$FIPS)
  overuse$ISO2      <- as.character(overuse$ISO2)
  overuse$ISO3      <- as.character(overuse$ISO3)
  overuse$UN        <- as.character(overuse$UN)
  overuse$NAME      <- as.character(overuse$NAME)
  overuse$REGION    <- as.character(overuse$REGION)
  overuse$SUBREGION <- as.character(overuse$SUBREGION)
  
  over.shp1[use.ci$start[i]:use.ci$end[i],] <- overuse
  
  rm(overuse); rm(over.pro)
  print(i)
}

# Flag "simple" when the points located on simple world map
land.notna               <- which(!is.na(over.shp1$ISO2))
overlay.land[land.notna] <- 'Simple'


## Country boundary detection
dif.ju    <- as.character(over.shp1$ISO2) == pro.shp@data$countryCode
NA.number <- which(is.na(dif.ju))
same.diff <- sort(union(which(!dif.ju),NA.number))
country.same[same.diff] <- 0
#country.same[NA.number] <- NA

# Natural world map
pro.shp2  <- pro.shp[same.diff,]
pro.shp2@proj4string <- worldmap2@proj4string
overuse.shp2         <- over(pro.shp2,worldmap2)

overland2 <- which(!is.na(overuse.shp2$ISO_A2))
overlay.land[same.diff[overland2]] <- str_replace_all(overlay.land[same.diff[overland2]],'0','NE')


dif.ju2    <- as.character(overuse.shp2$ISO_A2) == pro.shp2@data$countryCode
NA.number2 <- which(is.na(dif.ju2))
same.diff2 <- sort(union(which(!dif.ju2),NA.number2))
country.same[same.diff[which(dif.ju2)]] <- str_replace_all(country.same[same.diff[which(dif.ju2)]],'0','NE')

diff.use <- same.diff[same.diff2]

if(length(same.diff2) != 0)
{
  # GADM worldmap overlay
  pro.shp3 <- pro.shp2[same.diff2,]
  pro.shp3@proj4string <- worldmap@proj4string
  
  # overuse.shp3        <- over(pro.shp3,worldmap)
  
  use.ci3 <- sp.mt(1,nrow(pro.shp3),20000)
  
  o3 <- over(pro.shp3[1,],worldmap)
  over.shp3 <- data.frame(array(NA,dim = c(nrow(pro.shp3),ncol(o3))))
  names(over.shp3) <- names(o3)
  
  
  for(i in 1:nrow(use.ci3))
  {
    over.pro <- pro.shp3[use.ci3$start[i]:use.ci3$end[i],]
    overuse  <- over(over.pro,worldmap)
    
    overuse$NAME_ENGLI <- as.character(overuse$NAME_ENGLI)
    overuse$ISO2       <- as.character(overuse$ISO2)
    overuse$ISO        <- as.character(overuse$ISO)
    
    
    over.shp3[use.ci3$start[i]:use.ci3$end[i],] <- overuse
    
    rm(overuse); rm(over.pro)
    print(i)
    
  }
  
  # GADM world map results
  overland3 <- which(!is.na(over.shp3$ISO2))
  overlay.land[same.diff[same.diff2[overland3]]] <- str_replace_all(overlay.land[same.diff[same.diff2[overland3]]],'0','GADM')
  dif.ju3    <- as.character(over.shp3$ISO2) == pro.shp3@data$countryCode
  NA.number3 <- which(is.na(dif.ju3))
  same.diff3 <- sort(union(which(!dif.ju3),NA.number3))
  country.same[same.diff[same.diff2[which(dif.ju3)]]] <- str_replace_all(country.same[same.diff[same.diff2[which(dif.ju3)]]],'0','GADM')
  ####################################
  diff.use <- same.diff[same.diff2[same.diff3]]
  # save.image('file')
}

################################################################
## Correct coordinates
################################################################
adj.type     <- rep(NA,L.line) # transform types of coordinates
adj.lat      <- rep(NA,L.line) # latitude after correction
adj.lon      <- rep(NA,L.line) # longitude after correction
adj.correct  <- rep(NA,L.line) # Whether the coordinate need correct
adj.over     <- rep(NA,L.line) # After the correction process, where the coordinate located in?


## 7 types
coord.save   <- data.frame(x = rep(NA,length(diff.use)*7),y = rep(NA,length(diff.use)*7))
country.save <- rep(NA,length(diff.use)*7)
diff.save    <- rep(NA,length(diff.use)*7)

# country.same[diff.use[which(df.use$countryCode[diff.use] == "")]] <- NA
diff.use <- setdiff(diff.use, diff.use[which(df.use$countryCode[diff.use] == "")])

x <- df.use$decimalLongitude[diff.use]
y <- df.use$decimalLatitude[diff.use]
use.iso2 <- df.use$countryCode[diff.use] %>% as.character

# Transformation process and tests (Simple world map)
over.list1  <- Coord8.Out(coordx = x, coordy = y, worldmap_use = wrld_simpl)
match.trans <- Coord.match.trans(ISO.rec = use.iso2, match.overlist = over.list1,
                                 NameVec = 'ISO2')

# If the correction process were successed, we flag some variables
adj.correct[diff.use[which(!is.na(match.trans))]] <- 1
adj.over[diff.use[which(!is.na(match.trans))]] <- 'Simple'
adj.type[diff.use[which(!is.na(match.trans))]] <- match.trans[which(!is.na(match.trans))]



# Transformation process and tests (Natural Earth world map)
diff.use2 <- diff.use[which(is.na(match.trans))]
use.iso2.2 <- use.iso2[which(is.na(match.trans))]

x2 <- df.use$decimalLongitude[diff.use2]
y2 <- df.use$decimalLatitude[diff.use2]

over.list2 <- Coord8.Out(coordx = x2, coordy = y2, 
                         worldmap_use = worldmap2)
match.trans2 <- Coord.match.trans(ISO.rec = use.iso2.2, 
                                  match.overlist = over.list2,
                                  NameVec = 'ISO_A2')

table(match.trans2)
adj.correct[diff.use2[which(!is.na(match.trans2))]] <- 1
adj.over[diff.use2[which(!is.na(match.trans2))]] <- 'NE'
adj.type[diff.use2[which(!is.na(match.trans2))]] <- match.trans2[which(!is.na(match.trans2))]

# Transformation process and tests (GADM world map)
diff.use3 <- diff.use2[which(is.na(match.trans2))]
use.iso3 <- use.iso2.2[which(is.na(match.trans2))]

x3 <- df.use$decimalLongitude[diff.use3]
y3 <- df.use$decimalLatitude[diff.use3]

over.list3 <- Coord8.Out(coordx = x3, coordy = y3, worldmap_use = worldmap)
match.trans2 <- Coord.match.trans(ISO.rec = use.iso3, match.overlist = over.list3,
                                  NameVec = 'ISO')

# Flag relative variables in df.bio data.frame
df.bio$CountrySame[Line.use]   <- country.same
df.bio$Overland[Line.use]      <- overlay.land
df.bio$Adj.Coord[Line.use]     <- adj.correct
df.bio$Adj.overMap[Line.use]   <- adj.over
df.bio$Adj.TransTypr[Line.use] <- adj.type

# remove useless variables
rm(list = setdiff(ls(),c('df.bio','add.var','taxonbind',
                         'judge.taxon','flag.next','add.name')))



############################################
## Potential error test in space dimension
############################################
library(plyr)

# Read distribution data
tb.coord <- count(df.bio[,c('decimalLatitude','decimalLongitude')])
tb.coord$sumcoord <- tb.coord$decimalLatitude+tb.coord$decimalLongitude

# Institution data processing
# Data source: Zizka et al. (2019)
use.ins <- fread('L:/dropbox/Dropbox/Dropbox/Dropbox/data/institutions.csv')
loop.i  <- setdiff(1:nrow(tb.coord),which(is.na(tb.coord$sumcoord)))
pro1    <- which((tb.coord$decimalLatitude %in% use.ins$Lat) & 
                   (tb.coord$decimalLongitude %in% use.ins$Lon))

tb.coord$flag <- NA

for(i in pro1)
{
  sub1  <- subset(use.ins, Lat == tb.coord$decimalLatitude[i])
  j.num <- which(sub1$Lon == tb.coord$decimalLongitude[i])[1]
  
  if(length(j.num) != 0) {tb.coord$flag[i] <- sub1$Name[j.num]}
}

table(tb.coord$flag)

# Centroid data processing (GeoUse)
use.cent <- fread('YourPath/GeoNames_centroid.csv')
pro2     <- which((tb.coord$decimalLatitude %in% use.cent$lat) & 
                    (tb.coord$decimalLongitude %in% use.cent$long))
tb.coord$flagCent <- NA

for(i in pro2)
{
  sub1  <- subset(use.cent, lat == tb.coord$decimalLatitude[i])
  j.num <- which(sub1$long == tb.coord$decimalLongitude[i])[1]
  
  if(length(j.num) != 0) {tb.coord$flagCent[i] <- sub1$name[j.num]}
}

table(tb.coord$flagCent)

# Country data processing (MEE use)
library(CoordinateCleaner)
coord.inf <- countryref
loop.i  <- setdiff(1:nrow(tb.coord),which(is.na(tb.coord$sumcoord)))
pro3.1  <- which((tb.coord$decimalLatitude %in% coord.inf$centroid.lat) & 
                   (tb.coord$decimalLongitude %in% coord.inf$centroid.lon))
pro3.2  <- which((tb.coord$decimalLatitude %in% coord.inf$capital.lat) & 
                   (tb.coord$decimalLongitude %in% coord.inf$capital.lon))


tb.coord$flagMeeRef1 <- NA
tb.coord$flagMeeRef2 <- NA

for(i in pro3.1)
{
  sub1  <- subset(coord.inf, centroid.lat == tb.coord$decimalLatitude[i])
  j.num <- which(sub1$centroid.lon == tb.coord$decimalLongitude[i])[1]
  
  if(length(j.num) != 0) {tb.coord$flagMeeRef1[i] <- sub1$name[j.num];print(i)}
}


for(i in pro3.2[-1])
{
  sub1  <- subset(coord.inf, capital.lat == tb.coord$decimalLatitude[i])
  j.num <- which(sub1$capital.lon == tb.coord$decimalLongitude[i])[1]
  
  if(length(j.num) != 0) {tb.coord$flagMeeRef2[i] <- paste0(sub1$name[j.num],'_cap');print(i)}
}
table(tb.coord$flag)



# Institution and centroid data processing (to vector)
subcoord  <- subset(df.bio, select = c('scientificName','decimalLatitude',
                                       'decimalLongitude','ExGeo'))
flag.ins  <- rep(NA,nrow(df.bio))
flag.cent <- rep(NA,nrow(df.bio))

process1  <- which(!is.na(tb.coord$flag))
for(i in process1)
{
  rep.no <- which(subcoord$decimalLatitude == tb.coord$decimalLatitude[i] & 
                    subcoord$decimalLongitude == tb.coord$decimalLongitude[i])
  flag.ins[rep.no] <- tb.coord$flag[i]
  rm(rep.no)
}

process2 <- which(!is.na(tb.coord$flagCent))
for(i in process2)
{
  rep.no2 <- which(subcoord$decimalLatitude == tb.coord$decimalLatitude[i] & 
                     subcoord$decimalLongitude == tb.coord$decimalLongitude[i])
  flag.cent[rep.no2] <- tb.coord$flagCent[i]
  rm(rep.no2)
}

flag.ref1 <- rep(NA,nrow(df.bio))
process3.1 <- which(!is.na(tb.coord$flagMeeRef1))
for(i in process3.1)
{
  rep.no3.1 <- which(subcoord$decimalLatitude == tb.coord$decimalLatitude[i] & 
                     subcoord$decimalLongitude == tb.coord$decimalLongitude[i])
  flag.ref1[rep.no3.1] <- tb.coord$flagMeeRef1[i]
  rm(rep.no3.1)
}

flag.ref2 <- rep(NA,nrow(df.bio))
process3.2 <- which(!is.na(tb.coord$flagMeeRef2))
for(i in process3.2)
{
  rep.no3.2 <- which(subcoord$decimalLatitude == tb.coord$decimalLatitude[i] & 
                     subcoord$decimalLongitude == tb.coord$decimalLongitude[i])
  flag.ref2[rep.no3.2] <- tb.coord$flagMeeRef2[i]
  rm(rep.noe.2)
}

# Save results
table(flag.cent)
table(flag.ins)
table(flag.ref1)
table(flag.ref2)

df.bio$FlagCent  <- flag.cent # centroid
df.bio$FlagIns   <- flag.ins  # institution
df.bio$FlagRef1  <- flag.ref1 # adjusted coordinates
df.bio$FlagRef2  <- flag.ref2 # adjusted coordinates


#######################################
## Botanic Garden Filter
#######################################
library(data.table)
library(stringr)

garden.inf <- read.csv('YourPath/GardenInformation.csv', encongding='UTF-8')
NameGarden <- tolower(str_trim(as.character(garden.inf$`Institution Name`)))

## flag variable
flag.garden         <- rep(NA,nrow(df.bio))

# extract coordinates
garden.coord            <- garden.inf[,c("Latitude","Longitude")]
row.names(garden.coord) <- paste0('BGCI',1:nrow(garden.inf))
garden.uni.coord        <- unique(garden.coord)

na.row <- NULL
for(i in 1:nrow(garden.uni.coord))
{
  if( is.na(garden.uni.coord$Latitude[i]) | is.na(garden.uni.coord$Longitude[i]))
  {
    na.row <- c(na.row,i)
  }
}

garden.uni.coord <- garden.uni.coord[-na.row,]

# Adjusted coordinates results
num.adj <- which(df.bio$Adj.Coord == 1)
df.adj <- data.frame(num.row = num.adj, Type = df.bio$Adj.TransTypr[num.adj], 
                     Lon = df.bio$decimalLongitude[num.adj], 
                     Lat = df.bio$decimalLatitude[num.adj],
                     AdjLon = NA, AdjLat = NA)

for(i in 1:nrow(df.adj))
{
  coord.adj <- Coord8.single(x = df.adj$Lon[i], y = df.adj$Lat[i],type = df.adj$Type[i])
  df.adj[i,5:6] <- coord.adj
}

## Filter: coordinate, same with botanic garden in BGCI database
for(i in 1:nrow(garden.uni.coord))
{
  lati      <- garden.uni.coord$Latitude[i]; loni <- garden.uni.coord$Longitude[i]
  flagname  <- row.names(garden.uni.coord)[i]

  flag1.geo <- which(df.bio$decimalLatitude == lati & df.bio$decimalLongitude == loni & is.na(df.bio$Adj.Coord))
  flag2.geo <- df.adj$num.row[which(df.adj$AdjLat == lati & df.adj$AdjLon == loni)]
  flag.geo  <- c(flag1.geo,flag2.geo)
  
  if(length(flag.geo) != 0)
  {
    flag.garden[flag.geo] <- flagname
    print(i)
  } else {next}
  
}

table(flag.garden)

######
df.bio$FlagBGCI <- flag.garden
rm(list = setdiff(ls(),c('df.bio','add.var','taxonbind',
                         'judge.taxon','flag.next','add.name')))
# save.image('YourPath')


################################################################################
## Urban Flag
## 
## Note: In this part, we highly recommend you to use high performence 
##       computing service if your data volume is large. Therefore, we
##       gave the code to conduct this analysis using HPC and you can easily
##       transformed it if you only want to run this script on your PC.
#############################################################################

## High performance computing
# load packages
library(sp)
library(rgdal)
library(rworldxtra)
library(countrycode)
library(rworldmap)
library(biogeo)
library(maptools)
library(geonames)

# read GHS settlement grid data
# Code 1: Urban Center; Code 2: Urban cluster; Code 3: Rural
urban2    <- readOGR('YourPath/urban2.shp')
# Transform the projection into WGS84
urban2dis <- spTransform(urban2,CRS("+proj=longlat +datum=WGS84"))

# Creat a folder to save scripts for HPC
dir.create('YourPath/OverUrban')
setwd('YourPath/OverUrban')

# Select part of the important variables to save a smaller csv file
# For HPC use
sub.bio <- subset(df.bio,select = c("scientificName", "decimalLatitude", 
                                    "decimalLongitude","ExGeo"))
fwrite(sub.bio,'bind_sub.csv',row.names = F)

# Preparation
vec.Urban <- rep(NA,length= nrow(df.bio))
no.coord  <- which(df.bio$ExGeo == 1)
seq.ci    <- seq(1,length(vec.Urban),1000)
mt.ci     <- array(NA,dim=c(length(seq.ci),2))

for(i in 1:length(seq.ci))
{
  mt.ci[i,1] <- seq.ci[i]
  mt.ci[i,2] <- seq.ci[i]+999
  if(i == length(seq.ci)) { mt.ci[i,2] <- length(vec.Urban)}
}

## High performance computing
loop.mt  <- array(NA,dim=c(15,2))
loop.seq <- seq(1,nrow(mt.ci),2000)

for(i in 1:15)
{
  loop.mt[i,1] <- loop.seq[i]
  loop.mt[i,2] <- loop.seq[i]+1999
  if(i == 15) loop.mt[i,2] <- nrow(mt.ci)
}

start <- 1

# Creat different folders: we use 15 folders
for(i in 1:15)
{
  dir.create(paste0('sub',i))
  file.create(paste0('./sub',i,'/VecUrban',i,'.txt'))
  
}

# IMPORTANT: FOR hpc
save(mt.ci,urban2dis,vec.Urban,file='YourPath/useworkspace.RData')

# Creat different scripts files
for(m in 1:15)
{
  file.create(paste0('overUrban',m,'.R'))
  
  cat('library(sp);library(maptools);library(data.table)','\n',
      
      'load(\'./useworkspace.RData\')','\n',
      paste0('for(i in ',loop.mt[m,1],':',loop.mt[m,2],')'),'\n',
      '
      {
      ci.use <- mt.ci[i,1]:mt.ci[i,2]
      
      df.occ <- fread("./bind_sub.csv", 
      encoding = "UTF-8",skip=mt.ci[i,1],nrows = 1000)
      
      names(df.occ) <- c("scientificName", "decimalLatitude", "decimalLongitude", "ExGeo")
      have.coord    <- which(df.occ$ExGeo == 1)
      
      if(length(have.coord) == 0) next
      
      
      over.points <- df.occ[have.coord,]
      
      pro.shp <- SpatialPointsDataFrame(over.points[,c("decimalLongitude","decimalLatitude")],
      data = over.points[,c("scientificName")],
      proj4string = CRS("+proj=longlat +datum=WGS84"))
      pro.shp@proj4string <- urban2dis@proj4string
      
      overuse2 <- over(pro.shp,urban2dis)
      overlay2dis <- rep(0,length(have.coord))
      
      land.notna2  <- which(!is.na(overuse2$OBJECTID))
      overlay2dis[land.notna2] <- overuse2$OBJECTID[land.notna2]
      
      flag.vec <- ci.use[have.coord]
      vec.Urban[flag.vec] <- overlay2dis
      
      if(i %% 20 == 0)
      {','\n',
      
      paste0('save(vec.Urban,file = \'./sub',m,'/VecUrban',m,'.RData\')'),'\n',
      paste0('fwrite(data.frame(vec.Urban),row.names = F, \'./sub',m,'/VecUrban',m,'.csv\')'),'\n',
      '}','\n',
      
      'cat( paste0("i = ", i),"\\n",
      
      append = T, file = ',paste0('\'./sub',m,'/VecUrban',m,'.txt\')'),'\n',
      
      '
      rm(overlay2dis)
      }','\n',
  paste0('save.image(\'./UrbanSub',m,'.RData\')'),'\n',
  
  
  file = paste0('overUrban',m,'.R'))
}


#################################
## Run the beyond scripts 
################################

#######################################################################
## Processing of Urban overlay results from high performance computing
## over.urban: integer, polygon objectID
#########################################################################
dir.create('YourPath/OverUrbanResult')
setwd('YourPath/OverUrbanResult')

over.urban <- NULL
for(l in 1:15)
{
  load(paste0('UrbanSub',l,'.RData'))
  st <- (loop.mt[l,1]-1)*1000+1
  en <- (loop.mt[l,2])*1000
  if(l == 15) {en <- nrow(df.bio)}
  over.urban <- c(over.urban,vec.Urban[st:en])
  print(i)
}

## Urban overlay type
obj.center  <- urban2dis@data$OBJECTID[which(urban2dis@data$grid_code == 3)]
obj.cluster <- urban2dis@data$OBJECTID[which(urban2dis@data$grid_code == 2)]

## overUrbanType- NA,3,2: NA, center, cluster
overUrbanType <- rep(NA,length(over.urban))
overUrbanType[which(over.urban %in% obj.center)]  <- 3
overUrbanType[which(over.urban %in% obj.cluster)] <- 2

###########################################
## Adjust coordinates overlay and results
###########################################

# Extract records
no.sub <- which(df.bio$Adj.Coord == 1)
subuse <- df.bio[no.sub,]

## location.adj: 0, not used; 1, adj used
location.adj <- rep(0,nrow(df.bio))
location.adj[no.sub] <- 1

coord.adj <- data.frame(array(NA,dim = c(length(no.sub),2)))
names(coord.adj) <- c('lon','lat')

# Creat a dataframe to save records' coordinates
num.adj <- which(df.bio$Adj.Coord == 1)
df.adj <- data.frame(num.row = num.adj, Type = df.bio$Adj.TransTypr[num.adj], 
                     Lon = df.bio$decimalLongitude[num.adj], 
                     Lat = df.bio$decimalLatitude[num.adj],
                     AdjLon = NA, AdjLat = NA)

for(i in 1:nrow(df.adj))
{
  coord.adj <- Coord8.single(x = df.adj$Lon[i], y = df.adj$Lat[i],type = df.adj$Type[i])
  df.adj[i,5:6] <- coord.adj
}

seq.ci    <- seq(1,nrow(df.adj),1000)
mt.ci     <- array(NA,dim=c(length(seq.ci),2))

for(i in 1:length(seq.ci))
{
  mt.ci[i,1] <- seq.ci[i]
  mt.ci[i,2] <- seq.ci[i]+999
  if(i == length(seq.ci)) { mt.ci[i,2] <- nrow(df.adj)}
}

use.mt <- mt.ci
adj.loc.over <- rep(NA,nrow(df.adj))

# Urban area detection
for(i in 1:nrow(use.mt))
{
  ci.use  <- use.mt[i,1]:use.mt[i,2]
  df.occ  <- df.adj[ci.use,]
  
  pro.shp <- SpatialPointsDataFrame(df.occ[,c("AdjLon","AdjLat")],
                                    data = df.occ[,c("num.row","Type")],
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
  pro.shp@proj4string <- urban2dis@proj4string
  
  overuse2            <- over(pro.shp,urban2dis)
  adj.loc.over[ci.use]    <- overuse2$OBJECTID
  print(i)
}

over.type <- rep(NA,nrow(df.adj))
over.type[which(adj.loc.over %in% obj.center)]  <- 3
over.type[which(adj.loc.over %in% obj.cluster)] <- 2

## Adj.over.urban: integer, objectID
## Adj.overUrbanType: NA,3,2
Adj.over.urban    <- rep(NA,length(over.urban))
Adj.overUrbanType <- rep(NA,length(over.urban))

Adj.over.urban[no.sub]    <- adj.loc.over
Adj.overUrbanType[no.sub] <- over.type 

## add.space: data.frame, save all variables results
add.space <- data.frame(over.urban, overUrbanType,Adj.over.urban,Adj.overUrbanType)

## write results
df.bio[,names(add.space)] <- add.space


###########################################################################################################
## df.bio is the integrated dataset with different indexes presting the data cleaning results
###########################################################################################################


