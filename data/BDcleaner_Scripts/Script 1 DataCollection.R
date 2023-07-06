##############################################################################
## Script 1. Data collection and preprocessing
## Author: Jing Jin, Jun Yang*
## Note: This script contain the methods to collect part of data used in our 
##       manuscript.Some of them were collected used automated collection 
##       methods for scientific research purposes.
## Manuscript: BDcleaner: a workflow for cleaning taxonomic and geographic errors
##             in occurrence data archived in biodiversity databases
################################################################################


############################################
## Data 0. Global Tree Search database
############################################

## The most important dataset in our analysis
#  Tree name list data
#  Download website: https://tools.bgci.org/global_tree_search.php

# load packages
library(data.table)
library(stringi)
library(stringr)
library(magrittr)

# read GTS data
df.spe <- fread('YourPath/global_tree_search_trees_1_2.csv',header = T)[,1:2]

# Spilit the taxonname into genus name and specific epithet
pro1 <- df.spe$`Taxon name`  |> 
  str_trim |> str_split(pattern = " ") |> 
  unlist |> matrix(ncol=2,byrow = T) |> as.data.frame

# Add genus and species name into this dataframe
df.spe[,c('Genus','Species')] <- pro1


############################################
## Data 1. Botanical Garden Information
## Methods: Automatic data collection
## Important output: garden.inf data.frame
############################################
library(RCurl)
library(httr)
library(XML)

# Automatic collection of botanical garden information in BGCI garden search
# Read Garden name and correspongding websites
url.list    <- 'https://tools.bgci.org/garden_search.php?action=Find&ftrCountry=All&ftrKeyword=&x=50&y=21'
html.use    <- GET(url.list)
doc1        <- htmlParse(html.use)
garden.inf  <- readHTMLTable(doc1)[[3]]
garden.link <- getHTMLLinks(doc1)
link.no     <- which(str_detect(garden.link,'garden.php?'))
link.use    <- paste0('http://www.bgci.org/',garden.link[link.no])

# Read detail information of each botanical garden
path.list <- list()
for(i in 1:length(link.use))
{
  url1        <- link.use[i]
  doc         <- htmlParse(url1)
  
  xpath1      <- xpathSApply(doc = doc,path = "//div[@class='meta-heading']//li",xmlValue)
  save.no1    <- which(str_detect(xpath1,'Latitude')) # save latitude
  save.no2    <- which(str_detect(xpath1,'Longitude')) # save longitude
  save.no3    <- which(str_detect(xpath1,'Altitude')) # save altitude
  save.no     <- unique(c(save.no1,save.no2,save.no3))
  save.chr    <- xpath1[save.no]
  
  path.list[[i]] <- save.chr # save results in list
  
}

# String processing
use.list <- list()
for(i in 1:length(link.use))
{
  trans1 <- str_trim(path.list[[i]])
  if(length(trans1) == 0)
  {
    use.list[[i]] <- NULL
    next
  }
  trans2        <- str_split(trans1,pattern = '\t')[[1]]
  use.vec       <- trans2[which(trans2 != "")]
  use.list[[i]] <- use.vec
}


# Transform information format to data.frame
garden.inf$Latitude <- NA
garden.inf$Longitude <- NA

for(i in 1:length(link.use))
{
  if(is.null(use.list[[i]])) next
  t1 <- str_split(use.list[[i]],':')
  t2 <- lapply(t1, str_trim)
  
  for(k in 1:length(t2))
  {
    if(t2[[k]][1] != "Latitude" & t2[[k]][1] != "Longitude")
    {
      next
    } else
    {
      if(t2[[k]][1] == "Latitude")
      {
        garden.inf$Latitude[i]  <- t2[[k]][2]
      } else
      {
        garden.inf$Longitude[i] <- t2[[k]][2]
      }
    }
  }
  
}

# write the BGCI garden information
write.csv(garden.inf,'YourPath/GardenInformation.csv',row.names = F)

############################################
## Data 2. The Plant list
## Methods: Automatic data collection
## Important output: taxonbind data.frame
############################################
library(XML)
library(RCurl)

## TPL data website
URL <- 'http://www.theplantlist.org/1.1/browse/-/-/'
a <- getHTMLLinks(URL)
detect.use <- paste0('/1.1/browse/',LETTERS)

# Output available genus and corresponding family in TPL 
de.result <- NULL
for(i in 1:length(detect.use))
{
  de <- which(str_detect(a,detect.use[i]))
  de.result <- c(de.result,de)
}

uni.de <- sort(unique(de.result))
split.use <- str_split(a[uni.de],pattern = '/')
TPL.genus <- data.frame(GenusID = c(1:length(uni.de)), Family = rep(NA,length(uni.de)),
                        Genus = rep(NA,length(uni.de))) # save TPL information

for(i in 1:length(uni.de))
{
  sp.length <- length(split.use[[i]])
  TPL.genus$Family[i] <- str_trim(split.use[[i]][sp.length-2])
  TPL.genus$Genus[i]  <- str_trim(split.use[[i]][sp.length-1])
}

# Download species name list of each genus (save as *.csv file)
for(i in 1:nrow(TPL.genus))
{
  url.dl.genuse <- paste0('http://www.theplantlist.org/tpl1.1/search?q=',
                          TPL.genus$Genus[i],'&csv=true')
  download.file(url.dl.genuse,
                destfile = paste0('YourPath/',TPL.genus$Genus[i],'.csv'))
}


## TPL data processing
library(data.table)
setwd('YourPath')

# Bind all csv files into "taxonbind" data.frame
temp <- list.files(pattern = '*.csv')
TPLgenus <- unique(TPL.genus$Genus)
Genus.dl <- data.frame(GenusID = c(1:length(TPLgenus)), Genus = TPLgenus, 
                       Number = rep(NA,length(TPLgenus)))
file.name <- paste0(TPLgenus,'.csv')
for(i in 1:length(TPLgenus))
{
  temp2 <- fread(file = file.name[i],encoding = 'UTF-8')
  Genus.dl$Number[i] <- nrow(temp2)
}

taxonbind <- NULL
for(i in 1:length(temp))
{
  buse <- fread(file = temp[i],encoding = 'UTF-8')
  taxonbind <- rbind(taxonbind,buse)
}


## Taxonbind data.frame processing

# Bionominal name
taxonbind$latin2   <- paste(taxonbind$Genus,taxonbind$Species)

# Complete name
name.col <- c('Genus','Species','Infraspecific rank',"Infraspecific epithet","Authorship")
pro.bind  <- rep(NA,nrow(taxonbind))
pro.bind2 <- rep(NA,nrow(taxonbind))

for(i in 1:nrow(taxonbind))
{
  duse        <- as.data.frame(taxonbind[i,c('Genus','Species','Infraspecific rank',"Infraspecific epithet","Authorship")])
  no.paste    <- which(!is.na(duse[1,]) & duse[1,] != "")
  pro.bind[i] <- str_trim(paste(duse[1,no.paste], collapse = " "))
}

# consider special punctation
pro.bind2.1 <- str_replace_all(pro.bind,"\\(","")
pro.bind2.2 <- str_replace_all(pro.bind2.1,"\\)","")

# Author information
re.aut <- tolower(taxonbind$Authorship)
re.pro1 <- str_replace_all(re.aut,"\\(","")
re.pro2 <- str_replace_all(re.pro1,"\\)","")

# Output variables
# Bindname: intial species name, not removing special punctations
#           and consider both capitalization and lowercase
# bindname: Removing special punctations
#           Transform string to lowercase
taxonbind$AuthorInf   <- ProBlank(re.pro2)
taxonbind$bindname    <- ProBlank(tolower(pro.bind2.2))
taxonbind$Bindname    <- pro.bind2.2 

# save(taxonbind,file = 'YourPath')


#############################################
## Data 3. Tropicos, CoL, and GBIF name list
## Important output: garden.inf data.frame
#############################################

# COL AND TROPICOS Name processing

library("taxize")
col.name.list <- list()
tropicos.name <- list()

save.pt <- seq(0,nrow(df.spe),10000)
for(i in 1:nrow(df.spe))
{
  r1 <- try(synonyms(df.spe$SciName[i],db= 'col',rows = 1),silent = T)
  r2 <- synonyms(df.spe$SciName[i],db= 'tropicos',rows = 1)
  
  col.name.list[[i]] <- r1
  tropicos.name[[i]] <- r2
  
  rm(r1);rm(r2)
  
}

# null.col <- NULL
# null.tro <- NULL
# for(i in 1:nrow(df.spe))
# {
#   if(class(col.name.list[[i]][[1]]) == 'list') null.col <- c(null.col,i)
#   if(class(tropicos.name[[i]][[1]]) == 'list') null.tro <- c(null.tro,i)
# }

list.taxon.col <- list()
list.taxon.tro <- list()

for(i in 1:nrow(df.spe))
{
  col.use <- col.name.list[[i]][[1]]
  
  if(is.na(col.use)) next
  
  if(df.spe$SciName[i] == names(col.name.list[[i]])[1])
  {
    list.taxon.col[[i]] <- unique(c(df.spe$SciName[i], col.use$name))
  }
}

col.name2 <- tolower(str_trim(unlist(list.taxon.col)))
name.spe  <- unique(c(spe.name2,col.name2))

# GBIF name processing
library(RSQLite)
con <- dbConnect(dbDriver("SQLite"),dbname='J:\\GBIF_occ_all_180717\\gbif.sqlite')
dbListTables(con)
gbif.all   <- dbGetQuery(con, 'SELECT * FROM gbif')
gbif.plant <- subset(gbif.all, kingdom == 'Plantae' & taxonRank %in% 
                       c("species","subspecies","variety","form","unranked"))


gbif.plant.latinname <- tolower(str_trim(attain.latin(gbif.plant$scientificName)))
g1.row    <- which(gbif.plant.latinname %in% tolower(df.spe$SciName))
pro1      <- gbif.plant[g1.row,]
acc.id    <- setdiff(unique(pro1$acceptedNameUsageID),"")
taxon.id  <- pro1$taxonID
use.id    <- unique(c(acc.id,taxon.id))
row.out   <- which(gbif.plant$taxonID %in% use.id | gbif.plant$acceptedNameUsageID %in% use.id)

gbif.pro  <- gbif.plant[row.out,]
pro.latin <- str_trim(tolower(unique(attain.latin(gbif.pro$scientificName))))
name.spe2 <- unique(c(name.spe,pro.latin))
write.csv(gbif.pro ,'YourPath/GBIF_namelist_tree.csv',
          row.names = F)

tropicos.name <- list()
for(i in 1:nrow(df.spe))
{
  # r1 <- try(synonyms(df.spe$SciName[i],db= 'col',rows = 1),silent = T)
  r2 <- try(synonyms(df.spe$SciName[i],db= 'tropicos',rows = 1),silent = T)
  if(class(r2) == 'try-error')
  {
    repeat
    {
      Sys.sleep(120)
      r2 <- try(synonyms(df.spe$SciName[i],db= 'tropicos',rows = 1),silent = T)
      if(class(r2) != 'try-error') break
    }
  }
  tropicos.name[[i]] <- r2
  if(i %in% save.pt) save(tropicos.name,file = 'YourPath/temp_tropicos.RData')
}
save.image('YourPath/COL_TROName.RData')

list.taxon.tro <- list()
for(i in 1:nrow(df.spe))
{
  tro.use <- tropicos.name[[i]][[1]]
  
  if(is.na(tro.use)) next
  list.taxon.tro[[i]] <- unique(c(df.spe$SciName[i], setdiff(tro.use$scientificname,'no syns found')))
}

tro.name2 <- tolower(str_trim(unlist(list.taxon.tro)))
name.spe  <- unique(c(name.spe2,tro.name2))

gbif.pro <- fread('YourPath/GBIF_namelist_tree.csv',encoding = 'UTF-8')
row.out1 <- which(gbif.pro$taxonID %in% use.id & !(gbif.pro$acceptedNameUsageID %in% use.id))
row.out2 <- which(gbif.pro$acceptedNameUsageID %in% use.id)

gbif.list <- list()
gbif.spe.name <- tolower(str_trim(attain.latin(gbif.pro$scientificName)))
tree.latin    <- df.spe$SciName |> tolower

for(i in 1:nrow(df.spe))
{
  id.spe  <- which(gbif.spe.name == tree.latin[i])
  
  if(length(id.spe) > 0)
  {
    id.num  <- gbif.pro$taxonID[id.spe]
    id.acc  <- setdiff(gbif.pro$acceptedNameUsageID[id.spe],"")
    id.num2 <- unique(c(id.num,id.acc))
    row.use <- which(gbif.pro$acceptedNameUsageID %in% id.num2 | gbif.pro$taxonID %in% id.num2)
    
    name.use <- gbif.pro$scientificName[row.use]
    name.lat <- unique(gbif.spe.name[row.use])
    
    gbif.list[[i]] <- name.lat
  }
}

list.taxon.col  <- lapply(list.taxon.col, tolower)
list.taxon.gbif <- gbif.list
list.taxon.tro  <- lapply(list.taxon.tro, tolower)

list.bind <- list()

name.df <- df.spe$SciName |> tolower

for(i in 1:nrow(df.spe))
{
  if(i < 60095)  list.bind[[i]] <- unique(c(list.taxon.col[[i]],list.taxon.gbif[[i]],list.taxon.tro[[i]],name.df[i]))
  if(i == 60095) list.bind[[i]] <- unique(c(list.taxon.col[[i]],list.taxon.tro[[i]],name.df[i]))
  if(i > 60095)  list.bind[[i]] <- unique(c(list.taxon.col[[i]],list.taxon.tro[[i]],name.df[i]))
}




############################################
## Data 4. World map data
############################################

# world map 1: World map in maptools R package (Bivand & Lewin-Koh 2018)
library(maptools)
data("wrld_simpl")

# world map 2: Natural earth 10m resolution world map 
#             (https://www.naturalearthdata.com/downloads/10m-cultural-vectors/)

# world map 3: GADM world political boundary (https://gadm.org/data.html)


############################################
## Data 5. Centroids of administrative areas
############################################

# GeoNames global cities (https://www.geonames.org/);
# Global cities and province distribution (Zizka et al. 2019), please refer to this paper
# GADM world political boundary (https://gadm.org/data.html)


############################################
## Data 6. Coordinates of institutions
############################################

# Global database of biodiversity institutions (Zizka et al. 2019), please refer to this paper

############################################
## Data 7. ISO-3166 COUNTRY CODES
############################################

# ISO-3166 country code list online browser (https://www.iso.org/obp/ui/#search/code/).
# R package codelist
library(countrycode)
iso.code <- codelist 

############################################
## Data 8. GHS settlement grid data
############################################

# GHS settlement grid data (https://ghsl.jrc.ec.europa.eu/download.php)
# Note: we have trasformed the tiff file into shapefile in ArcGis.
#       Code 1: Urban Center; Code 2: Urban cluster; Code 3: Rural
#       we treat code 1 & 2 as urban classification