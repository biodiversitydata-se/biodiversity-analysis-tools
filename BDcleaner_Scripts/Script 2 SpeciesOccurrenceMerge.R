##############################################################################
## Script 2. Species occurrence data & data preprocessing
## Author: Jing Jin, Jun Yang*
## Manuscript: BDcleaner: a workflow for cleaning taxonomic and geographic errors
##             in occurrence data archived in biodiversity databases
################################################################################

###############################################
## Part 1. Data collection and processing
###############################################

##############################
## For GBIF
#############################
library(dismo)
library(rgbif)

genus.name   <- as.character(unique(df.spe$Genus))
df.spe$occ   <- rep(NA,nrow(df.spe))
fail.species <- NULL
big          <- NULL

for(i in 1:length(genus.name))
{
  row.use <- which(df.spe$Genus==genus.name[i])
  if(dir.exists(genus.name[i])==F) dir.create(genus.name[i])
  
  for(j in 1:length(row.use))
  {
    genus <- as.character(df.spe$Genus[row.use[j]])
    species <- as.character(df.spe$Species[row.use[j]])
    sp2 <- paste0(species,'*')
    
    down <- try(gbif(genus,sp2,download=T,ntries=20, geo=F))
    
    
    if(class(down)== "try-error"){
      fail.species <- c(fail.species,df.spe$ID[row.use[j]])
      
      cat(df.spe[row.use[j],1],'\n',
          append = T, file = 'YourPath/GBIF/fail.txt'
      )
      
    }else 
    {
      namefile <- paste(df.spe$Genus[row.use[j]],df.spe$Species[row.use[j]])
      try(fwrite(down,paste0('./',genus.name[i],'/',namefile,'.csv'),row.names = F),silent = T)
      save(down,file = paste0('./',genus.name[i],'/',namefile,'.RData'))
      
      if(is.null(down))
      {
        df.spe$occ[row.use[j]] <- 0
      } else
      {
        df.spe$occ[row.use[j]] <- nrow(down)
        if( nrow(down) == 18000)
        {
          big <- c(big,df.spe$ID[row.use[j]])
          
          cat(df.spe[row.use[j],1],'\n',
              append = T, file = 'YourPath/GBIF/big.txt'
          )
        }
      }
    }
    
  }
  
}


# For bid dataset
no.fail <- scan(file = 'big.txt')

lon.int <- c(seq(-180,0,45),seq(0,45,5)[-1],seq(45,180,15)[-1])
df.lon <- data.frame(start = rep(NA,length(lon.int)-1), end = rep(NA,length(lon.int)-1)) 
df.lon$start <- lon.int[1:length(lon.int)-1]
df.lon$end <- lon.int[-1]

down <- list()
down.le <- list()

failure <- NULL
flag <- 0 
for(j in no.fail)
{
  genus <- as.character(df.spe$Genus[j])
  sp    <- as.character(df.spe$Species[j])
  sp2   <- paste0(sp,'*')
  down[[j]]  <- list()
  down.le[[j]] <- rep(NA,22)
  
  for(p in 1:nrow(df.lon))
  {
    ext     <- extent(c(t(df.lon[p,]),-90,90))
    dl.try  <- try(gbif(genus,sp2,download=T,ntries=20,ext = ext))
    
    if(class(dl.try) != "try-error")
    {
      if(is.null(dl.try)) 
      {down.le[[j]][p] <- 0} else
      {
        down.le[[j]][p] <- nrow(dl.try)
        down[[j]][[p]] <- dl.try
        rm(dl.try)
      }
    } else
    {
      flag <- 1
    }
    
    if(flag == 1) break
  }
  
  flag <- 0
}


# Combine big dataset
setwd('YourPath')

for(j in no.fail)
{
  no.occ <- sum(down.le[[j]])
  
  if (is.na(no.occ))
  {next} else
  {
    no.df <- which(down.le[[j]]!=0)
    name.list <- sapply(down[[j]], names)
    
    name.init <- name.list[[no.df[1]]]
    
    for(m in no.df[-1])
    {
      name.init <- union(name.init,name.list[[m]])
    }
    
  }
  
  ## Total number of occ
  df.use <- data.frame(array(NA,dim=c(no.occ,length(name.init))))
  names(df.use) <- name.init
  
  cal.int <- down.le[[j]][no.df]
  cal.df  <- data.frame(start = rep(NA,length(no.df)), end = rep(NA,length(no.df)))
  for(p in 1:length(no.df))
  {
    cal.df$end[p] <- sum(cal.int[1:p])
    
    if(p == 1)
    {
      cal.df$start[p] <- 1
    } else
    {
      cal.df$start[p] <- cal.df$end[p-1]+1
    }
  }
  
  ## Save occ
  for(p in 1:length(no.df))
  {
    df.use[cal.df$start[p]:cal.df$end[p],c(name.list[[no.df[p]]])] <- down[[j]][[no.df[p]]][,c(name.list[[no.df[p]]])]
  }
  
  genus <- as.character(df.spe$Genus[j])
  
  if(dir.exists(genus)==F) dir.create(genus)
  
  namefile <- paste(df.spe$Genus[j],df.spe$Species[j])
  try(fwrite(df.use,paste0('./',genus,'/',namefile,'.csv'),row.names = F),silent = T)
  save(df.use,file = paste0('./',genus,'/',namefile,'.RData'))
  
}

sele.var <- c('scientificName','countryCode','decimalLatitude','decimalLongitude','year')

for(j in 1:length(genus.name))
{
  row.use <- which((df.spe$Genus==genus.name[j]) & df.spe$occ>0)
  if(length(row.use)==0) next
  
  filepath <- paste0('YourPath/',genus.name[j],'/')
  file.ex  <- df.spe$SciName[row.use]
  
  for(k in 1:length(file.ex))
  {
    load(paste0(filepath,file.ex[k],'.RData'))
    
    binduse  <- subset(df.try, select = sele.var)
    no.row   <- nrow(binduse)
    
    if(j == 1 & k == 1)
    {
      fwrite(binduse,'YourPath/GBIF.csv',row.names = F)
    } else
    {
      fwrite(binduse,'YourPath/GBIF.csv',append = T,row.names = F)
    }
    
  }
}




#############################
## For BIEN
#############################
library(BIEN)
library(data.table)

setwd('YourPath/BIENSpecies/')

genus.name   <- as.character(unique(df.spe$Genus))
df.spe$bienocc   <- rep(NA,nrow(df.spe))
fail.species.bien <- NULL


for(i in 1:length(genus.name))
{
  row.use <- which(df.spe$Genus==genus.name[i])
  if(dir.exists(genus.name[i])==F) dir.create(genus.name[i])
  
  for(j in 1:length(row.use))
  {
    genus <- as.character(df.spe$Genus[row.use[j]])
    species <- as.character(df.spe$Species[row.use[j]])
    sp2 <- paste0(species,'*')
    
    
    Bien_full <- try(BIEN_occurrence_species(species = df.spe$SciName[row.use[j]],
                                             cultivated = T,only.new.world = F,
                                             all.taxonomy = T,native.status = T,
                                             observation.type = T,
                                             political.boundaries = T),silent = T)
    
    
    if(class(Bien_full)== "try-error"){
      fail.species.bien <- c(fail.species.bien,df.spe$ID[row.use[j]])
      
      cat(df.spe[row.use[j],1],'\n',append = T, file = '../Bienfail.txt')
      
    }else 
    {
      namefile <- df.spe$SciName[row.use[j]]
      try(fwrite(Bien_full,paste0('./',genus.name[i],'/',
                                  namefile,'.csv'),row.names = F),silent = T)
      save(Bien_full,file = paste0('./',genus.name[i],'/',
                                   namefile,'.RData'))
      
      df.spe$bienocc[row.use[j]] <- nrow(Bien_full)
    }
    
  }
  
}

sele.var <- c('verbatim_scientific_name','country','latitude','longitude','date_collected')

for(j in 1:length(genus.name))
{
  row.use <- which((df.spe$Genus==genus.name[j]) & df.spe$bienocc>0)
  if(length(row.use)==0) next
  
  filepath <- paste0('YourPath/',genus.name[j],'/')
  file.ex  <- df.spe$SciName[row.use]
  
  for(k in 1:length(file.ex))
  {
    load(paste0(filepath,file.ex[k],'.RData'))
    
    binduse  <- subset(df.try, select = sele.var)
    no.row   <- nrow(binduse)
    
    if(j == 1 & k == 1)
    {
      fwrite(binduse,'YourPath/BIEN.csv',row.names = F)
    } else
    {
      fwrite(binduse,'YourPath/BIEN.csv',append = T,row.names = F)
    }
    
  }
}



#############################
## For iDigBio
#############################

library(ridigbio)

# setwd('H:\\Distribution\\iDigBio')
# genus.name <- unique(df.spe$Genus)
# for(i in 1:length(genus.name))
# {
#   dir.create(paste0('./',genus.name[i],'/'))
# }
fail.idig <- NULL

for(i in 1:length(genus.name))
{
  row.use <- which(df.spe$Genus==genus.name[i])
  
  for(j in 1:length(row.use))
  {
    use.ch <- str_trim(tolower(df.spe$SciName[row.use[j]]))
    down.occ <- try(idig_search_records(rq=list(scientificname=use.ch)),silent = T)
    
    if(class(down.occ)[1] == "try-error")
    {
      fail.idig <- c(fail.idig,df.spe$ID[row.use[j]])
      cat(df.spe[row.use[j],1],'\n',append = T, file = '../idigbiofail.txt')
      rm(down.occ);rm(use.ch)
    }else
    {
      namefile <- df.spe$SciName[row.use[j]]
      save(down.occ,file = paste0('./',genus.name[i],'/',namefile,'.RData'))
      rm(down.occ);rm(use.ch)
    }
    
  }
}

df.spe$occ_idigbio <- NA
df.spe$occ2_idigbio <- NA

setwd('YourPath\\iDigBio')
for(i in 1:length(genus.name))
{
  if(i == 1309) next
  row.use <- which(df.spe$Genus==genus.name[i])
  
  for(j in 1:length(row.use))
  {
    namefile <- df.spe$SciName[row.use[j]]
    file.use <- paste0('./',genus.name[i],'/',namefile,'.RData')
    if(file.exists(file.use))
    {
      load(file.use)
      down.occ2 <- unique(down.occ)
      
      if(nrow(down.occ)>0)
      {
        fwrite(down.occ2,'../idigbioDF_uninei.csv',append = T,row.names = F)
        df.spe$occ2_idigbio[row.use[j]] <- nrow(down.occ2)
      } else
      {
        df.spe$occ2_idigbio[row.use[j]] <- 0
      }
    }
    
  }
}


#############################
## For BISON
#############################
df.Bison <- fread('YourPath/bison_1489789055_csv/
                   bison.csv', 
                   encoding = 'UTF-8')
sp.sciname <- df.Bisonk$scientificName
sp.pro     <- str_trim(tolower(attain.latin(sp.sciname)))
use.row    <- which(sp.pro %in% name.spe)
df.Bison2  <- df.Bison[use.row,]

fwrite('YourPath/BISON2/bison2_tree.csv',row.names = F)

#############################
## For SpeciesLink
#############################
df.SpLink <- fread('YourPath/speciesLink_all_50886_20190421044811.txt/
                   speciesLink_all_50886_20190421044811.txt', 
                   encoding = 'UTF-8')
sp.sciname <- df.SpLink$scientificname
sp.pro     <- str_trim(tolower(attain.latin(sp.sciname)))
use.row    <- which(sp.pro %in% name.spe)
df.SpLink2 <- df.SpLink[use.row,]

fwrite(df.SpLink2,'YourPath/SpeciesLink/splink_tree.csv',row.names = F)



#############################
## For ALA
#############################
setwd('YourPath/records-2019-04-12')

# For all ALA records download from ALA website
temp <- list.files(pattern = '^records*')

for(i in 1:length(temp))
{
  
  df.ala <- fread(temp[i],encoding = 'UTF-8')
  
  name.use1 <- tolower(str_trim(attain.latin(df.ala$scientificName)))
  name.use2 <- tolower(str_trim(df.ala$species))
  
  num1 <- which(name.use1 %in% name.spe)
  num2 <- which(name.use2 %in% name.spe)
  row.df.ala <- sort(unique(c(num1,num2)))
  
  sub.ala <- df.ala[row.df.ala,]
  sub.ala2 <- unique(sub.ala)
  
  fwrite(sub.ala2,'YourPath/ALA/ala_tree_2names.csv',row.names = F, append = T)
  
  rm(df.ala); rm(sub.ala); rm(sub.ala2);rm(name.use1); rm(name.use2); rm(row.df.ala)
}
df.ala <- fread('YourPath/ALA/ala_tree_2names.csv',encoding = 'UTF-8')



#############################
## For RAINBIO
#############################
df.rainbio <- fread('YourPath/rainbio_published/published_database/RAINBIO.csv',
                    encoding = 'UTF-8')
df.rainbio$spe <- tolower(str_trim(df.rainbio$tax_sp_level))
sub.rainbio    <- subset(df.rainbio, spe %in% name.spe)
fwrite(sub.rainbio,'YourPath/RainBio_tree.csv',row.names = F)


#############################
## For BioTime
#############################
df.biotime <- fread('YourPath/BioTIMEQuery02_04_2018/BioTIMEQuery02_04_2018.csv',
                    encoding = 'UTF-8')
df.biotime$spe <- tolower(str_trim(df.biotime$GENUS_SPECIES))
sub.biotime    <- subset(df.biotime, spe %in% name.spe)
fwrite(sub.biotime,'YourPath/BioTime_tree.csv',row.names = F)






###############################################
## Part 2. Data integration
###############################################
library(countrycode) # for contry code list
library(passport) # for country code match use
df.code <- codelist

#########################
# GBIF
##########################
df.GBIF <- fread('YourPath/GBIF.csv',encoding = 'UTF-8') %>% unique()
db.GBIF <- subset(df.GBIF,select = c('scientificName','countryCode',
                                    'decimalLatitude','decimalLongitude','year'))
db.GBIF$FlagCountry <- 1
db.GBIF$ID     <- 1:nrow(db.GBIF)
db.GBIF$DB     <- 'GBIF'

#########################
# BIEN
##########################
df.BIEN <- fread('YourPath/BIEN.csv',encoding = 'UTF-8') %>% unique()
table.BIEN <- df.BIEN$country %>% table
name.BIEN  <- df.BIEN$country %>% table %>% names
match1     <- name2iso(name = name.idig,df.code = df.code)
code.BIEN  <- rep(NA,nrow(df.BIEN))
t.use2     <- df.BIEN$country
for(i in which(!is.na(match1))) {code.BIEN [which(t.use2 == name.BIEN[i])] <- match1[i]}
df.BIEN$countryCode <- code.BIEN


library(lubridate)
use.year <- df.BIEN$date_collected %>% as.Date %>% year
df.BIEN$year <- use.year


db.BIEN <- subset(df.BIEN,select = c('verbatim_scientific_name','countryCode',
                                     'latitude','longitude','year'))
names(db.BIEN) <- names(db.GBIF)
db.BIEN$FlagCountry <- 1
db.BIEN$ID     <- 1:nrow(db.BIEN)
db.BIEN$DB     <- 'BIEN'

#########################
# ALA
##########################
df.ALA     <- fread('YourPath/ALA/ala_tree_2names.csv', 
                    encoding = 'UTF-8')
db.ALA     <- subset(df.ALA,select = c('scientificName','countryCode',
                                       'decimalLatitude','decimalLongitude','year'))
ala.c1 <- df.ALA$countryCode
ala.c2 <- df.ALA$country

iso3166.code <- read.csv('ISO3166.csv',header = T) ## datasource: iso 3166 online browser (script1)
name.use <- ala.c2[which(ala.c1 == "" & ala.c2 != "")]
code.mat <- rep(NA,length(name.use))

pro1 <- name.use %>% unique %>% tolower %>% str_trim
num1 <- which(pro1 %in% tolower(iso3166.code$`English short name`))
name.mat <- rep(NA,length(pro1))
for(i in num1)
{
  loop1 <- which(tolower(iso3166.code$`English short name`) == pro1[i])
  name.mat[i] <- as.character(iso3166.code$`Alpha-2 code`[loop1])
}


pro2 <- pro1[-num1]
match3 <- name2iso(name = pro2,df.code = df.code)
match3[18] <- 'FM'
name.mat[-num1] <- match3
t.use <- name.use %>% tolower

for(i in 1:length(name.mat)) {code.mat[which(t.use == pro1[i])] <- name.mat[i]}
db.ALA$countryCode[which(ala.c1 == "" & ala.c2 != "")] <- code.mat
db.ALA$FlagCountry <- 1

db.ALA$ID     <- 1:nrow(db.ALA)
db.ALA$DB     <- 'ALA'

###########################
# BISON
###############################
library(data.table)
df.Bison   <- fread('YourPath/BISON/bison2_tree.csv',header=T)
db.Bison   <- subset(df.Bison,select=c('providedScientificName','countryCode',
                                       'decimalLatitude','decimalLongitude','year'))
db.Bison$FlagCountry <- 1
names(db.Bison) <- names(db.Bison)
db.Bison$ID     <- 1:nrow(db.Bison)
db.Bison$DB     <- 'Bison'

###########################
# iDigBio 
###############################
df.iDig    <- fread('YourPath/iDigBio/idigbioDF_uninei.csv', 
                    encoding = 'UTF-8')

db.idig    <- subset(df.iDig,select=c('scientificname','country','geopoint.lat',
                                      'geopoint.lon','datecollected'))

table.idig <- db.idig$country %>% table
name.idig  <- db.idig$country %>% table %>% names
match1     <- name2iso(name = name.idig,df.code = df.code)
code.idig  <- rep(NA,nrow(db.idig))
t.use2     <- db.idig$country

for(i in which(!is.na(match1))) {code.idig[which(t.use2 == name.idig[i])] <- match1[i]}

db.idig$countryCode <- code.idig
db.idig$FlagCountry <- 1
db.idig$FlagCountry[which(is.na(code.idig))] <- df.iDig$country[which(is.na(code.idig))]

# For year
library(lubridate)
use.year <- db.idig$datecollected %>% as.Date %>% year
db.idig$year <- use.year


db.iDigBio <- subset(db.idig, select = c("scientificname", "countryCode", 
                                         "geopoint.lat", "geopoint.lon", "year" , "FlagCountry"))
names(db.iDigBio) <- names(db.ALA)

db.iDigBio$ID     <- 1:nrow(db.iDigBio)
db.iDigBio$DB     <- 'iDigBio'


###########################
# SpeciesLink
###############################
df.SpeLink <- fread('YourPath/SpeciesLink/splink_tree.csv', encoding = 'UTF-8')
db.SpeLink <- subset(df.SpeLink,select=c('scientificname','scientificnameauthor','country',
                                         'latitude','longitude','yearcollected'))
table.spel <- db.SpeLink$country %>% table
name.spel  <- db.SpeLink$country %>% table %>% names
match2     <- name2iso(name = tolower(name.spel),df.code = df.code)
sum(table.spel[which(is.na(match2))])

code.spel <- rep(NA,nrow(db.SpeLink))
t.use3    <- db.SpeLink$country
for(i in which(!is.na(match2))) {code.spel[which(t.use3 == name.spel[i])] <- match2[i]}

db.SpeLink$countryCode <- code.spel
db.SpeLink$FlagCountry <- 1
db.SpeLink$FlagCountry[which(is.na(code.spel))] <- db.SpeLink$country[which(is.na(code.spel))]

v1 <- db.SpeLink$scientificname; v2 <- db.SpeLink$scientificnameauthor
Name.pro <-  str_trim(paste(v1,v2))
db.SpeLink$Name <- Name.pro

db.SpeciesLink <- subset(db.SpeLink, select = c("Name", "countryCode", "latitude", 
                                                "longitude", "yearcollected" , "FlagCountry"))
names(db.SpeciesLink) <- names(db.ALA)

db.SpeciesLink$ID     <- 1:nrow(db.SpeciesLink)
db.SpeciesLink$DB     <- 'SpeciesLink'


###########################
# BioTime
###############################
df.BioTime <- fread('YourPath/BioTime/BioTime_tree.csv', 
                    encoding = 'UTF-8')
db.BioT <- subset(df.BioTime,select=c('GENUS_SPECIES','LATITUDE',
                                      'LONGITUDE','YEAR')) # no countrycode
db.BioT$countryCode <- NA
db.BioT$FlagCountry <- 0
db.BioTime <- subset(db.BioT,select = c("GENUS_SPECIES", "countryCode", "LATITUDE",
                                        "LONGITUDE", "YEAR" , "FlagCountry"))
names(db.BioTime) <- names(db.ALA)
db.BioTime$ID     <- 1:nrow(db.BioTime)
db.BioTime$DB     <- 'BioTime'


###########################
# RainBio
###############################
df.rainbio <- fread('YourPath/RainBio/RainBio_tree.csv', encoding = 'UTF-8')
db.rain <- subset(df.rainbio,select=c('tax_sp_level','countryCode',
                                      'decimalLatitude','decimalLongitude')) # no year
db.rain$FlagCountry <- 1
db.rain$YEAR        <- '9999_rainbio'
db.RAINBIO <- subset(db.rain,select = c("tax_sp_level", "countryCode", "decimalLatitude", 
                                        "decimalLongitude", "YEAR" , "FlagCountry"))
names(db.RAINBIO) <- names(db.ALA)
db.RAINBIO$ID     <- 1:nrow(db.RAINBIO)
db.RAINBIO$DB     <- 'RAINBIO'


###############################
# Integrated dataset
###############################
list.database <- list(db.ALA,db.Bison,db.iDigBio,db.SpeciesLink,db.BioTime,db.RAINBIO,db.GBIF,db.BIEN)

bind.use <- NULL
for(i in 1:length(list.database)) {bind.use <- rbind(bind.use,list.database[[i]])}

# remove duplicate records
dup.num <- duplicated(bind.use[,-c(7:8)])
df.bio <- bind.use[which(!dup.num),]

save(df.bio,file = 'YourPath/dfBIO.RData')
