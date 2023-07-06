##############################################################################
## Script 4. Functions
## Author: Jing Jin, Jun Yang*
## Manuscript: BDcleaner: a workflow for cleaning taxonomic and geographic errors
##             in occurrence data archived in biodiversity databases
################################################################################


################################
## Taxon Filter Functions
###############################

## attain bionomial nomenclature
attain.latin <- function(vector)
{
  split      <- str_split(vector,pattern = " ")
  return.vec <- rep(NA,length(vector))
  for(i in 1:length(vector)) {return.vec[i] <- paste(split[[i]][1:2],collapse = " ")}
  return(return.vec)
}

# Fot TPL dataframe (taxonbind)
attain.inf <- function(vector)
{
  vec.use <- vector[,c(8,9,10)]
  vec2    <- paste0(vec.use,collapse = " ")
  vec.mod <- processAuthor(vec2)
  return(vec.mod)
}

# String processing of author information
processAuthor <- function(author)
{
  rep.s1 <- str_replace_all(author,"[\\(\\)]","")
  rep.s2 <- str_replace_all(rep.s1," ","")
  return(rep.s2)
}

## Accepted Name Function
acc.match <- function(tpl.judge)
{
  tpl.acc     <- subset(tpl.judge,`Taxonomic status in TPL` == 'Accepted')
  df.res      <- data.frame(TPL_Match = NA, TPL_Status = NA,TPL_AccID = NA, 
                            TPL_Conf  = NA, TPL_Line = NA)
  flag.status <- 0
  
  if(nrow(tpl.acc) != 0)
  {
    a.no         <- which(tpl.acc$conf == max(tpl.acc$conf))
    df.res[1,]   <- tpl.acc[a.no,c('latin2','Taxonomic status in TPL','Accepted ID','conf','Line')]
    df.res[,1:3] <- apply(df.res[,1:3],2,as.vector)
    
    flag.status <- 1
  }
  
  return(list(df.res,flag.status))
}

## Synonyn Name Function
syn.match <- function(tpl.judge, GBIF_Name, Latin2Name, taxonbind)
{
  tpl.syn     <- subset(tpl.judge,`Taxonomic status in TPL` == 'Synonym')
  flag.status <- 0
  df.res      <- data.frame(TPL_Match = NA, TPL_Status = NA,TPL_AccID = NA, 
                            TPL_Conf  = NA, TPL_Line = NA)
  if(nrow(tpl.syn) != 0)
  {
    
    aut.gbif1 <- str_trim(str_replace_all(GBIF_Name,Latin2Name,""))
    aut.GBIF  <- processAuthor(aut.gbif1)
    
    for(r in 1:nrow(tpl.syn)) {tpl.syn$authorInf[r] <- attain.inf(tpl.syn[r,])}
    
    LD.judge <- sort(adist(aut.GBIF, tpl.syn$authorInf),index.return = T)[[2]]
    
    for(s in LD.judge)
    {
      key.acc <- tpl.syn$`Accepted ID`[s]
      no.acc  <- which(taxonbind$ID == key.acc)
      
      if(length(no.acc) == 0) next
      
      if(taxonbind$latin2[no.acc] == Tree_name)
      {
        df.res[1,]   <- tpl.syn[s,c('latin2','Taxonomic status in TPL','Accepted ID','conf','Line')]
        df.res[,1:3] <- apply(df.res[,1:3],2,as.vector)
        
        flag.status  <- 2
        break
      }
    }
    
  }
  
  return(list(df.res,flag.status))
}

## Unresolved AND Misapplied Name Function
unrmis.match <- function(tpl.judge, GBIF_Name, Latin2Name, taxonbind)
{
  tpl.neg <- subset(tpl.judge,`Taxonomic status in TPL` %in% c("Unresolved","Misapplied"))
  df.res  <- data.frame(TPL_Match = NA, TPL_Status = NA,TPL_AccID = NA,
                        TPL_Conf  = NA, TPL_Line = NA)
  
  if(nrow(tpl.neg) != 0)
  {
    aut.gbif1 <- str_trim(str_replace_all(GBIF_Name,Latin2Name,""))
    aut.GBIF  <- processAuthor(aut.gbif1)
    flag <- 0
    for(r in 1:nrow(tpl.neg)) {tpl.neg$authorInf[r] <- attain.inf(tpl.neg[r,])}
    
    LD.judge <- sort(adist(aut.GBIF, tpl.neg$authorInf),index.return = T)[[2]]
    
    for(s in LD.judge)
    {
      key.acc.neg <- tpl.neg$`Accepted ID`[s]
      
      ## If the name don't have a synonym name
      if((str_length(key.acc.neg) == 0) | is.na(key.acc.neg)) 
      {
        if (tpl.neg$latin2[s] == Tree_name)
        {
          df.res[1,]   <- tpl.neg[s,c('latin2','Taxonomic status in TPL','Accepted ID','conf','Line')]
          df.res[,1:3] <- apply(df.res[,1:3],2,as.vector)
          flag <- 1
        } 
        
      } else
      {
        ## If the name have a synonym name
        no.acc.neg  <- which(taxonbind$ID == key.acc.neg)
        
        if(length(no.acc.neg) == 0) next 
        
        if(taxonbind$latin2[no.acc.neg] == Tree_name)
        {
          df.res[1,]   <- tpl.neg[s,c('latin2','Taxonomic status in TPL','Accepted ID','conf','Line')]
          df.res[,1:3] <- apply(df.res[,1:3],2,as.vector)
          flag <- 2
        }
        
      }
      
      if(flag != 0) break
      
      
    }
  }
  
  return(df.res)
}

## Unresolver_GTS Function
unracc.match <- function(df.spe,tpl.judge, GBIF_Name, Latin2Name, taxonbind)
{
  aut.gbif1 <- str_trim(str_replace_all(GBIF_Name,Latin2Name,""))
  aut.GBIF  <- processAuthor(aut.gbif1)
  df.res    <- data.frame(TPL_Match = NA, TPL_Status = NA,TPL_AccID = NA, 
                          TPL_Conf  = NA, TPL_Line = NA)
  
  flag <- 0
  
  for(t in 1:nrow(tpl.judge)) {tpl.judge$autInf[t] <- attain.inf(tpl.judge)}
  LD.judge <- sort(adist(aut.GBIF, tpl.judge$authorInf),index.return = T)[[2]]
  
  for(s in LD.judge)
  {
    acc.id <-  tpl.judge$`Accepted ID`[s]
    line.jud <- which(taxonbind$ID == acc.id)
    
    if(length(line.jud) == 0) next
    
    if(taxonbind$latin2[line.jud] == Tree_name)
    {
      df.res[1,]   <- tpl.judge[s,c('latin2','Taxonomic status in TPL','Accepted ID','conf','Line')]
      df.res[,1:3] <- apply(df.res[,1:3],2,as.vector)
      flag <- 1
      break
    }
  }
  
  return(df.res)
}

attain.spname <- function(judge.use)
{
  pro1   <- lapply(str_split(judge.use," "), str_trim)
  loop.j <- length(judge.use)
  
  detectuse <- c('f.','subsp.','var.','var')
  
  df.save <- data.frame(Genus.name = rep(NA,loop.j), Spe.name = rep(NA,loop.j), 
                        Author.name = rep(NA,loop.j), SpeType = rep(0,loop.j))
  for(j in 1:loop.j)
  {
    jud.pro <- which(pro1[[j]] %in% detectuse)
    
    if(length(jud.pro)!=0 )
    {
      
      df.save$Genus.name[j]  <- pro1[[j]][1]
      df.save$Spe.name[j]    <- str_trim(paste(pro1[[j]][2:(jud.pro[1]+1)],collapse = " "))
      df.save$Author.name[j] <- str_trim(paste(pro1[[j]][-c(1:(jud.pro[1]+1))],collapse = " "))
      df.save$SpeType[j]     <- 1
      
    }else
    {
      df.save$Genus.name[j]  <- pro1[[j]][1]
      df.save$Spe.name[j]    <- pro1[[j]][2]
      df.save$Author.name[j] <- str_trim(paste(pro1[[j]][-c(1,2)],collapse = " "))
    }
  }
  
  return(df.save)
}

ProBlank <- function(ch)
{
  pro1 <- str_replace(ch,"  "," ")
  pro2 <- str_replace(pro1,"   "," ")
  
  return(pro2)
}

####################################
## Geo Filter complex use functions
###################################

# Judge a vector: is.zero(x,y)
Coord.Notzero <- function(lat,lon)
{
  len <- length(lat)
  result <- rep(NA,len)
  
  loop.use <- c(1:len)
  
  no.NA <- union(which(is.na(lat)),which(is.na(lon)))
  loop  <- setdiff(loop.use,no.NA)
  
  if (length(loop) != 0)
  {
    
    for(i in loop)
    {
      if (lat[i] == 0 | lon[i] == 0)
      {
        result[i] <- 0
      } else {result[i] <-  1}
    }
  }
  
  return(result)
}

# Creat coordinates (8, +-, reverse)(transformation types of coordinates)
diff.coord <- function(x,y)
{
  result <- data.frame(xcoord = rep(NA,8), ycoord = rep(NA,8))
  result$xcoord <- c(rep(x,2),rep(0-x,2))
  result$ycoord[c(1,3)] <- rep(y,2)
  result$ycoord[c(2,4)] <- rep(0-y,2)
  result$xcoord[c(5:8)] <- result$ycoord[c(1:4)]
  result$ycoord[c(5:8)] <- result$xcoord[c(1:4)]
  
  return(result)
}


##############################
## Overlay Functions
##############################

## Return loop interval (split intervals)
sp.mt <- function(st, en, interval)
{
  seq.use <- seq(from = st,to = en,by = interval)
  
  result <- data.frame(array(NA,dim = c(length(seq.use),2)))
  names(result) <- c('start','end')
  
  for(i in 1:length(seq.use))
  {
    result[i,1] <- seq.use[i]
    
    if(i < length(seq.use))
    {
      result[i,2] <- seq.use[i]+interval-1
    } else
    { result[i,2] <- en}
  }
  
  return(result)
}


# output different correction types of coordinates
# and its' corresponding overlay results
Coord8.Out <- function(worldmap_use,coordx,coordy)
{
  #x-lon y-lat
  d1 <- cbind(coordx,-coordy)
  d2 <- cbind(-coordx,coordy)
  d3 <- cbind(-coordx,-coordy)
  d4 <- cbind(coordy,coordx)
  d5 <- cbind(coordy,-coordx)
  d6 <- cbind(-coordy,coordx)
  d7 <- cbind(-coordy,-coordx)
  
  d.list <- list(d1,d2,d3,d4,d5,d6,d7)
  rm(list = paste0('d',1:7))
  
  over.list <- list()
  
  for(d in 1:length(d.list))
  {
    caluse <- SpatialPoints(d.list[[d]])
    caluse@proj4string <- worldmap_use@proj4string
    overresult <- over(caluse,worldmap_use)
    over.list[[d]] <- overresult
    rm(list = c('caluse','overresult'))
  }
  
  return(over.list)
  
}

Coord8.single <- function(x,y,type=c(1:7))
{
  #x-lon y-lat
  
  d1 <- c(x,-y)
  d2 <- c(-x,y)
  d3 <- c(-x,-y)
  d4 <- c(y,x)
  d5 <- c(y,-x)
  d6 <- c(-y,x)
  d7 <- c(-y,-x)
  
  d.use <- get(paste0('d',type))
  return(d.use)
}

# Match the coodinates tansformation type
Coord.match.trans <- function(ISO.rec, match.overlist, NameVec)
{
  df.coord1  <- data.frame(array(NA,dim = c(nrow(match.overlist[[1]]),7)))
  names(df.coord1) <- paste0('Trans',1:7)
  
  for(d in 1:7) {df.coord1[,d] <- match.overlist[[d]][,NameVec] |> as.character}
  match.trans <- rep(NA,nrow(df.coord1))
  
  for(i in 1:nrow(df.coord1))
  {
    replace.use <- which(df.coord1[i,] == ISO.rec[i])
    if(length(replace.use) != 0) match.trans[i] <- replace.use[1]
  }
  
  return(match.trans)
  
}

attain.sciname <- function(ch.name)
{
  pro.ch1 <- str_trim(ch.name)
  pro.ch2 <- str_replace_all(pro.ch1,"\\(","")
  pro.ch3 <- str_replace_all(pro.ch2,"\\)","")
  return(pro.ch3)
}

match.fun <- function(vector,df.iso)
{
  
  ret.vec <- rep(NA,length(vector))
  for(i in 1:length(vector))
  {
    ret.vec[i] <-  df.iso$ISO2[which(df.iso$NAME == vector[i])]
  }
  
  return(ret.vec)
}

# Attain number of species
get.speNo <- function(name.use, list.spe)
{
  
  num <- 0
  
  for(i in 1:length(list.spe))
  {
    ju1 <- which(name.use %in% list.spe[[i]])
    if(length(ju1 != 0)) num <- num+1
  }
  
  
  return(num)
}


## Convert country name to ISO-3166 country code
name2iso <- function(name,df.code)
{
  num.match.list <- list()
  name       <- str_trim(tolower(name))
  parse.name <- parse_country(name, to = "iso2c", how = c("regex", "google", "dstk"),
                              language = c("en", "de"), factor = is.factor(x))
  name.loop <- name[which(is.na(parse.name))]
  name2     <- name.loop |> str_replace_all(pattern = '[:,.?!@;]$',replacement = "") |> str_replace_all(pattern = '^[:,.?!@;]',replacement = "")
  
  
  for(j in 1:length(name.loop))
  {
    num.match <- NULL
    for(i in 1:ncol(df.code))
    {
      num.match <- c(num.match,which(df.code[,i] == tolower(name.loop[j])))
    }
    
    if(length(num.match) == 0)
    {
      for(i in 1:ncol(df.code))
      {
        num.match <- c(num.match,which(df.code[,i] == name2[j]))
      }
    }
    
    num.match.list[[j]] <- num.match
  }
  
  
  num.judge <- lapply(num.match.list,unique)
  st1 <- NULL
  st2 <- NULL
  
  
  loop.vec <- rep(NA,length(name.loop))
  for(j in 1:length(name.loop))
  {
    ju.number <- length(num.judge[[j]])
    
    if(ju.number > 1)  st2 <- c(st2,j)
    if(ju.number == 0) st1 <- c(st1,j)
    
    if(ju.number == 1) loop.vec[j] <- toupper(df.code$iso2c[num.judge[[j]][1]])
  }
  
  parse.name[which(is.na(parse.name))] <- loop.vec
  
  return(parse.name)
}







