library(rvest)
library(dplyr)

## Scrapes relavent infomrational summaries of CRAN libraries from liabrary CRAN webpage 


#gets a list of packages from a previous run
 PackageVector <- read.csv("./CRAN_data/dependancies.csv")
 PackageVector <- PackageVector[grep("cran.r-project",PackageVector$webpage),2]
# 

PackageVectorChoice <- c("Hypervolume") # This is used to start if there has been no previous run can be any package

#sees if the package list has alread been summarise
PackagesDone <- read.csv("./CRAN_data/Summary.csv")
PackagesDone <- PackagesDone$package
if((length(PackagesDone)+length(PackageVectorChoice))!=0){
  PackageVector <- c(PackageVector,PackageVectorChoice)
PackageVector <- grep(paste0(PackagesDone,collapse = "|"),PackageVector,value = TRUE,invert = TRUE)

}

PackageVector

for(i in 1:length(PackageVector)){
pck <- PackageVector[i] #which package to look at
URL <- read_html(paste0("https://cran.r-project.org/web/packages/", #CRAN desciroption webiste
                        pck,
                        "/index.html")) 

pack_desc <- URL %>%    #ddexcription of packacge
    html_nodes("p")%>%
    html_text()
  

pack_table <- URL %>%      #extracts webtables of desciption downloads and dependancies
              html_table()

published <- pack_table %>%
              magrittr::extract(1)            #extracts package description out
published <- published[[1]]$X2[grep("Published",published[[1]]$X1)]

published <- data.frame(package = pck,     #name of package
                        type = "published", #published type date
                        date = published)   #date

pack_tablenodes <- URL %>%          #probably can be removed but not just yet
                  html_nodes("table")

#If there are dependancies (packages which use this package) get some info
if(length(pack_tablenodes) == 3){      
  dependancy_packs <- pack_tablenodes %>%                
    magrittr::extract2(3) %>%
                    html_table()
  dependancy_packs <- trimws( unlist( strsplit( unlist( dependancy_packs$X2), #list packages
                                             ","))) 
  dependancy_packsWebs<- pack_tablenodes %>%  #extract webapages for packages
    magrittr::extract2(3) %>%           #NB can include none CRAN packages
                    html_nodes("a") %>%
                    html_attr("href")
  
  #replace .. with CRAN if appropriate
  dependancy_packsWebs <- gsub("[.][.]/","https://cran.r-project.org/web/packages/",
                               dependancy_packsWebs) 
  dependancies <- data.frame(package = pck,                 #what package
                             dependancy = dependancy_packs, #what dependancies
                             webpage = dependancy_packsWebs)#dependancy website
} else {
  dependancies <- data.frame(package = pck,                 #what package
                             dependancy = NA, #what dependancies
                             webpage = NA)#dependancy website
}

webPages <- pack_tablenodes %>% # vignettes and manual websites
  magrittr::extract2(2) %>%
  html_nodes("a") %>%
  html_attr("href")

vignsL <- grep("vignettes",webPages,value=TRUE) #are there vignetted

if(length(vignsL)!=0){
  vigns <- grep("vignettes",webPages,value=TRUE) #if there are fignettes list them
  
  vigns <- data.frame(package = pck,             #what package
                    pdf = "vignette",            #is it an vignette
                    Webpage = vigns)             #webpage for vignette
  
  
  vigns <- rbind(vigns,data.frame(package = pck, #bind manual weblink
               pdf = "manual",
               Webpage = paste0(pck,".pdf")))
  
} else {
  #no Vigenette create summary for manual only
  vigns <- data.frame(package = pck,                       
                                  pdf = "manual",
                                  Webpage = paste0(pck,".pdf"))
                      }

vigns$Webpage <- paste0("https://cran.r-project.org/web/packages/", #turn vignette and manual link into a weblink
                        pck,
                        "/",vigns$Webpage)

ARCweb <- grep("Archive", webPages, value = TRUE) #is there an archive

if(length(ARCweb)!=0){
  
updates <- read_html(ARCweb)%>% #retrieve update history
  html_table()%>%
  magrittr::extract2(1)%>%
  select("Last modified")

updates[updates==""] <- NA  #remove rows tih no info
updates <- na.omit(updates)

updates <- data.frame(package = pck,    #name of package
                      type = "updates", #updates ttype
                      date = updates$`Last modified`)# when were the updates

updates <- rbind(published, updates) #bind to publihed data

} else {
  
  updates <- published #no updates or archived 

}
#the following are csv tables 

#SUMMARY.CSV summaries information that each package has
Summa <- read.csv("./CRAN_data/Summary.csv") 

write.csv(rbind(data.frame(package = pck,              #name of package
           vignette=length(grep("vignette",vigns$pdf)),#How many vignettes are there
           dependancies=nrow(dependancies),            #How many dependecies are listed
           updatesN = nrow(updates),                   #How many updates (this includes date of publication)
           lastUpdate = updates$date[nrow(updates)],   #when was last update
           description = gsub("\n","",pack_desc[1])),Summa), #package description
     paste0("./CRAN_data/Summary.csv"),
     row.names = FALSE)

#VIGNETTE.CSV info on vignettes and manuals available for each package
Vignets <- read.csv("./CRAN_data/vignettes.csv")
write.csv(rbind(vigns, 
                Vignets),
          "./CRAN_data/vignettes.csv",
          row.names = FALSE)

#DEPENDANCIES.CSV info on dependancies for each package
deps <- read.csv("./CRAN_data/dependancies.csv")
write.csv(rbind(deps,dependancies),"./CRAN_data/dependancies.csv",
          row.names = FALSE)

#uppdate information
UPPS <- read.csv("./CRAN_data/updates.csv")
write.csv(rbind(UPPS,updates),"./CRAN_data/updates.csv",
          row.names = FALSE)

print(pck)
}
