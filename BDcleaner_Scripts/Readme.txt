Manuscript:BDcleaner: a workflow for cleaning taxonomic and geographic errors in occurrence data archived in biodiversity databases; 
Author: Jing Jin and Jun Yang*

This folder contains five R scripts, which involve the collection and integration of tree occurrece data from 8 biodiversity databases, data cleaning and analysis in this study.If you have any problems with the code, please contact us and we will do our best to assist you.

###########################################
## Simple description of each R script
###########################################

(1) Script 1 Data Collection:

This script contains the code to collect the auxiliary data that needs to be used in the data cleaning process. Some of the data is easy to download (such as Global tree search database, GeoNames database), so we only provided the download URL. For that data that is difficult to download directly (such as TPL, BGCI garden search), we used automatic data collection process to complete.

(2) Script 2 SpeciesOccurrenceMerge:

This script contains the code to collect the tree species occurrece data from 8 different biodiversity databases. Besides, we also provided the code to integrate the multi-source databases following the methods in our manuscript.

(3) Script 3 DataCleaning:

This script is the core part of our study, which contains the data cleaning process of space and taxonomy dimension. 

(4) Script 4 Functions:

This script contains all functions used in our analysis. Please run this script before you use any script. 

(5) Script 5 Analysis and Plot:

This script contains the analysis process in our manuscript.








