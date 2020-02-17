# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
setwd("~/Dropbox/EnvAge/Final files 08:10:18")

#install dependency for metagear - EBImage
#source("https://bioconductor.org/biocLite.R") 
biocLite("EBImage")

#install.packages('revtools')
library(revtools)
#install.packages('topicmodels')
library(topicmodels)
#install.packages('plyr')
library(plyr)
#install.packages('reshape')
library(reshape)
#install.packages('metagear')
library(metagear)
#install.packages('randomForest')
library(randomForest)

# LOAD BIBLIOS AND CLEAN --------------------------------------------------

#import SCOPUS bib file and convert to df
scopus_biblio_1<-as.data.frame(read_bibliography('scopus_1to676_130818.bib')) #import biblio
final_scopus_biblio.080818<-scopus_biblio_1 #compile final biblio (in this case just one file)
final_scopus_biblio.080818$abstract[which(is.na(final_scopus_biblio.080818$abstract))]<-"" #turn NA abstracts to empty space
final_scopus_biblio.080818$source<-'scopus' #label with source for downstream tracking

#import WoS bib file and convert to df
wos_biblio_1<-as.data.frame(read_bibliography('wos_1to500_130818.bib'))
wos_biblio_2<-as.data.frame(read_bibliography('wos_501-832_130818.bib'))
colnames.dup<-count(c(colnames(wos_biblio_1),colnames(wos_biblio_2)))
matchingcols<-as.character(colnames.dup$x[colnames.dup$freq==2])
final_wos_biblio.080818<-rbind(wos_biblio_1[,matchingcols],wos_biblio_2[,matchingcols])
final_wos_biblio.080818$abstract[which(is.na(final_wos_biblio.080818$abstract))]<-""
final_wos_biblio.080818$source<-'wos'

#import pubmed nbib files and convert to df
pubmed_biblio_1<-as.data.frame(read_bibliography('pubmed_1to200_130818.nbib'))
pubmed_biblio_2<-as.data.frame(read_bibliography('pubmed_201to400_130818.nbib'))
pubmed_biblio_3<-as.data.frame(read_bibliography('pubmed_400to541_130818.nbib'))
colnames.dup<-count(c(colnames(pubmed_biblio_1),colnames(pubmed_biblio_2),colnames(pubmed_biblio_3)))
matchingcols<-as.character(colnames.dup$x[colnames.dup$freq==3])
final_pubmed_biblio.080818<-rbind(pubmed_biblio_1[,matchingcols],pubmed_biblio_2[,matchingcols],pubmed_biblio_3[,matchingcols])
final_pubmed_biblio.080818$abstract[which(is.na(final_pubmed_biblio.080818$abstract))]<-""
final_pubmed_biblio.080818$source<-'pubmed'

#bind all dataframes tgoether without pubmed to create master biblio
#colnames.dup.all<-count(c(colnames(final_scopus_biblio.080818),colnames(final_wos_biblio.080818)))
#matchingcols<-as.character(colnames.dup.all$x[colnames.dup.all$freq==2])
#masterbiblio.080818<-rbind(final_scopus_biblio.080818[,matchingcols],final_wos_biblio.080818[,matchingcols])
#count(is.na(masterbiblio.080818$abstract))

#bind all dataframes tgoether with pubmed to create master biblio
colnames.dup.all<-count(c(colnames(final_scopus_biblio.080818),colnames(final_wos_biblio.080818),colnames(final_pubmed_biblio.080818)))
matchingcols<-as.character(colnames.dup.all$x[colnames.dup.all$freq==3])
masterbiblio.080818<-rbind(final_scopus_biblio.080818[,matchingcols],final_wos_biblio.080818[,matchingcols],final_pubmed_biblio.080818[,matchingcols])
count(is.na(masterbiblio.080818$abstract))

# REVTOOLS PROCESSING -----------------------------------------------------

#remove duplicates
dupsfinal.080818<-find_duplicates(masterbiblio.080818)
masterbiblio.080818_dupsremoved<-extract_unique_references(dupsfinal.080818)
#remove duplicates double check by removing duplicated dois
doi.dups<-which(duplicated(masterbiblio.080818_dupsremoved$doi, incomparables = NA))
masterbiblio.080818_dupsremoved<-masterbiblio.080818_dupsremoved[-doi.dups,]

# TRAINING: FINDING RELEVANT LITERATURE ------------------------------------------

#GOAL: We have lots of literature, dominated by medical stuff. 
#Some of this medical stuff is about hospital air, so is *kind of* relevant.
#Let's use the abstract screener to separate up this literature into (don't use air terms as in all):

#YES IT'S RELEVANT! (YES): DIRECTLY Environmental (non-hospital) AMR e.g. case studies, disease descriptions
#MAYBE OF RELEVANCE (MAYBE): DIRECTLY About Environmental AMR, but within a hospital context.
#NO RELEVANCE (NO): Not DIRECTLY about environmental AMR AND/OR AIR

#check what effort files you have before removing
list.files(pattern = "effort")
#do.call(file.remove, list(list.files(pattern="effort")))
list.files(pattern = "effort")
set.seed(1235)
topic.temp<-masterbiblio.080818_dupsremoved[sample(1:nrow(masterbiblio.080818_dupsremoved)),]
set.seed(1235)
masterbiblio.080818_dupsremoved<-masterbiblio.080818_dupsremoved[sample(1:nrow(masterbiblio.080818_dupsremoved)),]

#set up stuff for metagear abstract screening - correcting and matching column titles and order etc
colnames(topic.temp)
topic.temp$pages<-"1-?"
LPAGES_UPAGES<-t(as.data.frame(strsplit(topic.temp$pages,'-')))
colnames(LPAGES_UPAGES)<-c('LPAGES','UPAGES')
topic.temp<-cbind(topic.temp,LPAGES_UPAGES)
metagear_topic.temp<-topic.temp[,c("author","year","title","journal","volume","LPAGES","UPAGES","doi","abstract")]
colnames(metagear_topic.temp)<-colnames(example_references_metagear)

#screen abstracts of and classify manually - remember to save as you go along in metagear console
effort_distribute(metagear_topic.temp, initialize = TRUE, reviewers = "envamr", save_split = TRUE)
abstract_screener('effort_envamr.csv', aReviewer = 'envamr')
effort_distribute(metagear_topic.temp, initialize = TRUE, reviewers = "test", save_split = TRUE)
abstract_screener('effort_test.csv', aReviewer = 'test')

#check file wrote OK and reload it as 'screenset'
list.files(pattern = "effort")
screenset<-read.csv('effort_envamr.csv')

#create sub-bibliographies of relevant and clinical data
masterbiblio.080818_dupsremoved$relevant<-screenset$INCLUDE
masterbiblio.080818_dupsremoved_relevant<-masterbiblio.080818_dupsremoved[masterbiblio.080818_dupsremoved$relevant=='YES',]
masterbiblio.080818_dupsremoved_clinical<-masterbiblio.080818_dupsremoved[masterbiblio.080818_dupsremoved$relevant=='MAYBE',]

# EXPLORE ENVIRONMENTS ----------------------------------------------------

set.seed(1234)
relevant_DTM<-as.data.frame(make_DTM(masterbiblio.080818_dupsremoved_relevant))
clinical_DTM<-as.data.frame(make_DTM(masterbiblio.080818_dupsremoved_clinical))

#manual identification of topics with abstract screener - also refined Yes/nos
title_abs<-paste(masterbiblio.080818_dupsremoved_relevant$title,masterbiblio.080818_dupsremoved_relevant$abstract)
title_abs_clin<-paste(masterbiblio.080818_dupsremoved_clinical$title,masterbiblio.080818_dupsremoved_clinical$abstract)

#make some vectors with words to hunt for in titles and abstracts to get sub-bibliographies of different topics
env_hosp<-c("hospital","ward","nosocomial")
env_dental<-c("dental","dentist")
env_urban<-c("city","cities","smog","Nanjing","Istanbul","park","parks")
env_rural<-c("rural")
env_pigfarm<-c("pig","pigs","swine")
env_poultryfarm<-c("poultry","chicken","chickens","turkey","turkeys")
env_cowfarm<-c("cow","cows","beef","cattle","dairy","calve","calves")
env_agriculture<-c("manure","slurry","solid waste","irrigation","reclaimed","crop","crops")
env_untouched<-c("glacier","lake","space station")
env_wastewater<-c("wastewater","waste water","WWTP","sewage")
env_indoor<-c("indoor")
env_outdoor<-c("outdoor")
env_residential<-c("residential","home","homes","household","households","flat","flats")
env_homeless<-c("homeless")
env_compost<-c("compost","composting")
env_office<-c("office")
env_factory<-c("factory")
env_arable<-c("vineyard","orchard","plantation","rice")
env_atmosphere<-c("atmosphere")
env_kindergarten<-c("kindergarten","day care","nursery")
env_stable<-c("hay","stable","stables","horse")
env_dairy<-c("dairy")
env_sawmill<-c("sawmill")
env_aircon<-c("air conditioning")
env_foodprep<-c("food prep")
env_church<-c("church")
env_aircraft<-c("aircraft")
env_pharma<-c("pharma")

#use vectors to subset biblios based on above classifications - not perfect but through trial and error have checked that this doesn't miss stuff out and is imperfect in that it includes a few irrelevant stuff to, which are removed manually later in Excel
biblio_hosp<-masterbiblio.080818_dupsremoved_clinical[unique(grep(paste(env_hosp,collapse="|"),title_abs_clin,ignore.case = TRUE)),]
biblio_dental<-masterbiblio.080818_dupsremoved_clinical[unique(grep(paste(env_dental,collapse="|"),title_abs_clin,ignore.case = TRUE)),]
biblio_urban<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_urban,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_rural<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_rural,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_pigfarm<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_pigfarm,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_poultryfarm<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_poultryfarm,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_cowfarm<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_cowfarm,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_agriculture<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_agriculture,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_untouched<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_untouched,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_wastewater<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_wastewater,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_indoor<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_indoor,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_outdoor<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_outdoor,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_residential<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_residential,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_homeless<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_homeless,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_compost<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_compost,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_office<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_office,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_factory<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_factory,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_arable<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_arable,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_atmosphere<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_atmosphere,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_kindergarten<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_kindergarten,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_dairy<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_dairy,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_sawmill<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_sawmill,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_aircon<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_aircon,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_foodprep<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_foodprep,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_aircraft<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_aircraft,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_church<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_church,collapse="|"),title_abs,ignore.case = TRUE)),]
biblio_pharma<-masterbiblio.080818_dupsremoved_relevant[unique(grep(paste(env_pharma,collapse="|"),title_abs,ignore.case = TRUE)),]

# BROWSE ABSTRACTS OF A SUB-TOPIC -----------------------------------------

set.seed(1235)
topic.temp<-biblio_pigfarm

#set up stuff for metagear abstract screening - correcting and matching column titles and order etc
set.seed(1235)
colnames(topic.temp)
topic.temp$pages<-"1-?"
LPAGES_UPAGES<-t(as.data.frame(strsplit(topic.temp$pages,'-')))
colnames(LPAGES_UPAGES)<-c('LPAGES','UPAGES')
topic.temp<-cbind(topic.temp,LPAGES_UPAGES)
metagear_topic.temp<-topic.temp[,c("author","year","title","journal","volume","LPAGES","UPAGES","doi","abstract")]
colnames(metagear_topic.temp)<-colnames(example_references_metagear)

#screen abstracts of a sub-bibliography e.g. pig
effort_distribute(metagear_topic.temp, initialize = TRUE, reviewers = "pig", save_split = TRUE)
abstract_screener('effort_pig.csv', aReviewer = 'pig')

#write it to a csv
#write.csv(biblio_wastewater, 'biblio_wastewater.csv')
