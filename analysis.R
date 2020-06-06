# START -------------------------------------------------------------------
rm(list=ls())
set.seed(1234)

# INSTALL NECESSARY PACKAGES ----------------------------------------------

#install.packages("plyr") - count function needed (not to be confused with base count function)
library(plyr)
#install.packages('meta') - for metaprop meta-analysis
library(meta)
#install.packages('scale') - for viewing colour palettes
library(scales)
#install.packages('randomForest') - for random forest analysis (randomForest function etc.)
library(randomForest)
#install.packages('caret') - for varImp function to display varImps of randomForest object
library(caret)
#install.packages('rpart') - for regression tree (rpart function)
library(rpart)
#install.packages('ggsci') - for regression tree (rpart function)
library(ggsci)

# LOAD FUNCTIONS ----------------------------------------------------------

#load function to calculate standard error (to be used in graph error bars)
se <- function(x) sd(x, na.rm = T)/sqrt(length(x))

# LOAD DATA ---------------------------------------------------------------

#load full bibliographic dataframe of all studies used in the lit review (though not necessarily quantitatively)
fullbiblio<-read.csv('analysis/input_data/biblio_all_final.csv')
#load antimicrobial susceptibility testing dataset
MAR.allstudies<-read.csv('analysis/input_data/AMSdataset_051018.csv', stringsAsFactors = F)

#i've classified each of the antibiotics into 5 levels of antibiotic classification, as best as I can.
#these lines strip out each of these classifications into separate vectors to keep them safe so we can work with them later is we want, before we revert to the highest level for the actual analysis
ABclass1<-as.factor(as.character(MAR.allstudies[1,20:ncol(MAR.allstudies)]))
ABclass2<-as.factor(as.character(MAR.allstudies[2,20:ncol(MAR.allstudies)]))
ABclass3<-as.factor(as.character(MAR.allstudies[3,20:ncol(MAR.allstudies)]))
ABclass4<-as.factor(as.character(MAR.allstudies[4,20:ncol(MAR.allstudies)]))

#re-load dataframe with just highest level antibiotic classification as the heading of the antibiotic columns 
#this line removes the unwanted antibiotic classification rows and the first column describing what each row is.
MAR.allstudies<-read.csv('analysis/input_data/AMSdataset_051018.csv', stringsAsFactors = F)[-c(1:4),-1]

#this line removes sample.name column and the wind direction column, which just describes the meaning of the distance.m column
MAR.allstudies<-MAR.allstudies[,-c(1,8)]

#these lines convert the number of strains tested (N) and number of antibiotics tested (B) columns to numeric class
MAR.allstudies$N<-as.numeric(MAR.allstudies$N)
MAR.allstudies$B<-as.numeric(MAR.allstudies$B)

#there was one empty row (froms gibbs2004's study), in which no strains were found and therefore none tested; we remove it here for the purposes of analysis
MAR.allstudies<-MAR.allstudies[MAR.allstudies$N!=0,]
#creates a vector of antimicrobial (biotic) susceptibility (AMS) testing results columns; acts as a shortcut in later code
MAR.allstudies.ABcols<-MAR.allstudies[,17:ncol(MAR.allstudies)]
#here we check whether the AMS columns are numeric - they're not due to the way we read in the dataframe (to avoid some bugs to do with aggregating the data later)
apply(MAR.allstudies.ABcols, 2, function(x) is.numeric(x))
#here we check whether the AMS columns are numeric - they're not...
MAR.allstudies.ABcols<-as.data.frame(apply(MAR.allstudies.ABcols, 2, function(x) as.numeric(x)))
#here we check how many lines are not equal to NA (i.e. we have datapoints for) for each antibiotic 
colSums(MAR.allstudies.ABcols,na.rm = T)
#here we make sure the number of antibiotics tested column is right by calculating, for each row, the row sum excluding NAs (i.e. not tested) and making the B column this
MAR.allstudies$B<-rowSums(apply(MAR.allstudies.ABcols, 2, function(x) !is.na(x)))

# GENERAL TRENDS (ALL STUDIES) ---------------------------------
#this section of the code calculates summary statistics for the '2.1: General Trends' section of the report, summarising coverage of different topics etc. in the literature

#check total number of studies excluding duplicates (i.e. studies covering multiple categories)
#row.id column tells us which row in original large bibliographic dataset it came from i.e. whether it is unique

fullbiblio_dupsremoved<-fullbiblio[match(unique(fullbiblio$row.id),fullbiblio$row.id),]
nrow(fullbiblio_dupsremoved)

#2.1.1 Taxonomic coverage
#how many studies of staphylococci/MRSA?
count(fullbiblio_dupsremoved$MRSA.SA)
#how many studies of Enterobacteriacae/ESBL?
count(fullbiblio_dupsremoved$ESBL.ENT)
#how many studies of fungi?
count(fullbiblio_dupsremoved$FUNGI)
#how many studies of viruses
count(fullbiblio_dupsremoved$VIRUS)

#how many studies were culture-based - really roundabout way of calculating it...
#count number of AMS-testing studies ('YES', 'NO', or NOT POSSIBLE' i.e. did testing but results were not reported well enough to do quantiative comparison) on bacteria and make a dataframe object of it
MAR.bacteria<-count(fullbiblio_dupsremoved$MABres)
#count number of AMS-testing studies on fungi in same way
MAR.fungi<-count(fullbiblio_dupsremoved$MAFres)
#join dataframes
MAR.both<-rbind(MAR.bacteria,MAR.fungi)
#aggregate the information
MAR.both<-aggregate(MAR.both$freq, list(MAR.both$x), sum)
#calculate how many studies equal 'yes' or 'not possible' i.e. are culture-based. Well that was simple wasn't it?
sum(MAR.both$x[MAR.both$Group.1%in%c('YES','NOT POSSIBLE')])
#calculate how many studies are culture-independent i.e. == gene, metagenome
count(fullbiblio_dupsremoved$wasMARtested)

#section 2.1.3 - Antimicrobials covered
#count how many studies doing antibiotic susceptibility testing
sum(MAR.bacteria$freq[MAR.bacteria$x%in%c('YES','NOT POSSIBLE')])
#count how many studies doing antifungal susceptibility testing
sum(MAR.fungi$freq[MAR.bacteria$x%in%c('YES','NOT POSSIBLE')])

MAR.bacteria_gene<-count(fullbiblio_dupsremoved$MABgeneres)
MAR.fungi_gene<-count(fullbiblio_dupsremoved$MAFgeneres)
#count how many studies doing antibiotic susceptibility testing
sum(MAR.bacteria_gene$freq[MAR.bacteria_gene$x%in%c('YES','NOT POSSIBLE')])
#count how many studies doing antifungal susceptibility testing
sum(MAR.fungi_gene$freq[MAR.bacteria_gene$x%in%c('YES','NOT POSSIBLE')])


# GENERAL TRENDS (24 AMS TESTING STUDIES) ---------------------------------

#this section of the code calculates summary statistics for the '2.1: General Trends' section of the report, summarising coverage of different topics etc. in the literature

#count how many antibiotics covered
ncol(MAR.allstudies.ABcols)
nlevels(ABclass1)
nlevels(ABclass2)
nlevels(ABclass3)
nlevels(ABclass4)

#create dataframes aggregated per class of antibiotic in case need to calculate resistance at lower classification resolution
MAR.allstudies.ABcols_class1<-MAR.allstudies.ABcols
colnames(MAR.allstudies.ABcols_class1)<-ABclass1
MAR.allstudies.ABcols_class1<-as.data.frame(sapply(unique(colnames(MAR.allstudies.ABcols_class1)), function(x) rowSums( MAR.allstudies.ABcols_class1[ , grep(x, names(MAR.allstudies.ABcols_class1)), drop=FALSE], na.rm = F)))
#recalculate B column at new class level
MAR.allstudies.ABcols_class1_B<-rowSums(apply(MAR.allstudies.ABcols_class1, 2, function(x) !is.na(x)))

MAR.allstudies.ABcols_class2<-MAR.allstudies.ABcols
colnames(MAR.allstudies.ABcols_class2)<-ABclass2
MAR.allstudies.ABcols_class2<-as.data.frame(sapply(unique(colnames(MAR.allstudies.ABcols_class2)), function(x) rowSums( MAR.allstudies.ABcols_class2[ , grep(x, names(MAR.allstudies.ABcols_class2)), drop=FALSE], na.rm = F)))
#recalculate B column at new class level
MAR.allstudies.ABcols_class2_B<-rowSums(apply(MAR.allstudies.ABcols_class2, 2, function(x) !is.na(x)))

MAR.allstudies.ABcols_class3<-MAR.allstudies.ABcols
colnames(MAR.allstudies.ABcols_class3)<-ABclass3
MAR.allstudies.ABcols_class3<-as.data.frame(sapply(unique(colnames(MAR.allstudies.ABcols_class3)), function(x) rowSums( MAR.allstudies.ABcols_class3[ , grep(x, names(MAR.allstudies.ABcols_class3)), drop=FALSE], na.rm = F)))
#recalculate B column at new class level
MAR.allstudies.ABcols_class3_B<-rowSums(apply(MAR.allstudies.ABcols_class3, 2, function(x) !is.na(x)))

MAR.allstudies.ABcols_class4<-MAR.allstudies.ABcols
colnames(MAR.allstudies.ABcols_class4)<-ABclass4
MAR.allstudies.ABcols_class4<-as.data.frame(sapply(unique(colnames(MAR.allstudies.ABcols_class4)), function(x) rowSums( MAR.allstudies.ABcols_class4[ , grep(x, names(MAR.allstudies.ABcols_class4)), drop=FALSE], na.rm = F)))
#recalculate B column at new class level
MAR.allstudies.ABcols_class4_B<-rowSums(apply(MAR.allstudies.ABcols_class4, 2, function(x) !is.na(x)))


#count those not equal to NA i.e. most studied for each AB column
sort(apply(!is.na(MAR.allstudies.ABcols), 2, function(x) sum(x)))
sort(apply(!is.na(MAR.allstudies.ABcols_class1), 2, function(x) sum(x)))
sort(apply(!is.na(MAR.allstudies.ABcols_class2), 2, function(x) sum(x)))
sort(apply(!is.na(MAR.allstudies.ABcols_class3), 2, function(x) sum(x)))
sort(apply(!is.na(MAR.allstudies.ABcols_class4), 2, function(x) sum(x)))

# Calculate MAR indices for later ----------------------------------------------------------------

dev.off()
#check if type column is a factor
levels(MAR.allstudies$type)
#coerce it to be a factor, change order of levels to increasing MAR order
MAR.allstudies$type<-factor(MAR.allstudies$type, levels=c("horse","urban","home","wwtp","poultry","cow","pig"))
#coerce distance.m column, N,B and AMS testing results columns to numeric - introducing NAs warning is fine (we need the NAs)
MAR.allstudies[,c(7,15:ncol(MAR.allstudies))]<-apply(MAR.allstudies[,c(7,15:ncol(MAR.allstudies))], 2, function(x) as.numeric(x))
#calculate MAR index for each row
MAR.allstudies$MAR.index<-rowSums(MAR.allstudies.ABcols, na.rm=T)/(MAR.allstudies$N*MAR.allstudies$B)
#view it to check it - within 0-1?
MAR.allstudies$MAR.index

# Pig CAFO AMR prevalence calculations ------------------------------------

#create subsets of data for staph/MRSA and enterobacteria (i.e. not staph)
staph.subset<-MAR.allstudies[grep('staph',MAR.allstudies$strain.type),]
entero.subset<-MAR.allstudies[-grep('staph',MAR.allstudies$strain.type),]

#create pig subset
pig.staph.subset<-staph.subset[staph.subset$type%in%c('pig'),]

#calculate tetracycline resistance for staph (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
pig.staph.tetracycline.mean_uw<-mean(pig.staph.subset$tetracycline/pig.staph.subset$N,na.rm=T)
pig.staph.tetracycline.mean_w<-metaprop(pig.staph.subset$tetracycline,pig.staph.subset$N,prediction = T)
pig.staph.tetracycline.mean_uw
pig.staph.tetracycline.mean_w

#calculate erythromycin resistance for staph (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
pig.staph.erythromycin.mean_uw<-mean(pig.staph.subset$erythromycin/pig.staph.subset$N,na.rm=T)
pig.staph.erythromycin.mean_w<-metaprop(pig.staph.subset$erythromycin,pig.staph.subset$N,prediction = T)
pig.staph.erythromycin.mean_uw
pig.staph.erythromycin.mean_w

#calculate tetracycline resistance for enterobacteria (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
pig.entero.subset<-entero.subset[entero.subset$type%in%c('pig'),]
pig.entero.tetracycline.mean_uw<-mean(pig.entero.subset$tetracycline/pig.entero.subset$N,na.rm=T)
pig.entero.tetracycline.mean_w<-metaprop(pig.entero.subset$tetracycline,pig.entero.subset$N,prediction = T)
pig.entero.tetracycline.mean_uw
pig.entero.tetracycline.mean_w

#calculate erythromycin resistance for enterobacteria (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
pig.entero.erythromycin.mean_uw<-mean(entero.subset$erythromycin/entero.subset$N,na.rm=T)
pig.entero.erythromycin.mean_w<-metaprop(pig.entero.subset$erythromycin,pig.entero.subset$N,prediction = T)
pig.entero.erythromycin.mean_uw
pig.entero.erythromycin.mean_w

# Poultry CAFO AMR prevalence calculations ------------------------------------

#Only 1 study (so no meta-nalysis/weighted means) for staph, none for enterobacteria.
#Just took results straight from paper

#Cow CAFO AMR prevalence calculations ------------------------------------
#3 studies for staph, one for entero (so no meta-analysis/weighted means

#create cow subset
cow.staph.subset<-staph.subset[staph.subset$type%in%c('cow'),]

#calculate ampicillin resistance for staph(ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
cow.staph.ampicillin.mean_uw<-mean(cow.staph.subset$ampicillin/cow.staph.subset$N,na.rm=T)
cow.staph.ampicillin.mean_w<-metaprop(cow.staph.subset$ampicillin,cow.staph.subset$N,prediction = T)
cow.staph.ampicillin.mean_uw
cow.staph.ampicillin.mean_w

#calculate penicillin resistance for staph(ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
cow.staph.penicillin.mean_uw<-mean(cow.staph.subset$penicillin/cow.staph.subset$N,na.rm=T)
cow.staph.penicillin.mean_w<-metaprop(cow.staph.subset$penicillin,cow.staph.subset$N,prediction = T)
cow.staph.penicillin.mean_uw
cow.staph.penicillin.mean_w

#calculate cefaclor resistance for staph(ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
cow.staph.cefaclor.mean_uw<-mean(cow.staph.subset$cefaclor/cow.staph.subset$N,na.rm=T)
cow.staph.cefaclor.mean_w<-metaprop(cow.staph.subset$cefaclor,cow.staph.subset$N,prediction = T)
cow.staph.cefaclor.mean_uw
cow.staph.cefaclor.mean_w

#WWTP AMR prevalence calculations ------------------------------------
#only 1 study for staph (so no meta-aalysis/weighted means), 3 for entero.

#create wttp subset
wwtp.entero.subset<-entero.subset[entero.subset$type%in%c('wwtp'),]

#calculate ceftazidime resistance for entero (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
wwtp.entero.ceftazidime.mean_uw<-mean(wwtp.entero.subset$ceftazidime/wwtp.entero.subset$N,na.rm=T)
wwtp.entero.ceftazidime.mean_w<-metaprop(wwtp.entero.subset$ceftazidime,wwtp.entero.subset$N,prediction = T)
wwtp.entero.ceftazidime.mean_uw
wwtp.entero.ceftazidime.mean_w

#calculate cefotaxime resistance for entero (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
wwtp.entero.cefotaxime.mean_uw<-mean(wwtp.entero.subset$cefotaxime/wwtp.entero.subset$N,na.rm=T)
wwtp.entero.cefotaxime.mean_w<-metaprop(wwtp.entero.subset$cefotaxime,wwtp.entero.subset$N,prediction = T)
wwtp.entero.cefotaxime.mean_uw
wwtp.entero.cefotaxime.mean_w

#Urban AMR prevalence calculations ------------------------------------
#4 staph studies, 0 entero (no meta/weighted means)

#create urban subset
urban.staph.subset<-staph.subset[staph.subset$type%in%c('urban'),]

#calculate tetracycline resistance for staph (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
urban.staph.tetracycline.mean_uw<-mean(urban.staph.subset$tetracycline/urban.staph.subset$N,na.rm=T)
urban.staph.tetracycline.mean_w<-metaprop(urban.staph.subset$tetracycline,urban.staph.subset$N,prediction = T)
urban.staph.tetracycline.mean_uw
urban.staph.tetracycline.mean_w

#calculate erythromycin resistance for staph (ignoring NAs) for pig CAFO studies, using unweighted mean and weighted means from metaprop (see proportion column, random effects model row)
urban.staph.erythromycin.mean_uw<-mean(urban.staph.subset$erythromycin/urban.staph.subset$N,na.rm=T)
urban.staph.erythromycin.mean_w<-metaprop(urban.staph.subset$erythromycin,urban.staph.subset$N,prediction = T)
urban.staph.erythromycin.mean_uw
urban.staph.erythromycin.mean_w

# SUPPLEMENTARY PLOTS OF RAW DATA -------------------------------------------------------------------

#display rick and morty colour palette
show_col(pal_rickandmorty()(12))
#pick 7 colours to colour different types of environment
type.cols<-pal_rickandmorty()(12)[c(9,4,5,6,8,2,7)]

#create a vector of the colour for each row of dataframe, according to type/environmental source of air
type.cols.perrow<-MAR.allstudies$type
levels(type.cols.perrow)<-pal_rickandmorty()(12)[c(9,4,5,6,8,2,7)]
type.cols.perrow<-as.character(type.cols.perrow)

#create a vector of the pch for each row of dataframe, according to study from which it came
study.pchs<-c(1:24)
study.pchs.perrow<-as.numeric(as.factor(MAR.allstudies$study))

#Appendix Figure 1 (vertical orientation): Plot showing means per type/environmental source of air with data plotted over, with point (pch) type according to study
dev.off()
plot.new()
par(mar=c(5, 10, 5, 0))
layout(matrix(c(1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE))
#layout.show(4)
boxplot(MAR.allstudies$MAR.index~MAR.allstudies$type, outpch=NA, horizontal=T, yaxt='n',col="grey", 
        MAR.allstudies$type, range=1, staplewex=0.1)
points(MAR.allstudies$MAR.index,jitter(as.numeric(as.factor(MAR.allstudies$type)),2),
       cex = 1.5, 
       #cex = MAR.allstudies$B*0.2, 
       pch = study.pchs.perrow,
       col = type.cols.perrow)
mtext(side=2,"Type of air (environment)", line=8)
mtext(side=2, cex=0.75, line=1,las=2, c("Horse Stable", "Urban","Indoor","WWTP","Poultry CAFO","Cow CAFO","Pig CAFO"),at=1:7)
mtext(side=1,line=3,"MAR index")
plot.new()
par(mar=c(0, 0, 0, 0))
legend("left", bty='n',title=expression(bold('Study')), title.adj=0.01,cex=1.20,legend=levels(as.factor(MAR.allstudies$study)),col = tapply(type.cols.perrow,MAR.allstudies$study,sample, size=1), pch=1:24)

#Appendix Figure 1 (horizontal orientation): Plot showing means per type/environmental source of air with data plotted over, with point (pch) type according to study
dev.off()
plot.new()
par(mar=c(10, 7, 3, 0))
layout(matrix(c(1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE))
#layout.show(4)
boxplot(MAR.allstudies$MAR.index~MAR.allstudies$type, outpch=NA, horizontal=F, xaxt='n',col="grey", 
        MAR.allstudies$type, range=1, staplewex=0.1)
points(MAR.allstudies$MAR.index~jitter(as.numeric(as.factor(MAR.allstudies$type)),2),
       cex = 1.5, 
       #cex = MAR.allstudies$B*0.2, 
       pch = study.pchs.perrow,
       col = type.cols.perrow)
mtext(side=1,"Type of air (environment)", line=8)
mtext(side=1, cex=0.75, line=1,las=2, c("Horse Stable", "Urban","Indoor", "WWTP","Poultry CAFO","Cow CAFO","Pig CAFO"),at=1:7)
mtext(side=2,line=3,"MAR index")
plot.new()
par(mar=c(0, 0, 0, 0))
legend("left", bty='n',title=expression(bold('Study')), title.adj=0.01,cex=0.95,legend=levels(as.factor(MAR.allstudies$study)),col = tapply(type.cols.perrow,MAR.allstudies$study,sample, size=1), pch=1:24)

#preparation for making Appendix Figure 2
#add colour and pch vectors to AMS dataframe and make new dataframe
MAR.allstudies<-cbind(MAR.allstudies,type.cols.perrow,study.pchs.perrow)
#subset out CAFO and WWTP data for which we have distance transects to make Appendix Figure 2
CAFOWWTP<-MAR.allstudies[MAR.allstudies$type%in%c('pig','poultry','cow','wwtp'),]
#check what distances are covered more than one datapoint/row of df
tapply(CAFOWWTP$distance.m,CAFOWWTP$type,paste)
#subset out only these distances
CAFOWWTP<-CAFOWWTP[CAFOWWTP$distance.m%in%c(-100,-25,0,25,100,150),]
#remove rows only testing single strains as these are targetted studies and do not numbers of organisms isolated at different distances
CAFOWWTP<-CAFOWWTP[CAFOWWTP$N!=1,]
#aggregate to calculate mean number of orgs per distance
CAFOWWTP.mean<-aggregate(CAFOWWTP$N,list(CAFOWWTP$distance.m),mean)
#aggregate to calculate standard error number of orgs per distance
CAFOWWTP.se<-aggregate(CAFOWWTP$N,list(CAFOWWTP$distance.m),se)
#rename columns
colnames(CAFOWWTP.mean)<-c('distance.m','N')
colnames(CAFOWWTP.se)<-c('distance.m','N')

#Appendix Figure 2: Plot showing mrelationship of number of organisms isolated to distance from point sources
dev.off()
plot.new()
#par(mar=c(5, 10, 5, 0))
layout(matrix(c(1, 1, 1, 1, 2, 2), nrow=1, byrow=TRUE))
#layout.show(4)
plot(CAFOWWTP.mean$N~as.numeric(CAFOWWTP.mean$distance.m), pch=19, cex=2, ylim=c(0,55), 
     col=c("#82491EFF","#E762D7FF","black","#E762D7FF","#E762D7FF","#E762D7FF"),
     xlab="Downwind distance from point-source(m)",
     ylab="Number of organisms isolated")
arrows(x0=as.numeric(CAFOWWTP.mean$distance.m), y0=CAFOWWTP.mean$N-CAFOWWTP.se$N,x1=as.numeric(CAFOWWTP.mean$distance.m), y1=CAFOWWTP.mean$N+CAFOWWTP.se$N,length=0.0,lwd=5,col=c("#82491EFF","#E762D7FF","black","#E762D7FF","#E762D7FF","#E762D7FF"))
points(CAFOWWTP$N~jitter(as.numeric(CAFOWWTP$distance.m)), cex=1.5, lwd=2, col=alpha(as.character(CAFOWWTP$type.cols.perrow),0.75), pch=as.numeric(CAFOWWTP$study.pchs.perrow))
#see studies
MAR.allstudies$study[as.numeric(rownames(CAFOWWTP))]
plot.new()
par(mar=c(0, 0, 0, 0))
legend("left", bty='n',title=expression(bold('Study')), title.adj=0.01,cex=1.5,legend=levels(as.factor(CAFOWWTP$study)),col = tapply(as.character(CAFOWWTP$type.cols.perrow),CAFOWWTP$study,sample, size=1), 
       pch=tapply(CAFOWWTP$study.pchs.perrow,CAFOWWTP$study,sample, size=1))

# META-ANALYSIS OF PROPORTIONS ---------------------------------------------------------

#multiply N x B in order to get maximum possible number of resistants, calculate how many resistants were observed in total
#aggregate this info by study (and output the type column also), calculating the mean number of observed and the mean number of total possible per study
#this is bit confusing, but essentially allows the meta-analysis of proportions to compare the MAR index between studies
perstudy.df<-aggregate(cbind(MAR.allstudies$N*MAR.allstudies$B,rowSums(MAR.allstudies.ABcols,na.rm = T)),
                       list(MAR.allstudies$study,MAR.allstudies$type),
                       mean)

#check dataframe
perstudy.df
#rename columns
colnames(perstudy.df)<-c('study','type','N','mean.res')

#create new colour vector equivelent to appendix figures by matching the order in ther perstudy.df dataframe to the order of the levels in the MAR.allstudies dataframe, 
#and then sampling 1 colour from each study
colours.match<-match(perstudy.df$study,names(tapply(type.cols.perrow,MAR.allstudies$study,sample, size=1)))
colours.type<-tapply(type.cols.perrow,MAR.allstudies$study,sample, size=1)[colours.match]

#perform meta-analysis of proportions, on rounded mean observations of resistance (n/observations) vs rounded mean number of max possible observations (events)
#again, to reiterate, this essentially performs a meta-analysis of MAR indices
perstudy.meta<-metaprop(round(perstudy.df$mean.res),round(perstudy.df$N),perstudy.df$study,byvar=perstudy.df$type)
perstudy.meta

#Figure 2_1
dev.off()
plot.new()
pdf(file = "Figure2_1.pdf", width = 10, height = 15) 
forest(perstudy.meta,
       xlim = c(0,1),
       xlab = 'MAR Index (weighted mean)',
       overall = F,
       col.study = colours.type,
       col.square  = colours.type)
dev.off()
#print summary stats for putting into report (Test for subgroup differences (random effects model):)
summary(perstudy.meta)

# Appendix B: SOME EXTRA ANALYSIS ON TO VALIDATE META-ANALYSIS ----------------------------------------------------------------

#create dataframe of potential explanatory variables
explans<-cbind(as.factor(MAR.allstudies$type),as.numeric(MAR.allstudies$distance.m),as.factor(MAR.allstudies$air.method), as.factor(MAR.allstudies$strain.type), as.numeric(MAR.allstudies$N),as.numeric(MAR.allstudies$B))
colnames(explans)<-c('type','distance','airmethod','strain','number of strains tested (N)','number of antibiotics tested (B)')

#generalised mixed binomial model to check for other effects that may be overlapping with type/environmental source of air
library(lme4)
model<-glmer(MAR.allstudies$MAR.index~explans[,1]+(1|explans[,2])+(1|explans[,3])+(1|explans[,4])+(1|explans[,5])+(1|explans[,6]), family=binomial)
options(scipen=1000)
summary(model)
plot(model)

#random forest to check for other effects that may be overlapping with type/environmental source of air

set.seed(1234)
ranf<-randomForest(explans,MAR.allstudies$MAR.index, ntree=10000, importance = T)
dev.off()
plot(ranf) #check error stabilised
ranf #summary of random forest results
varImp(ranf) #variable importances
varImpPlot(ranf, type=1) #variable importance plot

#regression tree to check for any effect of the identify of antibiotics tested (glmer and randomForest can't deal with NAs), that may be overlapping with type/environmental source of air
set.seed(1234)
tree<-rpart(MAR.allstudies$MAR.index~.,cbind(MAR.allstudies$type,MAR.allstudies[,15:73]))
plot(tree,margin = 0.10)
text(tree)
summary(tree)

# END ---------------------------------------------------------------------

