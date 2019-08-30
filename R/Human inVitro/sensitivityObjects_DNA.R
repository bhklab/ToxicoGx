library(dplyr)
library(PharmacoGx)

phenoData <- readRDS("data/phenoData.rds")

################################################################################################################

## (START) CREATING SENSITIVITY$INFO ##
sensInfo <- subset(phenoData,select=c(UID,cellid, drugid, duration, individual_id, concentration))
#set rownames to UIDs
rownames(sensInfo) <- sensInfo$UID
sensInfo <- subset(sensInfo,select=-c(UID))
#rename columns
names(sensInfo)<-c("cellid","drugid","duration_h","replicate", "Dose.uM")

################################################################################################################

## CREATE SENSITIVITY$PROFILES (DNA only) ##
#Select only for entries where DNA info is present -> DNApresent
#For DNA-based viability info, all dosages were in uM
DNApresent<-subset(phenoData,phenoData$DNA... != "NA",select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, DNA...,UID))
#The data.frame finalCombined will start with all the samples that do not have DNA data
finalCombined<-subset(phenoData,is.na(phenoData$DNA...),select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, DNA...,UID))
finalCombined$slope_recomputed<-NA
finalCombined$auc_recomputed<-NA
#Get all unique drug names in DNApresent (for the for loop later)
uniqueDrugNames<-unique(DNApresent$drugid)

#Create an empty data frame with all required columns
emptyFrame<-subset(phenoData,select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, DNA...,UID))
emptyFrame$slope_recomputed<-NA
emptyFrame$auc_recomputed <- NA
emptyFrame<-emptyFrame[FALSE,]

#The data.frame combinedSamples will start empty, and accumulate data in the for loop
combinedSamples<-emptyFrame

#for loop:
for (drug in uniqueDrugNames){ #iterate through each unique drug (set of experiments)
  combinedSamples<-emptyFrame #reset combinedSamples at the beginning of each loop
  drugSamples<-subset(DNApresent,DNApresent$drugid == drug) #drugSamples contains all the experiments for one drug
  
  controls<-subset(drugSamples,drugSamples$dose_level == "Control") #controls contains all the control experiments for one drug
  controls$slope_recomputed <-NA #give controls an empty SLOPE column of NA's
  controls$auc_recomputed <-NA
  finalCombined<-rbind(finalCombined,controls) #rbind controls to master profiles list (control experiments are accounted for; DONE)
  
  drugSamples<-subset(drugSamples,drugSamples$dose_level != "Control") #subset control experiments out of the experiments for the drug
  drugSamples$slope_recomputed <-NA #give drugSamples an empty SLOPE column of NA's
  drugSamples$auc_recomputed <- NA
  
  #divide drugSamples by time -> 2, 8, 24h and replicate number -> 1, 2
  samples2_1<-subset(drugSamples,drugSamples$duration == "2 hr" & drugSamples$individual_id == 1)
  samples2_2<-subset(drugSamples,drugSamples$duration == "2 hr" & drugSamples$individual_id == 2)
  samples8_1<-subset(drugSamples,drugSamples$duration == "8 hr" & drugSamples$individual_id == 1)
  samples8_2<-subset(drugSamples,drugSamples$duration == "8 hr" & drugSamples$individual_id == 2)
  samples24_1<-subset(drugSamples,drugSamples$duration == "24 hr" & drugSamples$individual_id == 1)
  samples24_2<-subset(drugSamples,drugSamples$duration == "24 hr" & drugSamples$individual_id == 2)
  
  #for each set of time-divided & duplicate-divided samples, check if at least 3 dose levels are present;
  #YES -> computeSlope & computeAUC
  #NO -> nothing (leave as NA)
  if (NROW(samples2_1) >= 3){
    samples2_1$slope_recomputed = computeSlope(samples2_1$concentration,samples2_1$DNA...)
    samples2_1$auc_recomputed = computeAUC(samples2_1$concentration, samples2_1$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  if (NROW(samples2_2) >= 3){
    samples2_2$slope_recomputed = computeSlope(samples2_2$concentration,samples2_2$DNA...)
    samples2_2$auc_recomputed = computeAUC(samples2_2$concentration, samples2_2$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  if (NROW(samples8_1) >= 3){
    samples8_1$slope_recomputed = computeSlope(samples8_1$concentration,samples8_1$DNA...)
    samples8_1$auc_recomputed = computeAUC(samples8_1$concentration, samples8_1$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  if (NROW(samples8_2) >= 3){
    samples8_2$slope_recomputed = computeSlope(samples8_2$concentration,samples8_2$DNA...)
    samples8_2$auc_recomputed = computeAUC(samples8_2$concentration, samples8_2$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  if (NROW(samples24_1) >= 3){
    samples24_1$slope_recomputed = computeSlope(samples24_1$concentration,samples24_1$DNA...)
    samples24_1$auc_recomputed = computeAUC(samples24_1$concentration, samples24_1$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  if (NROW(samples24_2) >= 3){
    samples24_2$slope_recomputed = computeSlope(samples24_2$concentration,samples24_2$DNA...)
    samples24_2$auc_recomputed = computeAUC(samples24_2$concentration, samples24_2$DNA...,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
  }
  
  #combinedSamples is an rbind of each time-divided sample after above processing
  #combinedSamples is reset to an empty data.frame each iteration
  combinedSamples<-rbind(samples2_1,samples2_2)
  combinedSamples<-rbind(combinedSamples,samples8_1)
  combinedSamples<-rbind(combinedSamples,samples8_2)
  combinedSamples<-rbind(combinedSamples,samples24_1)
  combinedSamples<-rbind(combinedSamples,samples24_2)
  
  #finalCombined is a master profiles list - an rbinded accumulation  of all combinedSamples data.frames
  finalCombined<-rbind(finalCombined,combinedSamples)
}

#Final formatting & arrange rows in ascending order
sensProf<-subset(finalCombined, select=c(UID,slope_recomputed,auc_recomputed,samplename))
sensProf <- arrange(sensProf, samplename)
rownames(sensProf)<-sensProf$UID
sensProf<-subset(sensProf,select=-c(samplename,UID))
sensProf<-as.matrix(sensProf)

################################################################################################################

## CREATE SENSITIVITY$RAW (DNA only) ##
## MUST RUN (START) SENSITIVITY$INFO BEFORE THIS ##
#number of unique concentrations = number of columns; set this number as conc_tested
conc_tested<-c("Control","Low","Middle","High")

#Create temportary data.frames to use for processing
#Dose Info from sensInfo & reformatting
allDose <- subset(phenoData, select=c(UID, concentration,dose_level))
rownames(allDose)<-allDose$UID
allDose<-subset(allDose,select=-c(UID))
#Viability Info (& dose info) from phenoData (DNA) & reformatting
allViability <- subset(phenoData, select=c(UID, concentration, DNA...,dose_level))
rownames(allViability)<-allViability$UID
allViability<-subset(allViability,select=-c(UID))

#Create empty dose array & format
doseArray <- array(NA, dim = c(NROW(phenoData$UID),NROW(conc_tested)))
rownames(doseArray)<-phenoData$UID
colnames(doseArray)<-conc_tested #for now, colnames will be the dosages instead of "doses.1", "doses.2", etc.
#Create empty viability array & format
viabilityArray <- array(NA, dim=c(NROW(phenoData$UID),NROW(conc_tested)))
rownames(viabilityArray)<-phenoData$UID
colnames(viabilityArray)<-conc_tested

#for loops
for (i in 1:NROW(allDose)){ #for loop iterates as long as the number of rows (#UID's)
  #go down the row of doses by iteration
  dose_row_name <- rownames(allDose)[i] #dose_row_name is the UID at that iteration
  dosage <- allDose[i,1] #dosage is the dose (concentration) at that iteration
  dose_level <- allDose[i,2] #dose_level is the dose level (control, low, middle, high) at that iteration
  doseArray[dose_row_name,dose_level]<-dosage #index into the correct spot using the above variables, and assign 'dosage' there
}
for (i in 1:NROW(allViability)){ #for loop iterates as long as the number of rows (#UID's)
  dose_row_name <- rownames(allViability)[i] #UID
  dosage <- allViability[i,1] #concentration
  viability<-allViability[i,2] #viability is on the same row, one column over
  dose_level <- allViability[i,3]
  viabilityArray[dose_row_name,dose_level]<-viability #index into the correct spot using the above variables, and assign 'viability' there
}

#Rename columns
colnames(doseArray)<-c("doses1","doses2","doses3","doses4")
colnames(viabilityArray)<-c("doses1","doses2","doses3","doses4")

#Combine doseArray and viabilityArray -> 3D array
library(abind)
doseViabilityArray<-abind(doseArray, viabilityArray, along=3)
dimnames(doseViabilityArray)[[3]]<-c("Dose","Viability")

sensRaw <- doseViabilityArray

################################################################################################################

## FINISH SENSITIVITY$INFO ## -> need to be like sensitivity$raw's dose array
#Convert to data.frame
doseFrame <- as.data.frame(doseArray)
#doseFrame has different colnames
colnames(doseFrame)<-c("dose1.uM","dose2.uM","dose3.uM","dose4.uM")
#Delete dose column
sensInfo<-subset(sensInfo,select=-c(Dose.uM))
sensInfo<-merge(sensInfo,doseFrame,by=0)
#Final formatting & arrange rows in ascending order
sensInfo <- arrange(sensInfo, cellid)
rownames(sensInfo)<-sensInfo$Row.names
sensInfo<-subset(sensInfo,select=-c(Row.names))

################################################################################################################

#Write objects into data/ file
saveRDS(sensInfo,"data/sensitivityInfo_DNA.rds")
saveRDS(sensProf,"data/sensitivityProfiles_DNA.rds")
saveRDS(sensRaw,"data/sensitivityRaw_DNA.rds")

#phenoData : phenoData from phenowLot.csv
#UID : data.frame of unique ID's (2,573 entries)
#sensitivity$info : sensInfo (2,573 entries)
#sensitivity$profiles : sensProf (2,573 entries)
#sensitivity$raw : sensRaw (2,573 entries)