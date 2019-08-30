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
conc_tested<-NROW(unique(sensInfo$Dose.uM))
#sort these doses from lowest->highest (optional; idk this just looks better)
sortedUniqueDose <- as.character(sort.int(unique(sensInfo$Dose.uM),decreasing=FALSE))

#Create temportary data.frames to use for processing
#Dose Info from sensInfo & reformatting
allDose <- subset(phenoData, select=c(UID, concentration))
rownames(allDose)<-allDose$UID
allDose<-subset(allDose,select=-c(UID))
#Viability Info (& dose info) from phenoData (DNA) & reformatting
allViability <- subset(phenoData, select=c(UID, concentration, DNA...))
rownames(allViability)<-allViability$UID
allViability<-subset(allViability,select=-c(UID))

#Create a temporary data.frame to store sorted UID's
#sortedPhenoUID<-subset(phenoData, select=c(UID,samplename))
#sortedPhenoUID2<-arrange(sortedPhenoUID,samplename)

#Create empty dose array & format
doseArray <- array(NA, dim = c(NROW(phenoData$UID),conc_tested))
rownames(doseArray)<-phenoData$UID
colnames(doseArray)<-sortedUniqueDose #for now, colnames will be the dosages instead of "doses.1", "doses.2", etc.
#Create empty viability array & format
viabilityArray <- array(NA, dim=c(NROW(phenoData$UID),conc_tested))
rownames(viabilityArray)<-phenoData$UID
colnames(viabilityArray)<-sortedUniqueDose

#for loops
for (i in 1:NROW(allDose)){ #for loop iterates as long as the number of rows (#UID's)
  #go down the row of doses by iteration
  dose_row_name <- rownames(allDose)[i] #dose_row_name is the UID at that iteration
  dosage <- allDose[i,] #dosage is the dose (concentration) at that iteration
  doseArray[dose_row_name,as.character(dosage)]<-dosage #index into the correct spot using the above variables, and assign 'dosage' there
}
for (i in 1:NROW(allViability)){ #for loop iterates as long as the number of rows (#UID's)
  dose_row_name <- rownames(allViability)[i] #UID
  dosage <- allViability[i,1] #concentration
  viability<-allViability[i,2] #viability is on the same row, one column over
  viabilityArray[dose_row_name,as.character(dosage)]<-viability #index into the correct spot using the above variables, and assign 'viability' there
}

#Rename columns
colnames(doseArray)<-c("doses1","doses2","doses3","doses4","doses5","doses6","doses7","doses8","doses9","doses10","doses11","doses12","doses13","doses14","doses15","doses16","doses17","doses18","doses19","doses20","doses21","doses22","doses23","doses24","doses25","doses26","doses27","doses28","doses29","doses30","doses31","doses32","doses33","doses34","doses35","doses36","doses37","doses38","doses39","doses40","doses41","doses42","doses43","doses44","doses45","doses46","doses47","doses48","doses49","doses50","doses51","doses52","doses53","doses54","doses55","doses56","doses57","doses58","doses59","doses60","doses61","doses62","doses63","doses64","doses65","doses66","doses67","doses68","doses69","doses70","doses71","doses72","doses73","doses74","doses75","doses76","doses77","doses78","doses79","doses80","doses81","doses82","doses83","doses84","doses85","doses86","doses87","doses88","doses89","doses90","doses91","doses92","doses93","doses94","doses95","doses96","doses97","doses98","doses99","doses100","doses101","doses102","doses103","doses104","doses105","doses106","doses107","doses108")
colnames(viabilityArray)<-c("doses1","doses2","doses3","doses4","doses5","doses6","doses7","doses8","doses9","doses10","doses11","doses12","doses13","doses14","doses15","doses16","doses17","doses18","doses19","doses20","doses21","doses22","doses23","doses24","doses25","doses26","doses27","doses28","doses29","doses30","doses31","doses32","doses33","doses34","doses35","doses36","doses37","doses38","doses39","doses40","doses41","doses42","doses43","doses44","doses45","doses46","doses47","doses48","doses49","doses50","doses51","doses52","doses53","doses54","doses55","doses56","doses57","doses58","doses59","doses60","doses61","doses62","doses63","doses64","doses65","doses66","doses67","doses68","doses69","doses70","doses71","doses72","doses73","doses74","doses75","doses76","doses77","doses78","doses79","doses80","doses81","doses82","doses83","doses84","doses85","doses86","doses87","doses88","doses89","doses90","doses91","doses92","doses93","doses94","doses95","doses96","doses97","doses98","doses99","doses100","doses101","doses102","doses103","doses104","doses105","doses106","doses107","doses108")

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
colnames(doseFrame)<-c("dose1.uM","dose2.uM","dose3.uM","dose4.uM","dose5.uM","dose6.uM","dose7.uM","dose8.uM","dose9.uM","dose10.uM","dose11.uM","dose12.uM","dose13.uM","dose14.uM","dose15.uM","dose16.uM","dose17.uM","dose18.uM","dose19.uM","dose20.uM","dose21.uM","dose22.uM","dose23.uM","dose24.uM","dose25.uM","dose26.uM","dose27.uM","dose28.uM","dose29.uM","dose30.uM","dose31.uM","dose32.uM","dose33.uM","dose34.uM","dose35.uM","dose36.uM","dose37.uM","dose38.uM","dose39.uM","dose40.uM","dose41.uM","dose42.uM","dose43.uM","dose44.uM","dose45.uM","dose46.uM","dose47.uM","dose48.uM","dose49.uM","dose50.uM","dose51.uM","dose52.uM","dose53.uM","dose54.uM","dose55.uM","dose56.uM","dose57.uM","dose58.uM","dose59.uM","dose60.uM","dose61.uM","dose62.uM","dose63.uM","dose64.uM","dose65.uM","dose66.uM","dose67.uM","dose68.uM","dose69.uM","dose70.uM","dose71.uM","dose72.uM","dose73.uM","dose74.uM","dose75.uM","dose76.uM","dose77.uM","dose78.uM","dose79.uM","dose80.uM","dose81.uM","dose82.uM","dose83.uM","dose84.uM","dose85.uM","dose86.uM","dose87.uM","dose88.uM","dose89.uM","dose90.uM","dose91.uM","dose92.uM","dose93.uM","dose94.uM","dose95.uM","dose96.uM","dose97.uM","dose98.uM","dose99.uM","dose100.uM","dose101.uM","dose102.uM","dose103.uM","dose104.uM","dose105.uM","dose106.uM","dose107.uM","dose108.uM")
#Delete dose column
sensInfo<-subset(sensInfo,select=-c(Dose.uM))
sensInfo<-merge(sensInfo,doseFrame,by=0)
#Final formatting & arrange rows in ascending order
sensInfo <- arrange(sensInfo, cellid)
rownames(sensInfo)<-sensInfo$Row.names
sensInfo<-subset(sensInfo,select=-c(Row.names))

################################################################################################################

#Write objects into data/ file
saveRDS(sensInfo,"data/sensitivityInfo.rds")
saveRDS(sensProf,"data/sensitivityProfiles_biorep.rds")
saveRDS(sensRaw,"data/sensitivityRaw.rds")

#phenoData : phenoData from phenowLot.csv
#UID : data.frame of unique ID's (2,573 entries)
#sensitivity$info : sensInfo (2,573 entries)
#sensitivity$profiles : sensProf (2,573 entries)
#sensitivity$raw : sensRaw (2,573 entries)
