## Modularize functions for creating toxicoSet -> in vitro human, rat

# Function to map TG-GATEs compounds with Lab compound annotations in TG-GATEs phenoData
# OUTPUT: preformatted phenoData with correct drug mapping, number & order of columns for in vitro TG-GATEs tSets
create_phenoData <- function(species=c("Human","Rat")){
  library(dplyr)
  library(gdata)
  
  #load master phenoData file from TG-GATEs
  #Master phenoData file from TG-GATEs: #14 from https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html
  all_attribute <- read.delim("data/Open-tggates_AllAttribute.tsv",stringsAsFactors = F)
  
  #subset out unwanted samples (rows)
  all_attribute <- subset(all_attribute, all_attribute$BARCODE != "No ChipData" & 
                            all_attribute$TEST_TYPE == "in vitro" &
                            all_attribute$DOSE_UNIT != "ƒÊg/kg")
  
  #subset out unwanted columns (all NA; for in vivo studies)
  all_attribute <- Filter(function(x) !all(is.na(x)), all_attribute)
  
  
  ## Manual mapping for drug curation
  drug_curation <- unique(all_attribute[, "COMPOUND_NAME", drop=F]) #get all unique drugs from TG-GATEs in vitro experiments
  drug_curation$unique.drugid <- drug_curation$COMPOUND_NAME #add column for lab drugid mapping; start off with all unique.drugid == tggates.drugid's
  
  lab_curation <- read.csv("data/drugs_with_ids.csv", stringsAsFactors = F) #load lab annotation file
  
  #exact/case match with lab annotations
  drug_curation$unique.drugid <- sapply(drug_curation$COMPOUND_NAME, function(x) if (tolower(x) %in% tolower(lab_curation$unique.drugid)) (lab_curation[which(tolower(x) == tolower(lab_curation$unique.drugid)),"unique.drugid"]) else x)
  #synonyms - manually curate
  drug_curation[drug_curation$COMPOUND_NAME == "coumarin", "unique.drugid"] <- "2H-1-Benzopyran-2-one"
  drug_curation[drug_curation$COMPOUND_NAME == "cycloheximide", "unique.drugid"] <- "cicloheximide"
  drug_curation[drug_curation$COMPOUND_NAME == "cyclosporine A", "unique.drugid"] <- "ciclosporin"
  drug_curation[drug_curation$COMPOUND_NAME == "ibuprofen", "unique.drugid"] <- "2-(4-Isobutylphenyl)propanoic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "indomethacin", "unique.drugid"] <- "indometacin"
  drug_curation[drug_curation$COMPOUND_NAME == "methimazole", "unique.drugid"] <- "thiamazole"
  drug_curation[drug_curation$COMPOUND_NAME == "nitrofurazone", "unique.drugid"] <- "nitrofural"
  drug_curation[drug_curation$COMPOUND_NAME == "phenylanthranilic acid", "unique.drugid"] <- "N-phenylanthranilic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "WY-14643", "unique.drugid"] <- "pirinixic acid"
  drug_curation[drug_curation$COMPOUND_NAME == "cefalotin", "unique.drugid"] <- "cephalothin"
  drug_curation[drug_curation$COMPOUND_NAME == "valproic acid", "unique.drugid"] <- "Valproic.acid"
  drug_curation[drug_curation$COMPOUND_NAME == "chlorpheniramine", "unique.drugid"] <- "chlorphenamine"
  drug_curation[drug_curation$COMPOUND_NAME == "cephalothin", "unique.drugid"] <- "cefalotin"
  
  drug_curation <- drug_curation[,c(2,1)] #reorder columns
  names(drug_curation)[2] <- "tggates.drugid" #rename column
  rownames(drug_curation) <- drug_curation$unique.drugid #rename rows
  
  
  ################# BACK TO PHENODATA ####################
  ## Add necessary columns
  ##Batchid & conversions
  prl<-"C:/Strawberry/perl/bin/perl5.30.0.exe"
  #Human
  if (species == "Human"){
    batch <- read.xls("C:/Users/Owner/Desktop/ToxicoGx - New/data/nar-02356-data-e-2014-File006.xlsx", sheet = 4, header=FALSE, as.is=TRUE, perl=prl)
    batch <- batch[-c(1,2),]
    batch <- batch[,c(1,4)]
    names(batch) <- c("BARCODE", "CELL_NAME_TYPE_ID")
    batch <- subset(batch, batch$BARCODE != "No ChipData")
    all_attribute <- merge(all_attribute, batch, by.x = "BARCODE", by.y = "BARCODE")
    
    conv <- readRDS("data/conversions_human.rds")
  } else if (species == "Rat"){ #Rat
    batch <- read.xls("C:/Users/Owner/Desktop/ToxicoGx - New/data/nar-02356-data-e-2014-File006.xlsx", sheet = 5, header=FALSE, as.is=TRUE, perl=prl)
    batch <- batch[-c(1:4),]
    batch <- batch[,c(1,9)]
    colnames(batch) <- c("drugid", "batchid")
    colnames(batch) <- as.factor(colnames(batch))
    batch$drugid[batch$drugid=="TNFalpha"]<- "TNFƒ¿"
    
    all_attribute <- subset(all_attribute, all_attribute$ARR_DESIGN == "Rat230_2" & all_attribute$TEST_TYPE == "in vitro")
    all_attribute <- merge(all_attribute, batch, by.x = "COMPOUND_NAME", by.y = "drugid")
    
    conv <- readRDS("data/conversions_rat.rds")
  }
  #Include converted doses
  need_conversion <- subset(all_attribute,all_attribute$DOSE_UNIT != "ƒÊM",select = -c(DOSE,DOSE_UNIT))
  converted <- merge(need_conversion,conv,by.x="BARCODE",by.y="BARCODE")
  all_attribute <- subset(all_attribute,all_attribute$DOSE_UNIT == "ƒÊM")
  all_attribute$RECOMP_DOSE <- all_attribute$DOSE
  all_attribute <- rbind(all_attribute,converted)
  ##Merge drug mapping to lab annotations to phenoData
  all_attribute <- merge(all_attribute, drug_curation, by.x = "COMPOUND_NAME", by.y = "tggates.drugid")
  #CEL file name
  all_attribute$celfilename <- paste0(as.character(all_attribute$BARCODE),".CEL")
  # Label control or perturbation xptype
  all_attribute$xptype <- NA
  all_attribute$xptype[all_attribute$DOSE_LEVEL != "Control"] <- "perturbation"
  all_attribute$xptype[all_attribute$DOSE_LEVEL == "Control"] <- "control"
  #Substring barcode & assign as cellid
  all_attribute$BARCODE<-substr(all_attribute$BARCODE,3,12)
  all_attribute$cellid <- all_attribute$BARCODE
  #fix 
  all_attribute$unique.drugid[all_attribute$unique.drugid=="TNFƒ¿"]<- "TNFa"
  #take out " hr" from duration column entries
  all_attribute$SACRI_PERIOD <- gsub(" hr","",all_attribute$SACRI_PERIOD)
  #UID
  all_attribute$UID<-paste("drugid_",all_attribute$COMPOUND.Abbr.,"_",all_attribute$unique.drugid,"_",all_attribute$BARCODE,"_",all_attribute$DOSE,"_",all_attribute$SACRI_PERIOD,"hr_rep",all_attribute$INDIVIDUAL_ID, sep="")
  
  #reorder columns
  #all_attribute <- all_attribute[,c(1,2,3,4,5,21,16,17,22,18,15,26,23,8,9,19,20,27,10,11,12,6,7,24,25,13,14)]
  all_attribute <- all_attribute[,c(2,3,4,5,6,22,17,18,23,19,16,27,24,1,9,10,20,21,28,11,12,13,7,8,25,26,14,15)] 
  
  #Rename columns
  colnames(all_attribute) <- c("samplename","chiptype","exp_id","group_id","individual_id","batchid",
                               "concentration_old","concentration_old_units","concentration","dose_level",
                               "duration","cellid","drugid","tggates_drugid","drugid_abbr","drugid_no","DNA","LDH",
                               "UID","species","test_type","sex_type","organ_id","material_id",
                               "celfilename","xptype","STRAIN_TYPE","ADM_ROUTE_TYPE")
  
  #replace "ƒÊ" with "µ"
  all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊM"] <- "µM"
  all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊg/mL"] <- "µg/mL"
  
  #Rearrange order of rows by ascending order of barcodes
  all_attribute <- dplyr::arrange(all_attribute, samplename)
  #Assign rownames as samplename
  rownames(all_attribute)<-all_attribute$samplename
  
  return(all_attribute)
}

create_exprsData <- function(species=c("Human","Rat"), phenoData){
  library(affy)
  library(PharmacoGx)
  if (species == "Human"){
    # install.packages("data/hgu133plus2hsensgcdf_22.0.0.tar.gz", repos = NULL, type = "source")
    library("hgu133plus2hsensgcdf")
  } else if (species == "Rat"){
    # install.packages("data/rat2302rnensgcdf_23.0.0.tar.gz", repos = NULL, source = 'source')
    library("rat2302rnensgcdf")
  }
  celfn <- list.celfiles("data/CELfiles", full.names = TRUE)
  
  ###########################################  NORMALIZATION  ############################################
  # eset <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)
  ########################################################################################################
  
  ########################################  ALREADY NORMALIZED  ##########################################
  eset <- readRDS("C:/Users/Owner/Desktop/ToxicoGx - New/data/esetNORMALIZED_ONLY.rds")
  ########################################################################################################
  
  storageMode(eset)<-"environment"
  #rename eset@assayData column names
  colnames(eset@assayData$exprs)<-substr(colnames(eset@assayData$exprs),3,12)
  colnames(eset@assayData$se.exprs)<-substr(colnames(eset@assayData$se.exprs),3,12)
  rownames(eset@protocolData@data)<-substr(rownames(eset@protocolData@data),3,12)
  #subset eset@assayData : 2605 columns/samples -> 2573
  #missingCEL is a data.frame containing the info for 2573 samples
  missingCEL <- phenoData[,c("samplename"), drop=F]
  #subsetting samples
  eset@assayData$exprs<-subset(eset@assayData$exprs,select=missingCEL$samplename)
  eset@assayData$se.exprs<-subset(eset@assayData$se.exprs,select=missingCEL$samplename)
  eset@phenoData@data <- subset(eset@phenoData@data, rownames(eset@phenoData@data) %in% paste("00",missingCEL$samplename,".CEL", sep=""))
  eset@protocolData@data<-subset(eset@protocolData@data,rownames(eset@protocolData@data) %in% missingCEL$samplename)
  #subsetting probes
  eset<-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  #replace _at
  rownames(eset)<-gsub("_at","",rownames(eset))
  #lock eset@assayData environment again
  storageMode(eset)<-"lockedEnvironment"
  
  annotation(eset)<-"rna"
  
  return(eset)
}

create_featureData <- function(species=c("Human","Rat"), eset){
  library("xml2")
  library(biomaRt)
  
  if (species == "Human"){ ensembl_data <- "hsapiens_gene_ensembl" }
  else if (species == "Rat"){ ensembl_data <- "rnorvegicus_gene_ensembl" }
  CELgenes <- rownames(eset@assayData$exprs)
  
  ensembl<-useMart("ensembl", dataset = ensembl_data, host="uswest.ensembl.org",ensemblRedirect = FALSE)
  results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=CELgenes,mart=ensembl)
  uniqueBiomaRt<-results[!duplicated(results$ensembl_gene_id),]
  
  if(species == "Rat"){ return(uniqueBiomaRt) }
  
  labAnnot <- read.csv("data/annot_ensembl_all_genes.csv",header=TRUE)[,-1, drop=F] #read in lab's gene annotation file
  CELnotbiomaRt<-unique(CELgenes)[!(unique(CELgenes) %in% uniqueBiomaRt$ensembl_gene_id)] #in CELgenes but not in biomaRt output (66)
  newLabAnnot<-subset(labAnnot,labAnnot$gene_id %in% CELnotbiomaRt) #in CELgenes but not in biomaRt output but in lab annotation file (51)
  
  names(uniqueBiomaRt)<-c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id") #rename biomaRt output columns
  
  #Merge biomaRt output and lab annotations
  finalFeature<-rbind(uniqueBiomaRt,newLabAnnot)
  
  #For the remaining genes not mapped by biomaRt or in the lab annotation file:
  leftoverGenes<-as.data.frame(CELgenes[!(CELgenes %in% finalFeature$gene_id)]) #15 genes remaining
  names(leftoverGenes) <- "gene_id"
  leftoverGenes$gene_name <- NA
  leftoverGenes$gene_biotype <- NA
  leftoverGenes$EntrezGene.ID <- NA
  leftoverGenes$transcript_name <- NA
  leftoverGenes$transcript_id <- NA
  
  finalFeature<-rbind(finalFeature,leftoverGenes)
  
  names(finalFeature)[1] <- "Symbol"
  finalFeature$BEST <- NA
  finalFeature<-arrange(finalFeature,finalFeature$gene_id)
  rownames(finalFeature)<-finalFeature$gene_id
  
  return(finalFeature)
}

create_sensitivityProfiles <- function(phenoData){
  ## CREATE SENSITIVITY$PROFILES ##
  #Select only for entries where viability info is present -> viabpresent
  #For viab-based viability info, all dosages were in uM
  viabpresent<-subset(phenoData,phenoData$Viability != "NA",select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, Viability, UID))
  #The data.frame finalCombined will start with all the samples that do not have viab data
  finalCombined<-subset(phenoData,is.na(phenoData$Viability),select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, Viability, UID))
  finalCombined$slope_recomputed<-NA
  finalCombined$auc_recomputed<-NA
  #Get all unique drug names in viabpresent (for the for loop later)
  uniqueDrugNames<-unique(viabpresent$drugid)
  
  #Create an empty data frame with all required columns
  emptyFrame<-subset(phenoData,select=c(samplename, exp_id, group_id, individual_id, drugid, duration, concentration, dose_level, Viability,UID))
  emptyFrame$slope_recomputed<-NA
  emptyFrame$auc_recomputed <- NA
  emptyFrame<-emptyFrame[FALSE,]
  
  #The data.frame combinedSamples will start empty, and accumulate data in the for loop
  combinedSamples<-emptyFrame
  
  #for loop:
  for (drug in uniqueDrugNames){ #iterate through each unique drug (set of experiments)
    combinedSamples<-emptyFrame #reset combinedSamples at the beginning of each loop
    drugSamples<-subset(viabpresent,viabpresent$drugid == drug) #drugSamples contains all the experiments for one drug
    
    controls<-subset(drugSamples,drugSamples$dose_level == "Control") #controls contains all the control experiments for one drug
    controls$slope_recomputed <-NA #give controls an empty SLOPE column of NA's
    controls$auc_recomputed <-NA
    finalCombined<-rbind(finalCombined,controls) #rbind controls to master profiles list (control experiments are accounted for; DONE)
    
    drugSamples<-subset(drugSamples,drugSamples$dose_level != "Control") #subset control experiments out of the experiments for the drug
    drugSamples$slope_recomputed <-NA #give drugSamples an empty SLOPE column of NA's
    drugSamples$auc_recomputed <- NA
    
    #divide drugSamples by time -> 2, 8, 24h and replicate number -> 1, 2
    samples2_1<-subset(drugSamples,drugSamples$duration == "2" & drugSamples$individual_id == 1)
    samples2_2<-subset(drugSamples,drugSamples$duration == "2" & drugSamples$individual_id == 2)
    samples8_1<-subset(drugSamples,drugSamples$duration == "8" & drugSamples$individual_id == 1)
    samples8_2<-subset(drugSamples,drugSamples$duration == "8" & drugSamples$individual_id == 2)
    samples24_1<-subset(drugSamples,drugSamples$duration == "24" & drugSamples$individual_id == 1)
    samples24_2<-subset(drugSamples,drugSamples$duration == "24" & drugSamples$individual_id == 2)
    
    #for each set of time-divided & duplicate-divided samples, check if at least 3 dose levels are present;
    #YES -> computeSlope & computeAUC
    #NO -> nothing (leave as NA)
    if (NROW(samples2_1) >= 3){
      samples2_1$slope_recomputed = computeSlope(samples2_1$concentration,samples2_1$Viability)
      samples2_1$auc_recomputed = computeAUC(samples2_1$concentration, samples2_1$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
    if (NROW(samples2_2) >= 3){
      samples2_2$slope_recomputed = computeSlope(samples2_2$concentration,samples2_2$Viability)
      samples2_2$auc_recomputed = computeAUC(samples2_2$concentration, samples2_2$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
    if (NROW(samples8_1) >= 3){
      samples8_1$slope_recomputed = computeSlope(samples8_1$concentration,samples8_1$Viability)
      samples8_1$auc_recomputed = computeAUC(samples8_1$concentration, samples8_1$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
    if (NROW(samples8_2) >= 3){
      samples8_2$slope_recomputed = computeSlope(samples8_2$concentration,samples8_2$Viability)
      samples8_2$auc_recomputed = computeAUC(samples8_2$concentration, samples8_2$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
    if (NROW(samples24_1) >= 3){
      samples24_1$slope_recomputed = computeSlope(samples24_1$concentration,samples24_1$Viability)
      samples24_1$auc_recomputed = computeAUC(samples24_1$concentration, samples24_1$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
    if (NROW(samples24_2) >= 3){
      samples24_2$slope_recomputed = computeSlope(samples24_2$concentration,samples24_2$Viability)
      samples24_2$auc_recomputed = computeAUC(samples24_2$concentration, samples24_2$Viability,conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
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
  
  return (sensProf)
}

create_sensitivityRaw <- function(phenoData){
  #number of unique concentrations = number of columns; set this number as conc_tested
  conc_tested<-c("Control","Low","Middle","High")
  
  #Create temportary data.frames to use for processing
  #Dose Info from sensInfo & reformatting
  allDose <- subset(phenoData, select=c(UID, concentration,dose_level))
  rownames(allDose)<-allDose$UID
  allDose<-subset(allDose,select=-c(UID))
  #Viability Info (& dose info) from phenoData & reformatting
  allViability <- subset(phenoData, select=c(UID, concentration, Viability,dose_level))
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
  
  return(sensRaw)
}

create_sensitivityInfo <- function(phenoData, doseArray){
  sensInfo <- subset(phenoData,select=c(UID,cellid, drugid, duration, individual_id, concentration))
  #set rownames to UIDs
  rownames(sensInfo) <- sensInfo$UID
  sensInfo <- subset(sensInfo,select=-c(UID))
  #rename columns
  names(sensInfo)<-c("cellid","drugid","duration_h","replicate", "Dose.uM")
  
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
}

create_curationDrug <- function(phenoData){
  curationDrug <- unique(subset(phenoData, select=c(drugid, tggates_drugid)))
  rownames(curationDrug) <- curationDrug$drugid
  names(curationDrug) <- c("unique.drugid", "tggates.drugid")
  
  return(curationDrug)
}

create_curationCell <- function(phenoData){
  curationCell <- subset(phenoData, select=c(samplename, cellid))
  names(curationCell) <- c("unique.cellid", "tggates.cellid")
  
  return(curationCell)
}

create_curationTissue <- function(phenoData){
  curationTissue <- subset(phenoData, select=c(organ_id))
  curationTissue$organ_id <- "liver"
  curationTissue$tggates.tissueid <- "liver"
  names(curationTissue)[1] <- "unique.tissueid"
  
  return(curationTissue)
}

create_drug <- function(phenoData){
  drug <- unique(subset(phenoData, select=c(drugid, drugid_abbr, drugid_no)))
  rownames(drug) <- drug$drugid
  
  return(drug)
}

create_cell <- function(phenoData){
  cell <- subset(phenoData,select=c(samplename, organ_id, material_id, species, test_type, batchid))
  names(cell)<-c("cellid","tissueid","materialid", "species","testType","batchid")
  cell$tissueid<-"liver"
  
  return(cell)
}

getTGGATEs <- function(species=c("Human","Rat"),
                       type=c("DNA", "LDH")){
  library(dplyr)
  library(PharmacoGx)
  
  #get all elements of eset
  message("Creating phenoData object...")
  phenoData <- create_phenoData(species)
  message("phenoData object created!")
  message("Creating eset object...")
  eset <- create_exprsData(species, phenoData)
  message("eset object created!")
  message("Creating featureData object...")
  featureData <- create_featureData(species, eset)
  message("featureData object created!")
  
  message("Putting the eset together...")
  #put the eset together
  storageMode(eset)<-"environment"
  pData(eset) <- phenoData
  fData(eset) <- featureData
  storageMode(eset)<-"lockedEnvironment"
  message("Done!")
  
  message(paste("Requested ",type, sep = ""))
  if (type == "DNA"){ 
    phenoData <- subset(phenoData, select=-c(LDH)) 
  } else { 
    phenoData <- subset(phenoData, select=-c(DNA))
  }
  names(phenoData)[which(names(phenoData) == type)] <- "Viability"
  
  message("Creating Sensitivity Profiles...")
  sensitivityProfiles <- create_sensitivityProfiles(phenoData)
  message("Sensitivity Profiles object created!")
  message("Creating Sensitivity Raw...")
  sensitivityRaw <- create_sensitivityRaw(phenoData)
  message("Sensitivity Raw object created!")
  message("Creating Sensitivity Info...")
  sensitivityInfo <- create_sensitivityInfo(phenoData, doseArray=as.array(sensitivityRaw[,,1]))
  message("Sensitivity Info object created!")
  
  message("Creating curation objects...")
  curationDrug <- create_curationDrug(phenoData)
  curationCell <- create_curationCell(phenoData)
  curationTissue <- create_curationTissue(phenoData)
  message("Done!")
  message("Creating cell, drug objects...")
  drug <- create_drug(phenoData)
  cell <- create_cell(phenoData)
  message("Done!")
  
  message("Putting toxicoSet together...")
  TGGATES_human <- PharmacoSet(paste("TGGATES_human",type,sep=""),
                                  molecularProfiles=list("rna"=eset),
                                  cell=cell,
                                  drug=drug,
                                  sensitivityInfo=sensitivityInfo,
                                  sensitivityRaw=sensitivityRaw,
                                  sensitivityProfiles=sensitivityProfiles,
                                  curationDrug=curationDrug,
                                  curationCell=curationCell,
                                  curationTissue=curationTissue,
                                  datasetType=c("both"),
                                  verify = TRUE)
  message("Done!")
  return(TGGATES_human)
}

#Example- creating toxicoSet for Human, DNA data:
TGGATES_humanDNA <- getTGGATEs(species = "Human", type = "DNA")
