#Only includes samples where dosage of drug was measured in microMolar, OR the drug is gentamicin.
#ExpressionSet includes 2382 samples for Human, 3276 samples for Rat.

create_phenoData <- function(species=c("Human","Rat")){
  #load master phenoData file from TG-GATEs
  #Master phenoData file from TG-GATEs: #14 from https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html
  all_attribute <- read.delim("data/Open-tggates_AllAttribute.tsv",stringsAsFactors = F)
  
  #subset out unwanted samples (rows)
  all_attribute <- subset(all_attribute, all_attribute$BARCODE != "No ChipData" &
                            all_attribute$TEST_TYPE == "in vitro" &
                            (all_attribute$DOSE_UNIT == "ƒÊM" | all_attribute$COMPOUND_NAME == "gentamicin"))
  
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
  
  #Human
  if (species == "Human"){
    batch <- read_xlsx("data/nar-02356-data-e-2014-File006.xlsx", sheet = "Sup_table.2")
    batch <- batch[-1,]
    batch <- batch[,c(1,4)]
    names(batch) <- c("BARCODE", "CELL_NAME_TYPE_ID")
    batch <- subset(batch, batch$BARCODE != "No ChipData")
    all_attribute <- merge(all_attribute, batch, by.x = "BARCODE", by.y = "BARCODE")
    names(all_attribute)[names(all_attribute) == "CELL_NAME_TYPE_ID"] <- "batchid"
    
  } else if (species == "Rat"){ #Rat
    batch <- read_xlsx("data/nar-02356-data-e-2014-File006.xlsx", sheet = "Sup_table.3")
    batch <- batch[-c(1:3),]
    batch <- batch[,c(1,9)]
    colnames(batch) <- c("drugid", "batchid")
    colnames(batch) <- as.factor(colnames(batch))
    batch$drugid[batch$drugid=="TNFalpha"]<- "TNFƒ¿"
    
    all_attribute <- subset(all_attribute, all_attribute$ARR_DESIGN == "Rat230_2" & all_attribute$TEST_TYPE == "in vitro")
    all_attribute <- merge(all_attribute, batch, by.x = "COMPOUND_NAME", by.y = "drugid")
    
    conv <- readRDS("data/conversions_gentamicin.rds")
    #Include converted doses
    need_conversion <- subset(all_attribute,all_attribute$DOSE_UNIT != "ƒÊM",select = -c(DOSE,DOSE_UNIT))
    converted <- merge(need_conversion,conv,by.x="BARCODE",by.y="BARCODE")
    all_attribute <- subset(all_attribute,all_attribute$DOSE_UNIT == "ƒÊM")
    all_attribute <- rbind(all_attribute,converted)
  }
  all_attribute <- all_attribute[, names(all_attribute) != "DOSE_UNIT"]
  ##Merge drug mapping to lab annotations to phenoData
  all_attribute <- merge(all_attribute, drug_curation, by.x = "COMPOUND_NAME", by.y = "tggates.drugid")
  #CEL file name
  all_attribute$celfilename <- paste0(as.character(all_attribute$BARCODE),".CEL")
  # Label control or perturbation xptype
  all_attribute$xptype <- NA
  all_attribute$xptype[all_attribute$DOSE_LEVEL != "Control"] <- "perturbation"
  all_attribute$xptype[all_attribute$DOSE_LEVEL == "Control"] <- "control"
  #Substring barcode & assign cellids
  all_attribute$BARCODE<-substr(all_attribute$BARCODE,3,12)
  all_attribute$cellid <- "Hepatocyte"
  # all_attribute$cellid <- all_attribute$BARCODE
  #fix
  all_attribute$unique.drugid[all_attribute$unique.drugid=="TNFƒ¿"]<- "TNFa"
  #take out " hr" from duration column entries
  all_attribute$SACRI_PERIOD <- gsub(" hr","",all_attribute$SACRI_PERIOD)
  #UID
  all_attribute$UID<-paste("drugid_",all_attribute$COMPOUND.Abbr.,"_",all_attribute$unique.drugid,"_",all_attribute$cellid,"_",all_attribute$SACRI_PERIOD,"hr_rep",all_attribute$INDIVIDUAL_ID, sep="")
  
  #reorder columns
  all_attribute <- all_attribute[,c("BARCODE","ARR_DESIGN","EXP_ID","GROUP_ID","INDIVIDUAL_ID",
                                    "batchid","DOSE","DOSE_LEVEL","SACRI_PERIOD","cellid",
                                    "unique.drugid","COMPOUND_NAME","COMPOUND.Abbr.","COMPOUND_NO",
                                    "DNA...","LDH...","UID","SPECIES","TEST_TYPE","SEX_TYPE",
                                    "ORGAN_ID","MATERIAL_ID","celfilename","xptype","STRAIN_TYPE",
                                    "ADM_ROUTE_TYPE")]

  #Rename columns
  colnames(all_attribute) <- c("samplename","chiptype","exp_id","group_id","individual_id","batchid",
                               "concentration","dose_level","duration","cellid","drugid",
                               "tggates_drugid","drugid_abbr","drugid_no","DNA","LDH",
                               "UID","species","test_type","sex_type","organ_id","material_id",
                               "celfilename","xptype","STRAIN_TYPE","ADM_ROUTE_TYPE")
  
  #replace "ƒÊ" with "µ"
  # all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊM"] <- "µM"
  # all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊg/mL"] <- "µg/mL"
  
  #Rearrange order of rows by ascending order of barcodes
  all_attribute <- dplyr::arrange(all_attribute, samplename)
  #Assign rownames as samplename
  rownames(all_attribute)<-all_attribute$samplename
  
  return(all_attribute)
}

create_exprsData <- function(species=c("Human","Rat"), phenoData){
  if (species == "Human"){
    # install.packages("data/hgu133plus2hsensgcdf_22.0.0.tar.gz", repos = NULL, type = "source")
    library("hgu133plus2hsensgcdf")
    cdf <- "hgu133plus2hsensgcdf"
  } else if (species == "Rat"){
    # install.packages("data/rat2302rnensgcdf_23.0.0.tar.gz", repos = NULL, source = 'source')
    library("rat2302rnensgcdf")
    cdf <- "rat2302rnensgcdf"
  }
  celfn <- paste("CELfiles - ",species,"/", phenoData[,"celfilename"], sep="")
  
  ###########################################  NORMALIZATION  ############################################
  # eset <- just.rma(filenames = celfn, verbose = TRUE, cdfname = cdf)
  ########################################################################################################
  
  ########################################  ALREADY NORMALIZED  ##########################################
  eset <- readRDS(paste("data/eset_",species,"_",nrow(phenoData),".rds", sep = ""))
  ########################################################################################################
  
  storageMode(eset)<-"environment"
  #missingCEL is a data.frame containing the barcodes for all present samples
  # missingCEL <- phenoData[,c("celfilename"), drop=F]
  #subsetting samples
  # eset <- eset[,sampleNames(eset) %in% missingCEL$celfilename]
  #subsetting probes
  eset<-subset(eset, substr(rownames(eset@assayData$exprs), 0, 4) != "AFFX")
  #replace _at
  # rownames(eset)<-gsub("_at","",rownames(eset))
  #rename??
  colnames(eset@assayData$exprs)<-substr(colnames(eset@assayData$exprs),3,12)
  colnames(eset@assayData$se.exprs)<-substr(colnames(eset@assayData$se.exprs),3,12)
  rownames(eset@protocolData@data)<-substr(rownames(eset@protocolData@data),3,12)
  #lock eset@assayData environment again
  storageMode(eset)<-"lockedEnvironment"
  
  annotation(eset)<-"rna"
  
  return(eset)
}

create_featureData <- function(species=c("Human","Rat"), eset){
  if (species == "Human"){
    ensembl_data <- "hsapiens_gene_ensembl"
  } else if (species == "Rat"){
    ensembl_data <- "rnorvegicus_gene_ensembl"
  }
  CELgenes <- rownames(eset@assayData$exprs)
  
  ensembl<-useMart("ensembl", dataset = ensembl_data, host="uswest.ensembl.org",ensemblRedirect = FALSE)
  results <- getBM(attributes=c("external_gene_name","ensembl_gene_id","gene_biotype","entrezgene_id","external_transcript_name","ensembl_transcript_id"), filters = "ensembl_gene_id",values=gsub("_at","",CELgenes),mart=ensembl)
  uniqueBiomaRt<-results[!duplicated(results$ensembl_gene_id),]
  names(uniqueBiomaRt)<-c("gene_name", "gene_id", "gene_biotype", "EntrezGene.ID", "transcript_name", "transcript_id")
  
  if(species == "Rat"){
    names(uniqueBiomaRt)[1] <- "Symbol"
    uniqueBiomaRt$BEST <- NA
    uniqueBiomaRt<-arrange(uniqueBiomaRt,uniqueBiomaRt$gene_id)
    uniqueBiomaRt$gene_id <- paste(uniqueBiomaRt$gene_id,"_at", sep = "")
    rownames(uniqueBiomaRt)<-uniqueBiomaRt$gene_id
    
    return(uniqueBiomaRt)
  }
  
  labAnnot <- read.csv("data/annot_ensembl_all_genes.csv",header=TRUE,stringsAsFactors = F)[,-1, drop=F] #read in lab's gene annotation file
  CELnotbiomaRt<-unique(CELgenes)[!(unique(CELgenes) %in% uniqueBiomaRt$gene_id)] #in CELgenes but not in biomaRt output (66)
  newLabAnnot<-subset(labAnnot,labAnnot$gene_id %in% CELnotbiomaRt) #in CELgenes but not in biomaRt output but in lab annotation file (51)
  
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
  finalFeature$gene_id <- paste(finalFeature$gene_id,"_at", sep = "")
  rownames(finalFeature)<-finalFeature$gene_id
  
  return(finalFeature)
}

create_sensitivityProfiles <- function(phenoData){
  # create empty data frame of unique UID's
  sensProf <- data.frame(row.names = unique(phenoData$UID),
                         slope_recomputed = double(length(unique(phenoData$UID))),
                         auc_recomputed = double(length(unique(phenoData$UID))),
                         stringsAsFactors = F)
  sensProf[sensProf == 0] <- NA
  
  #subset out phenoData
  for (i in unique(phenoData$UID)){
    #samples of an experimental condition that are not control
    samples <- phenoData[phenoData$UID == i & phenoData$dose_level != "Control", , drop = F]
    if (nrow(samples) >= 3){
      sensProf[i,"slope_recomputed"] <- PharmacoGx::computeSlope(samples$concentration,samples$Viability)
      sensProf[i,"auc_recomputed"] <- PharmacoGx::computeAUC(samples$concentration, samples$Viability,
                                                             conc_as_log = FALSE, viability_as_pct = TRUE, area.type = "Actual")
    }
  }
  
  return(sensProf)
}

create_sensitivityRaw <- function(phenoData){
  conc_tested<-c("Control","Low","Middle","High")
  #Create temportary data.frames to use for processing
  #Dose Info from sensInfo & reformatting
  allDose <- subset(phenoData, select=c(UID, concentration,dose_level))
  #Viability Info (& dose info) from phenoData & reformatting
  allViability <- subset(phenoData, select=c(UID, concentration, Viability,dose_level))
  
  #Create empty dose array & format
  doseArray <- array(NA, dim = c(NROW(unique(phenoData$UID)),NROW(conc_tested)))
  rownames(doseArray)<-unique(phenoData$UID)
  colnames(doseArray)<-conc_tested #for now, colnames will be the dosages instead of "doses.1", "doses.2", etc.
  #Create empty viability array & format
  viabilityArray <- array(NA, dim=c(NROW(unique(phenoData$UID)),NROW(conc_tested)))
  rownames(viabilityArray)<-unique(phenoData$UID)
  colnames(viabilityArray)<-conc_tested
  
  for (i in 1:NROW(allDose)){ #for loop iterates as long as the number of rows (#UID's)
    #go down the row of doses by iteration
    dose_row_name <- allDose[i,"UID"] #dose_row_name is the UID at that iteration
    dosage <- allDose[i,"concentration"] #dosage is the dose (concentration) at that iteration
    dose_level <- allDose[i,"dose_level"] #dose_level is the dose level (control, low, middle, high) at that iteration
    doseArray[dose_row_name,dose_level]<-dosage #index into the correct spot using the above variables, and assign 'dosage' there
  }
  for (i in 1:NROW(allViability)){ #for loop iterates as long as the number of rows (#UID's)
    dose_row_name <- allViability[i,"UID"] #UID
    dosage <- allViability[i,"concentration"] #concentration
    viability<-allViability[i,"Viability"] #viability is on the same row, one column over
    dose_level <- allViability[i,"dose_level"]
    viabilityArray[dose_row_name,dose_level]<-viability #index into the correct spot using the above variables, and assign 'viability' there
  }
  
  colnames(doseArray)<-c("Control","doses1","doses2","doses3")
  colnames(viabilityArray)<-c("Control","doses1","doses2","doses3")
  
  #Combine doseArray and viabilityArray -> 3D array
  doseViabilityArray<-abind(doseArray, viabilityArray, along=3)
  dimnames(doseViabilityArray)[[3]]<-c("Dose","Viability")
  
  sensRaw <- doseViabilityArray
  
  return(sensRaw)
}

create_sensitivityInfo <- function(phenoData, doseArray){
  sensInfo <- subset(phenoData,select=c(UID, cellid, drugid, duration, individual_id))
  rownames(sensInfo) <- c()
  sensInfo <- unique(sensInfo)
  #set rownames to UIDs
  rownames(sensInfo) <- sensInfo$UID
  sensInfo <- subset(sensInfo,select=-c(UID))
  #rename columns
  names(sensInfo)<-c("cellid","drugid","duration_h","replicate")
  
  temp <- data.frame(row.names = unique(phenoData$UID),
                     Control = character(length(unique(phenoData$UID))),
                     Low = character(length(unique(phenoData$UID))),
                     Middle = character(length(unique(phenoData$UID))),
                     High = character(length(unique(phenoData$UID))),
                     stringsAsFactors = F)
  temp[temp == ""] <- NA
  
  for (i in unique(phenoData$UID)){
    c <- phenoData[phenoData$UID == i, "dose_level"]
    if ("Control" %in% c){ temp[i, "Control"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Control"] }
    if ("Low" %in% c) { temp[i, "Low"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Low"] }
    if ("Middle" %in% c) { temp[i, "Middle"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "Middle"] }
    if ("High" %in% c) { temp[i, "High"] <- phenoData$samplename[phenoData$UID == i & phenoData$dose_level == "High"] }
  }
  
  sensInfo <- merge(sensInfo, temp, by.x = 0, by.y = 0, sort = F)
  rownames(sensInfo)<-sensInfo$Row.names
  sensInfo<-subset(sensInfo,select=-c(Row.names))
  
  return (sensInfo)
}

create_curationDrug <- function(phenoData){
  curationDrug <- unique(subset(phenoData, select=c(drugid, tggates_drugid)))
  rownames(curationDrug) <- curationDrug$drugid
  names(curationDrug) <- c("unique.drugid", "tggates.drugid")
  
  return(curationDrug)
}

create_curationCell <- function(phenoData){
  curationCell <- unique(subset(phenoData, select=c(cellid)))
  curationCell$tggates.cellid <- curationCell$cellid
  names(curationCell) <- c("unique.cellid", "tggates.cellid")
  rownames(curationCell) <- curationCell$unique.cellid
  
  return(curationCell)
}

create_curationTissue <- function(phenoData){
  curationTissue <- unique(subset(phenoData, select=c(organ_id)))
  curationTissue$tggates.tissueid <- "Liver"
  names(curationTissue)[1] <- "unique.tissueid"
  rownames(curationTissue) <- "Hepatocyte"
  
  return(curationTissue)
}

create_drug <- function(phenoData){
  drug <- unique(subset(phenoData, select=c(drugid, drugid_abbr, drugid_no)))
  rownames(drug) <- drug$drugid
  
  return(drug)
}

create_cell <- function(phenoData){
  cell <- unique(subset(phenoData,select=c(cellid, organ_id, material_id, species, test_type)))
  names(cell)<-c("cellid","tissueid","materialid", "species","testType")
  cell$tissueid<-"Liver"
  rownames(cell) <- cell$cellid
  
  return(cell)
}

getTGGATEs <- function(species=c("Human","Rat"),
                       type=c("DNA", "LDH")){
  library(dplyr)
  library(PharmacoGx)
  library(gdata)
  library(readxl)
  library("xml2")
  library(biomaRt)
  library(affy)
  library(abind)
  
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
  pData(eset) <- phenoData
  fData(eset) <- featureData
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
  TGGATES <- ToxicoSet(paste("TGGATES ",species," ",type, sep = ""),
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
  return(TGGATES)
}

# EXAMPLE -
tggates <- getTGGATEs(species = "Human", type = "DNA")
