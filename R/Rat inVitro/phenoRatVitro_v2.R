ls()
remove(list = ls())

library(dplyr)

#Load Open-tggates_AllAttribute.tsv from LSDB database
all_attribute <- read.delim("data/Open-tggates_AllAttribute.tsv",stringsAsFactors = FALSE)

View(all_attribute)

#View(unique(all_attribute$DOSE_UNIT))
#Subset out rows with "No ChipData", rat samples, (kept dosage units are not "f g/kg" -> should result in 2573 entries // & all_attribute$DOSE_UNIT !%in% "mg/kg", "f g/kg","%" conversion will be replaced by na
rat_attribute <- subset(all_attribute,all_attribute$BARCODE != "No ChipData" & all_attribute$SPECIES == "Rat" & !(all_attribute$DOSE_UNIT %in% c( "mg/kg", "f g/kg","%" ))) # took out : & all_attribute$DOSE_UNIT != "f g/kg"

#Subset out unnecessary columns
rat_attribute <- subset(rat_attribute,select=-c(SIN_REP_TYPE,ANIMAL_AGE.week.,TERMINAL_BW.g.,LIVER.g.,KIDNEY_TOTAL.g.,
                                                KIDNEY_R.g.,KIDNEY_L.g.,RBC.x10_4.ul.,Hb.g.dL.,Ht...,MCV.fL.,MCH.pg.,
                                                MCHC...,Ret...,Plat.x10_4.uL.,WBC.x10_2.uL.,Neu...,Eos...,Bas...,
                                                Mono...,Lym...,PT.s.,APTT.s.,Fbg.mg.dL.,ALP.IU.L.,TC.mg.dL.,TG.mg.dL.,
                                                PL.mg.dL.,TBIL.mg.dL.,DBIL.mg.dL.,GLC.mg.dL.,BUN.mg.dL.,CRE.mg.dL.,
                                                Na.meq.L.,K.meq.L.,Cl.meq.L.,Ca.mg.dL.,IP.mg.dL.,TP.g.dL.,RALB.g.dL.,
                                                A.G,AST.IU.L.,ALT.IU.L.,LDH.IU.L.,GTP.IU.L.))

#Substring barcode
rat_attribute$BARCODE<-substr(rat_attribute$BARCODE,3,12)
rat_attribute$COMPOUND_NAME[rat_attribute$COMPOUND_NAME=="TNFfø"]<- "TNFa"
#Add necessary columns:
#Batch ID
#cellLots<-read.csv("data/Cell Lots 1.csv",stringsAsFactors = FALSE) # for human 
#View(cellLots)
#all_attribute <- merge(all_attribute, cellLots, by.x = "BARCODE", by.y = "BARCODE")# for human

View(rat_attribute)

# Column Batch ID #
# Read mapping values
batchId <- read.csv("https://raw.githubusercontent.com/cjtangyx/ToxicoGx/master/Rat%20inVitro/data/phenoRatVitro.csv",stringsAsFactors = FALSE)
#Rename columns for convenient mapping
colnames(batchId) <-as.character(c('drugid', 'drugid_abbr', 'drugid_no', 'medium_type' ,'cell','low_dose','med_dose','high_dose','batch_id'))
colnames(batchId) <- as.factor(colnames(batchId))
# Eliminate unnecessary rows and columns
batchId <- batchId[-c(1:3),] 
batchId <- batchId[,-c(10:10)]
batchId$drugid[batchId$drugid=="TNFalpha"]<- "TNFa"

# Merge batchId with ratAttributes to create new column in dataframe
rat_attribute <- merge(rat_attribute, batchId, by.x = "COMPOUND_NAME", by.y = "drugid")
View(rat_attribute)
length(unique(rat_attribute$COMPOUND_NAME))
colnames(rat_attribute) <- as.factor(colnames(rat_attribute))


#Converted concentration
#Load conversions  

conversion <- read.csv("data/conversions_rat_v1.csv")
Conv_adj <- subset(conversion, !(conversion$X %in% c("3017922030",	"3017923001",	"3017923004",	"3017923005")))
View(Conv_adj)

#subset out rows that need compound dosages converted:
need_conversion <- subset(rat_attribute,rat_attribute$COMPOUND_NAME %in% c("LPS","TNFa","gentamicin","tunicamycin","phalloidin"))
need_conversion <- subset(need_conversion,select=-c(DOSE,DOSE_UNIT))
converted <- merge(need_conversion,Conv_adj,by.x="BARCODE",by.y="X")
converted<- subset(converted[,-c(28)])
View(converted)

rat_attribute_noneedconv <- subset(rat_attribute,!(rat_attribute$COMPOUND_NAME %in% c("LPS","TNFa","gentamicin","tunicamycin","phalloidin")))
rat_attribute_noneedconv$RECOMP_DOSE <- rat_attribute_noneedconv$DOSE
View(rat_attribute_noneedconv)

all_rat_attribute <- rbind(rat_attribute_noneedconv,converted)
View(all_rat_attribute)
write.csv(all_rat_attribute,"all_rat_attribute_v1.csv")
View(unique(all_rat_attribute$COMPOUND_NAME))
#cellid (barcode)
all_rat_attribute$cellid <- all_rat_attribute$BARCODE
#UID
all_rat_attribute$UID<-paste0("drugid_",all_rat_attribute$COMPOUND.Abbr.,"_",all_rat_attribute$COMPOUND_NAME,"_",all_rat_attribute$BARCODE,"_",all_rat_attribute$DOSE,"_",gsub(" hr","",all_rat_attribute$SACRI_PERIOD),"hr_rep",all_rat_attribute$INDIVIDUAL_ID)
#CEL file name
all_rat_attribute$celfilename <- paste0("00",as.character(all_rat_attribute$BARCODE),".CEL")
#xptype
all_rat_attribute$xptype <- "both"
View(all_rat_attribute)
write.csv(all_rat_attribute,"all_rat_attribute_v2.csv")
#Rearrange columns
all_rat_attribute <- all_rat_attribute[,c(2,3,4,5,6,29,	17,18,30,19,16,31,1,9,10,	20,	21,32,11,12,13,7,8,33,34,14,15)]


#Rename columns
colnames(all_rat_attribute) <- c("samplename","chiptype","exp_id","group_id","individual_id","batchid",
                             "concentration_old","concentration_old_units","concentration","dose_level",
                             "duration","cellid","drugid","drugid_abbr","drugid_no","DNA","LDH",
                             "UID","species","test_type","sex_type","organ_id","material_id",
                             "celfilename","xptype","STRAIN_TYPE","ADM_ROUTE_TYPE")

#Rearrange order of rows by ascending order of barcodes
#all_ratattribute <- arrange(all_rat_attribute, samplename)
#View(all_ratattribute)
#Change rownames to barcodes
#rownames(all_ratattribute)<-all_ratattribute$samplename

View(all_rat_attribute)

curationdrug<- read.csv("curationDrug.csv")
all_rat_attribute <- merge(all_rat_attribute, curationdrug, by.x = "drugid", by.y = "tggates.drugid")
View(all_rat_attribute)

#drop drugid column and keep unique.drugid column
ratAttributes_r1<-all_rat_attribute[,c(2:29)] #drop col compound name
all_ratattribute<-ratAttributes_r1[,-c(27)] #drop col compound name
colnames(all_ratattribute)[27]<-"drugid" # change the last col name to drugid from unique.drugid
View(all_ratattribute)

#replace "f " with "µ"
all_ratattribute$concentration_old_units[all_ratattribute$concentration_old_units == "f M"] <- "µM"
all_ratattribute$concentration_old_units[all_ratattribute$concentration_old_units == "fEg/mL"] <- "µg/mL"

View(all_ratattribute)

#Rearrange order of rows by ascending order of barcodes
all_ratattribute <- arrange(all_ratattribute, samplename)
View(all_ratattribute)
#Change rownames to barcodes
rownames(all_ratattribute)<-all_ratattribute$samplename

write.csv(all_ratattribute,"all_ratattribute_v4.csv")
saveRDS(all_ratattribute,"rds/phenoData_rat.rds")
pheno<-readRDS("rds/phenoData_rat.rds")
View(pheno)
