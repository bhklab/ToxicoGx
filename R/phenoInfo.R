library(dplyr)

#Load Open-tggates_AllAttribute.tsv from LSDB database
all_attribute <- read.delim("data/Open-tggates_AllAttribute.tsv",stringsAsFactors = FALSE)

#Subset out rows with "No ChipData", rat samples, or dosage units are not "ƒÊg/kg" -> should result in 2573 entries
all_attribute <- subset(all_attribute,all_attribute$BARCODE != "No ChipData" & all_attribute$SPECIES == "Human" & all_attribute$DOSE_UNIT != "ƒÊg/kg")

#Subset out unnecessary columns
all_attribute <- subset(all_attribute,select=-c(SIN_REP_TYPE,ANIMAL_AGE.week.,TERMINAL_BW.g.,LIVER.g.,KIDNEY_TOTAL.g.,
                                                KIDNEY_R.g.,KIDNEY_L.g.,RBC.x10_4.ul.,Hb.g.dL.,Ht...,MCV.fL.,MCH.pg.,
                                                MCHC...,Ret...,Plat.x10_4.uL.,WBC.x10_2.uL.,Neu...,Eos...,Bas...,
                                                Mono...,Lym...,PT.s.,APTT.s.,Fbg.mg.dL.,ALP.IU.L.,TC.mg.dL.,TG.mg.dL.,
                                                PL.mg.dL.,TBIL.mg.dL.,DBIL.mg.dL.,GLC.mg.dL.,BUN.mg.dL.,CRE.mg.dL.,
                                                Na.meq.L.,K.meq.L.,Cl.meq.L.,Ca.mg.dL.,IP.mg.dL.,TP.g.dL.,RALB.g.dL.,
                                                A.G,AST.IU.L.,ALT.IU.L.,LDH.IU.L.,GTP.IU.L.))

#Substring barcode
all_attribute$BARCODE<-substr(all_attribute$BARCODE,3,12)

#Add necessary columns:
#Batch ID
cellLots<-read.csv("data/Cell Lots 1.csv",stringsAsFactors = FALSE)
all_attribute <- merge(all_attribute, cellLots, by.x = "BARCODE", by.y = "BARCODE")
#Converted concentration
#Load conversions
conv <- readRDS("data/conversions.rds")
#subset out rows that need compound dosages converted:
need_conversion <- subset(all_attribute,all_attribute$COMPOUND_NAME %in% c("LPS","TNFƒ¿","aflatoxin B1","hepatocyte growth factor, human","interferon alpha, human","interleukin 1 beta, human","interleukin 6, human","phalloidin","transforming growth factor beta 1","tunicamycin"))
need_conversion <- subset(need_conversion,select=-c(DOSE,DOSE_UNIT))
converted <- merge(need_conversion,conv,by.x="BARCODE",by.y="X")

all_attribute <- subset(all_attribute,!(all_attribute$COMPOUND_NAME %in% c("LPS","TNFƒ¿","aflatoxin B1","hepatocyte growth factor, human","interferon alpha, human","interleukin 1 beta, human","interleukin 6, human","phalloidin","transforming growth factor beta 1","tunicamycin")))
all_attribute$RECOMP_DOSE <- all_attribute$DOSE

all_attribute <- rbind(all_attribute,converted)
#cellid (barcode)
all_attribute$cellid <- all_attribute$BARCODE
#UID
all_attribute$UID<-paste0("drugid_",all_attribute$COMPOUND.Abbr.,"_",all_attribute$COMPOUND_NAME,"_",all_attribute$BARCODE,"_",all_attribute$DOSE,"_",gsub(" hr","",all_attribute$SACRI_PERIOD),"hr_rep",all_attribute$INDIVIDUAL_ID)
#CEL file name
all_attribute$celfilename <- paste0("00",as.character(all_attribute$BARCODE),".CEL")
#xptype
all_attribute$xptype <- "both"

#Rearrange columns
all_attribute <- all_attribute[,c(1,2,3,4,5,22,17,18,23,19,16,24,8,9,10,20,21,25,11,12,13,6,7,26,27,14,15)]

#Rename columns
colnames(all_attribute) <- c("samplename","chiptype","exp_id","group_id","individual_id","batchid",
                             "concentration_old","concentration_old_units","concentration","dose_level",
                             "duration","cellid","drugid","drugid_abbr","drugid_no","DNA...","LDH...",
                             "UID","species","test_type","sex_type","organ_id","material_id",
                             "celfilename","xptype","STRAIN_TYPE","ADM_ROUTE_TYPE")

#Rearrange order of rows by ascending order of barcodes
all_attribute <- arrange(all_attribute, samplename)

#Change rownames to barcodes
rownames(all_attribute)<-all_attribute$samplename

#replace "ƒÊ" with "µ"
all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊM"] <- "µM"
all_attribute$concentration_old_units[all_attribute$concentration_old_units == "ƒÊg/mL"] <- "µg/mL"

saveRDS(all_attribute,"data/phenoData.rds")
