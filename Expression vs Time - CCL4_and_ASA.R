#Load tSet
tset <- readRDS("data/TGGATES.rds")

##########################################################################################################################################

#Plot Expression vs Time -> Carbon Tetrachloride (CCL4) for ENSG00000140465
#Subset only for samples tested with carbon tetrachloride from pData
CCL4 <- subset(pData(tset@molecularProfiles$rna),
               drugid == "carbon tetrachloride",
               select=c(samplename, dose_level, individual_id))
#Assign expression values of tSet to assay
assay <- exprs(tset@molecularProfiles$rna)
#Select for expression values of tSet which correspond to samples tested with carbon tetrachloride
CCL4$expression <- assay["ENSG00000140465",as.character(CCL4$samplename)]
#Select for expression values that correspond to the correct treatment conditions
ctrlA <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Control" & CCL4$individual_id=="1")])
ctrlB <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Control" & CCL4$individual_id=="2")])
lowA <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Low" & CCL4$individual_id=="1")])
lowB <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Low" & CCL4$individual_id=="2")])
medA <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Middle" & CCL4$individual_id=="1")])
medB <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "Middle" & CCL4$individual_id=="2")])
highA <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "High" & CCL4$individual_id=="1")])
highB <- na.omit(CCL4$expression[select=c(CCL4$dose_level == "High" & CCL4$individual_id=="2")])
#Plotting
#Set x-limit, y-limit, and time points
xlimit=range(0:30)
ylimit=range(5:11)
time=c(2,8,24)
#Plot using matplot
matplot(x = time, 
        y = matrix(c(ctrlA,ctrlB, lowA, lowB, medA, medB, highA, highB), ncol=8), 
        col=c("red","red","green","green","blue","blue","cyan","cyan"), 
        pch=c(21,24,3,4,23,25,22,8), 
        type=c("b"), 
        lty=c(1,2,3,4,5,4,1,2),
        bg=c("red","red","green","green","blue","blue","cyan","cyan"),
        xlim=xlimit,
        ylim=ylimit,
        main="cytochrome P450 (ENSG00000140465) - carbon tetrachloride (CCL4)",
        xlab="Time",
        ylab="Expression",
        xaxt='n')
#Label x-axis
axis(1,at=c(2,8,24),labels=c(2,8,24))
#Set Legend
legend(26,10,legend=c("ctrlA","ctrlB","lowA","lowB","medA","medB","highA","highB"),
       col=c("red","red",'green','green','blue','blue','cyan','cyan'),
       pch=c(21,24,3,4,23,25,22,8),
       pt.bg=c("red","red",'green','green','blue','blue','cyan','cyan'),
       lty=c(1,2,3,4,5,4,1,2),
       cex=0.75)

##########################################################################################################################################

#Plot Expression vs Time -> Aspirin (ASA) for ENSG00000169213
#Subset only for samples tested with aspirin from pData
aspirin <- subset(pData(tset@molecularProfiles$rna),
                  drugid == "aspirin",
                  select=c(samplename, dose_level, individual_id))
#Assign expression values of tSet to assay
assay <- exprs(tset@molecularProfiles$rna)
#Select for expression values of tSet which correspond to samples tested with aspirin
aspirin$expression <- assay["ENSG00000169213",as.character(aspirin$samplename)]
#Select for expression values that correspond to the correct treatment conditions
ctrlA <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Control" & aspirin$individual_id=="1")])
ctrlB <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Control" & aspirin$individual_id=="2")])
lowA <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Low" & aspirin$individual_id=="1")])
lowB <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Low" & aspirin$individual_id=="2")])
medA <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Middle" & aspirin$individual_id=="1")])
medB <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "Middle" & aspirin$individual_id=="2")])
highA <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "High" & aspirin$individual_id=="1")])
highB <- na.omit(aspirin$expression[select=c(aspirin$dose_level == "High" & aspirin$individual_id=="2")])
#Plotting
#Set x-limit, y-limit, and time points
xlimit=range(0:30)
ylimit=range(8.4:9.4)
time=c(2,8,24)
#Plot using matplot
matplot(x = time, 
        y = matrix(c(ctrlA,ctrlB, lowA, lowB, medA, medB, highA, highB), ncol=8), 
        col=c("red","red","green","green","blue","blue","cyan","cyan"), 
        pch=c(21,24,3,4,23,25,22,8), 
        type=c("b"), 
        lty=c(1,2,3,4,5,4,1,2),
        bg=c("red","red","green","green","blue","blue","cyan","cyan"),
        xlim=xlimit,
        ylim=ylimit,
        main="RAB3B, member RAS oncogene family (ENSG00000169213) - aspirin (ASA)",
        xlab="Time",
        ylab="Expression",
        xaxt='n')
#Label x-axis
axis(1,at=c(2,8,24),labels=c(2,8,24))
#Set legend
legend(26,9.2,legend=c("ctrlA","ctrlB","lowA","lowB","medA","medB","highA","highB"),
       col=c("red","red",'green','green','blue','blue','cyan','cyan'),
       pch=c(21,24,3,4,23,25,22,8),
       pt.bg=c("red","red",'green','green','blue','blue','cyan','cyan'),
       lty=c(1,2,3,4,5,4,1,2),
       cex=0.75)