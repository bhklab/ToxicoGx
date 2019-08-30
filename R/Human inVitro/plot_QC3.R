#module 144
drug <- "azathioprine"
conc <- 72.8

drug <- "diclofenac"
conc <- 400
# 
 drug <- "flutamide"
 conc <- 50
# 
 drug <- "isoniazid"
 conc <- 10000

#genes to plot
genes <-c("ENSG00000198431","ENSG00000138678", "ENSG00000170293","ENSG00000023909", "ENSG00000131389", "ENSG00000142657","ENSG00000187134","ENSG00000133597","ENSG00000144820", "ENSG00000077279")

#apply function to extract exprs matrix
values <- lapply(genes, function(gene){
#subset pehnodata for desired drugs
drug_subset <- subset(pData(tset@molecularProfiles$rna),drugid == drug,select=c(samplename, dose_level, individual_id,concentration))
#subset for ony high conc
drug_subset_high <- subset(drug_subset, concentration == conc)
#extracting exprs
assay <- exprs(tset@molecularProfiles$rna)
#subsetting exprs matrix
drug_subset$expression <- assay[gene,as.character(drug_subset$samplename)]
drug_subset_high$expression <- assay[gene,as.character(drug_subset_high$samplename)]
#ctrl rep
ctrlA <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="1")])
ctrlB <- na.omit(drug_subset$expression[select=c(drug_subset$dose_level == "Control" & drug_subset$individual_id=="2")])

highA <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="1")])
highB <- na.omit(drug_subset_high$expression[select=c(drug_subset_high$dose_level == "High" & drug_subset_high$individual_id=="2")])

ctrl <- rowMeans(cbind(ctrlA, ctrlB))
high <- rowMeans(cbind(highA, highB))


normalised_vehicle <- (high-ctrl)*100
return(normalised_vehicle)
})

values <- as.data.frame(do.call(rbind,values))

colnames(values) <- c(2,8,24)
rownames(values) <- genes


time <- c(2,8,24)
# rescale y axis
ylimit=range(0,100,200,300,400,500)
xlimit=range(0,4,8,12,16,20,24)
 
# graph

matplot(x = time, y = t(values)+100, col=c("violet", "green","purple","orange","turquoise", "green", "grey","pink","black","blue"), 
        pch=rep(21,ncol(values)), type=c("b"), lty=rep(1,ncol(values)),
        bg=c("violet", "green","purple","orange","turquoise", "green", "grey","pink","black","blue"),
        xlim=xlimit,ylim=ylimit,main="Isoniazid",
        xlab="Time",ylab="mRNA Level (% Vehicle)")
legend(23.2,519,legend=c("TXRND1","GPAT3","CMTM8","GCLM","SLC6A6", "PGD","AKR1C1","ADCK2","ADGRG7","DCX"), 
       col = c("violet", "green","purple","orange","turquoise", "green", "grey","pink","black","blue"),
       pch=rep(21,ncol(values)), pt.bg = c("violet", "green","purple","orange","turquoise", "green", "grey","pink","black","blue"),
       text.col = c("violet", "green","purple","orange","turquoise", "green", "grey","pink","black","blue"),
       lty=rep(1,ncol(values)),lwd = c(2,2,2,2,2,2,2,2,2), cex=0.45,xjust = 0.5)

