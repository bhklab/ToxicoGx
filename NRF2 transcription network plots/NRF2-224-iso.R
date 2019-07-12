#module 224
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
genes <-c("ENSG00000156535","ENSG00000087842", "ENSG00000065600","ENSG00000065833", "ENSG00000133393", "ENSG00000213190","ENSG00000116701")

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
ylimit=range(0,100,200,300,330)
xlimit=range(0,4,8,12,16,20,24)

# graph

legendnames <- c("CD109", 'PIR','TMEM206','ME1','FOPNL','MLLT11','NCF2','LACC1')
colours <- c("purple", "red","green","violet","yellow", "turquoise", "lime green","grey")
matplot(x = time, y = t(values)+100, col=colours, 
        pch=rep(21,ncol(values)), type=c("b"), lty=rep(1,ncol(values)), lwd=rep(3,ncol(values)),
        bg=colours,
        xlim=xlimit,ylim=ylimit,main="Isoniazid (Mod 224)",
        xlab="Time",ylab="mRNA Level (% Vehicle)")
legend(23.1,343,legend=legendnames, 
       col = colours,
       pch=rep(21,ncol(values)), pt.bg = colours,
       text.col = 'black',
       lty=rep(1,ncol(values)),lwd = c(3,3,3,3,3,3,3,3), cex=0.54,xjust = 0.5)

