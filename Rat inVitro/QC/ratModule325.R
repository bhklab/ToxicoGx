#module 224
drug <- "azathioprine"
conc <- 3.60

drug <- "diclofenac"
conc <- 400
# 
drug <- "flutamide"
conc <- 75
# 
drug <- "isoniazid"
conc <- 10000

#genes to plot
genes <-c("ENSRNOG00000018054","ENSRNOG00000022189", "ENSRNOG00000012772","ENSRNOG00000009934", "ENSRNOG00000001701")

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

legendnames <- c('CBR3','EMC3','NQO1','HTAT1P2','F2RL2')
colours <- c("purple", "red","violet","yellow", "turquoise")
matplot(x = time, y = t(values)+100, col=colours, 
        pch=rep(21,ncol(values)), type=c("b"), lty=rep(1,ncol(values)), lwd=rep(3,ncol(values)),
        bg=colours,
        xlim=xlimit,ylim=ylimit,main="Rat Isoniazid (Mod 325)",
        xlab="Time",ylab="mRNA Level (% Vehicle)")
legend(23.1,520,legend=legendnames, 
       col = colours,
       pch=rep(21,ncol(values)), pt.bg = colours,
       text.col = 'black',
       lty=rep(1,ncol(values)),lwd = c(3,3,3,3,3,3), cex=0.54,xjust = 0.78)

