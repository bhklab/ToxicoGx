#' Plos One QC1

# Read tset
#library(PharmacoGx)
tset <- readRDS("rds/TGGATES_ratDNA.rds")

# gene id = ENSRNOG00000008001	
# drug = phenytoin	

# expression dataframe
table <- tset@molecularProfiles$rna@assayData$exprs

#subset dataframe for probe
singleGeneRow <- table[rownames(table)=="ENSRNOG00000012772", , drop = FALSE]

ccl4 <- rownames(p[p$drugid_abbr=="CCL4",])
aspirin <- rownames(p[p$drugid=="aspirin",])
valproicacid <- rownames(p[p$drugid=="Valproic.acid",])
azathioprine <- rownames(p[p$drugid=="azathioprine",])
cicloheximide <- rownames(p[p$drugid=="cicloheximide",])
etoposide <- rownames(p[p$drugid=="etoposide",])

drug <- etoposide #CHANGE THIS#########

# subset dataframe columns specific to drug
drugTimePoints <- subset(singleGeneRow, select = c(drug), drop = F)

# group 3 pts into a line by indexing
# subset control
control <- drugTimePoints[,c(1,9,17), drop = FALSE]
# subset low 
low <- drugTimePoints[,c(3,11,19), drop = FALSE]
# subset mid
mid <- drugTimePoints[,c(5,13,21), drop = FALSE]
# subset high
high <- drugTimePoints[,c(7,15,23), drop = FALSE]
# subset control
control1 <- drugTimePoints[,c(2,10,18), drop = FALSE]
# subset low 
low1 <- drugTimePoints[,c(4,12,20), drop = FALSE]
# subset mid
mid1 <- drugTimePoints[,c(6,14,22), drop = FALSE]
# subset high
high1 <- drugTimePoints[,c(8,16,24), drop = FALSE]

time <- c(2,8,24)
# rescale y axis
down <- 10
up <- 12
ylimit = range(down,up)
#plot
plot(time, control,lty = 2, pch = 16, col = "red", ylim = ylimit, xlim= c(2,26.5),type='n', axes=T, main = "Etoposide on Rat Gene NQO1 (ENSRNOG00000012772)", cex.main = 0.8, xlab = "Time", ylab = "Expression", ann = T, xaxt = 'n', yaxt = 'n')

# axes
# x axis
axis(side = 1, at = c(2,8,24), labels = c(2,8,24), tick = T, lty = "solid",lwd.ticks = 1, gap.axis = 0)
# y axis
axis(side = 2, at = c(down,up), labels = c(down,up), tick = T, lty = 1, col.ticks = T)

# replicate A
points(time, control,col="red", pch = 16)
lines(time, control, lty = 2, col ="red", lwd = 2)
points(time, low,col="green3", pch = 3)
lines(time, low, lty = 2, col ="green3", lwd = 2)
points(time, mid,col="blue", pch = 23, bg = "blue")
lines(time, mid, lty = 2, col ="blue", lwd = 2)
points(time, high,col="cyan", pch = 22, bg = "cyan")
lines(time, high, lty = 2, col ="cyan", lwd = 2)
# replicate B
points(time, control1, type="o",col="red", pch = 17)
lines(time, control1, type="o", lty=1, col="red", lwd = 2)
points(time, low1,col="green3", pch = 4)
lines(time, low1, lty = 1, col ="green3", lwd = 2)
points(time, mid1,col="blue", pch = 25, bg = "blue")
lines(time, mid1, lty = 1, col ="blue", lwd = 2)
points(time, high1,col="cyan", pch = 21, bg = "cyan")
lines(time, high1, lty = 1, col ="cyan", lwd = 2)

# legend
legend(23.5,10.8,legend=c("0 Ctrl CtrlA","0 Ctrl CtrlB","1 Low LowA","1 Low LowB","2 Mid MidA","2 Mid MidB","3 High HighA","3 High HighB"), col = c("red","red","green3","green3","blue","blue","cyan","cyan"),pch = c(16,17,3,4,23,25,22,21), pt.bg = c("red","red","green3","green3","blue","blue","cyan","cyan"),text.col = "black",lty = c(2,1,2,1,2,1,2,1,2),lwd = c(2,2,2,2,2,2,2,2), cex=0.28)

