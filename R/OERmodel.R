# OER function for a given hypoxic concentration by the user

###TODO:: Figure out why package is called inside itself?
#library(RadioGx)

OERmodel <- function(){

  pO2 <- readline(prompt="Enter the oxygen concentration (between 0.1 and 10): ")
  pO2 <- as.numeric(pO2)
  if (is.na(pO2)) {
    print("Error")}

  if (pO2<0.1 | pO2>10){
    print("Please enter a value between 0.1 and 10")
  }

  OER_m = 3
  K_m = 3
  a = ((OER_m*pO2)+K_m)/(pO2+K_m)
  OMF = (1/OER_m)*a

  D <- as.numeric(c("0","1","2","3","4","5","6","8","10"))
  SF1 = exp(-0.3*D*OMF-(0.03*D*D*OMF))
  #RadioGx::computeAUC(D,SF1)
  pdf("HyxpoxiaPlot.pdf")
  RadioGx::doseResponseCurve(Ds=list("Hypoxia" = D),
                              SFs=list("Hypoxia" = SF1), plot.type="Actual",
                              legends.label = NULL,title = "Effect of Hypoxia",
                              cex = 1.55,cex.main = 1.75,lwd = 2)
  dev.off()

}

###TODO::Make sure this function works
#OERmodel()
