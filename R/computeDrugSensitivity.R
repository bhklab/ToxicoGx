## This function computes intersected concentration range between a list of concentration ranges
.getCommonConcentrationRange <- function(doses)
{
  min.dose <- 0
  max.dose <- 10^100
  for(i in 1:length(doses))
  {
    min.dose <- max(min.dose, min(as.numeric(doses[[i]]), na.rm = TRUE), na.rm = TRUE)
    max.dose <- min(max.dose, max(as.numeric(doses[[i]]), na.rm = TRUE), na.rm = TRUE)
  }
  
  common.ranges <- list()
  for(i in 1:length(doses))
  {
    common.ranges[[i]] <- doses[[i]][
      which.min(abs(as.numeric(doses[[i]])-min.dose)):max(
        which(abs(as.numeric(doses[[i]]) - max.dose)==min(abs(as.numeric(doses[[i]]) - max.dose), na.rm=TRUE)))]
  }
  return(common.ranges)
}


updateMaxConc <- function(tSet){
  tSet@sensitivity$info$max.conc <- apply(tSet@sensitivity$raw[,,"Dose"], 1, max, na.rm=TRUE)
  return(tSet)
}
