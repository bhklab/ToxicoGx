.makeGrittyGuess <- function(lower_bounds, upper_bounds, D, SF) {
  dFromD2 <- abs(D - 2)
  dFromSF10 <- abs(exp(SF) - 0.1)
  D10Ind <- which.min(dFromSF10)[1]
  SF2Ind <- which.min(dFromD2)[1]

  if (D10Ind == SF2Ind) {
    dFromD2[SF2Ind] <- dFromD2[SF2Ind] + dFromD2[1] #make it not the smallest anymore
    SF2Ind <- which.min(dFromD2)[1]
  }

  DSF2 <- D[SF2Ind]
  SF2 <- SF[SF2Ind]
  D10 <- D[D10Ind]
  SFD10 <- SF[D10Ind]

  if(any(c(DSF2, D10, D[D10Ind] - D[SF2Ind])==0)){
    return(c(0.25, 0.25))
  }

  gritty_guess <- pmin(pmax(c((SF[D10Ind] * D[SF2Ind] ^ 2 - SF[SF2Ind] * D[D10Ind] ^ 2) / D[SF2Ind] / D[D10Ind] / (D[D10Ind] - D[SF2Ind]),
                              (SF[D10Ind] * D[SF2Ind] - SF[SF2Ind] * D[D10Ind]) / D[SF2Ind] / D[D10Ind] / (D[SF2Ind] - D[D10Ind])), lower_bounds), upper_bounds) # assumes the SF2Indth point is SF2 and D10Indth point is D10 and imputes alpha, beta from that assumption unless either would thus be out of bounds


  return(gritty_guess)
}
