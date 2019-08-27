.linearQuadraticInv <- function(SF, pars, SF_as_log = TRUE) {
  
  if (!SF_as_log) {
    if (SF < 0) {
      stop("Survival fraction must be nonnegative.")
    } else {
      SF <- log(SF)
    }
  }

  if (SF > 0) {
    stop("Positive log survival fraction ", SF,  "cannot be reached at any dose of radiation with linear quadratic paramaters alpha, beta > 0.")
  } else {
    if (pars[[2]] == 0) {
      if (pars[[1]] == 0) {
        if (SF == 1) {
          return(0)
        } else {
          stop(paste0("Survival fraction ", SF, " cannot be reached at any dose of radiation with linear-quadratic parameters alpha = ", pars[[1]], " and beta = ", pars[[2]], "."))
        }
      } else {
        return(-SF / pars[[1]])
      }
    } else {
      return((sqrt(pars[[1]] ^ 2 - 4 * pars[[2]] * SF) - pars[[1]]) / 2 / pars[[2]])
    }
  }
}
