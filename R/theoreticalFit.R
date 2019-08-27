# library("caTools")
# library("xlsx")
# 
# aucs <- read.xlsx2('XRT_CTD2_Dose_Response.xlsx', 1, stringsAsFactors = FALSE)
# 
# doses <- c(1, 2, 3, 4, 5, 6, 8, 10)
# cols <- matrix("uninitialized", nrow = 1, ncol = length(doses))
# 
# for (i in seq_along(doses)) {
#   cols[i] <- paste0("X", as.character(doses[i]))
#   if (doses[i] == 1) {
#     cols[i] <- paste0(cols[i], ".1.")
#   }
# }
# 
# aucLQs <- matrix(NA, nrow = 1, ncol = nrow(aucs))
# alphas <- matrix(NA, nrow = 1, ncol = nrow(aucs))
# betas <- matrix(NA, nrow = 1, ncol = nrow(aucs))
# a <- matrix(NA, nrow = 1, ncol = length(doses))
# 
# for (i in seq_along(aucs)) {
#   a <- apply(aucs[i,cols], 2, as.numeric)
#   nlsFit <- nls(a ~ exp(-alphas * doses - betas * doses ^ 2), #alpha, beta E (0, 0.3]; beta < alpha
#                 start = list(alphas = 0.01, betas = 0.001))
#   alphas[i] = coef(nlsFit)[1]
#   betas[i] = coef(nlsFit)[2]
#   aucLQs[i] <- trapz(doses, exp(-alphas[i] * doses - betas[i] * doses ^ 2))
# }
# 
# names(aucLQs) <- aucs$cell_line
# write.csv(aucLQs, file = "aucLQs.csv")