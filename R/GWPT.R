#' @title Genome-Wide Pleiotropy Test
#' @description
#' This function performs a genome-wide pleiotropy test (GWPT) after Mendelian randomization. It offers an option for a two-mixture model, where the residual is chosen as the smaller one resulting from the two causal effect estimates from two mixtures.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param theta The causal effect estimate.
#' @param theta.cov The covariance matrix of the causal effect estimate.
#'
#' @return A list with two components:
#' \item{BETA}{The estimated residual values.}
#' \item{SE}{The standard errors of the residual estimates.}
#'
#' @export

GWPT <- function(by, byse, bX, bXse, Rxy, theta, theta.cov) {

bZ <- cbind(bX,by)
bZse <- cbind(bXse,byse)
vartheta <- c(-theta,1)
residual <- bZ %*% vartheta

var_residual <- residual * 0
Covtheta <- diag(length(vartheta) - 1) * 0

if (any(theta != 0)) {
  Covtheta[which(theta != 0), which(theta != 0)] <- theta.cov[which(theta != 0), which(theta != 0)]
}
Covtheta <- Matrix::bdiag(Covtheta,0)

pb <- txtProgressBar(min = 0, max = length(var_residual), style = 3)

for (i in 1:nrow(bZ)) {
  setTxtProgressBar(pb, i)
  beta <- as.vector(bZ[i, ])
  se <- as.vector(bZse[i, ])
  G <- t(Rxy * se) * se
  var_residual[i] <- sum(vartheta * (G %*% vartheta)) + sum(beta * (Covtheta %*% beta))
}

close(pb)

# Result list
A <- list()
A$BETA <- residual
A$SE <- sqrt(var_residual)

return(A)
}
