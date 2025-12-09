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
#' @import CppMatrix
#' @return A list with two components:
#' \item{BETA}{The estimated residual values.}
#' \item{SE}{The standard errors of the residual estimates.}
#'
#' @export

GWPT <- function(by, byse, bX, bXse, Rxy, theta, theta.cov) {

bZ <- cbind(bX,by)
bZse <- cbind(bXse,byse)
vartheta <- c(-theta,1)
residual <- matrixVectorMultiply(bZ,vartheta)

var_residual <- residual * 0
Covtheta <- diag(length(vartheta) - 1) * 0

if (any(theta != 0)) {
  Covtheta[which(theta != 0), which(theta != 0)] <- theta.cov[which(theta != 0), which(theta != 0)]
}
Covtheta <- Matrix::bdiag(Covtheta,0)
Covtheta=as.matrix(Covtheta)

bZse_eff <- sweep(bZse, 2, vartheta, `*`)
G1=rowSums(bZse_eff*matrixMultiply(bZse_eff,Rxy))
G2=rowSums(bZ*matrixMultiply(bZ,Covtheta))
var_residual=G1+G2

# Result list
A <- list()
A$BETA <- residual
A$SE <- sqrt(var_residual)

return(A)
}
