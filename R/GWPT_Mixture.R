#' @title Genome-Wide Pleiotropy Test for mixture model
#' @description
#' This function performs a genome-wide pleiotropy test (GWPT) after Mendelian randomization. It offers an option for a two-mixture model, where the residual is chosen as the smaller one resulting from the two causal effect estimates from two mixtures.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param theta1 The causal effect estimate of the first mixture.
#' @param theta.cov1 The covariance matrix of the causal effect estimate of the first mixture.
#' @param theta2 The causal effect estimate of the second mixture.
#' @param theta.cov2 The covariance matrix of the causal effect estimate of the second mixture.
#' @param LD.block A vector of indices of LD blocks.
#' @import data.table
#' @import CppMatrix
#'
#' @return A list with two components:
#' \item{BETA}{The estimated residual values.}
#' \item{SE}{The standard errors of the residual estimates.}
#'
#' @export

GWPT_Mixture <- function(by, byse, bX, bXse, Rxy, theta1, theta.cov1, theta2, theta.cov2, LD.block) {
bZ <- cbind(bX,by)
bZse <- cbind(bXse,byse)
vartheta1 <- c(-theta1,1)
vartheta2 <- c(-theta2,1)
residual1 <- matrixVectorMultiply(bZ,vartheta1)
residual2 <- matrixVectorMultiply(bZ,vartheta2)
var_residual1 <- var_residual2 <- residual1 * 0
Covtheta1 <- diag(length(vartheta1) - 1) * 0
Covtheta1[which(theta1 != 0), which(theta1 != 0)] <- theta.cov1[which(theta1 != 0), which(theta1 != 0)]
Covtheta1 <- Matrix::bdiag(Covtheta1,0)
Covtheta2 <- diag(length(vartheta2) - 1) * 0
Covtheta2[which(theta2 != 0), which(theta2 != 0)] <- theta.cov2[which(theta2 != 0), which(theta2 != 0)]
Covtheta2 <- Matrix::bdiag(Covtheta2,0)
Covtheta2=as.matrix(Covtheta2)
Covtheta1=as.matrix(Covtheta1)

bZse_eff1 <- sweep(bZse, 2, vartheta1, `*`)
G11=rowSums(bZse_eff1*matrixMultiply(bZse_eff1,Rxy))
G21=rowSums(bZ*matrixMultiply(bZ,Covtheta1))
var_residual1=G11+G21
bZse_eff2 <- sweep(bZse, 2, vartheta2, `*`)
G12=rowSums(bZse_eff2*matrixMultiply(bZse_eff2,Rxy))
G22=rowSums(bZ*matrixMultiply(bZ,Covtheta2))
var_residual2=G12+G22

dt <- data.table(
idx = seq_len(length(residual1)),
residual1 = residual1,
residual2 = residual2,
var_residual1 = var_residual1,
var_residual2 = var_residual2,
LD.block = LD.block
)
dt_block <- dt[, .(
count1 = sum(abs(residual1) < abs(residual2)),
count2 = sum(abs(residual1) >= abs(residual2))
), by = LD.block]
dt_block[, cluster := ifelse(count1 >= count2, 1, 2)]
dt <- merge(dt, dt_block[, .(LD.block, cluster)], by = "LD.block", all.x = TRUE)
dt[cluster == 1, `:=`(
residual = residual1,
var_residual = var_residual1
)]
dt[cluster == 2, `:=`(
residual = residual2,
var_residual = var_residual2
)]
setorder(dt, idx)

A=data.frame(SNP=rownames(bX),BETA=dt$residual,SE=sqrt(dt$var_residual),Cluster=dt$cluster,LD.block=LD.block)
return(A)
}
