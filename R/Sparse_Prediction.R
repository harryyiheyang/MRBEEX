#' @title Sparse Prediction of xQTL Effects using SuSiE
#' @description
#' This function performs informative xQTL selection and sparse prediction using SuSiE.
#'
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param xQTL.max.L The maximum number of single effects used to estimate the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres The minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to \code{0.95}.
#' @param unmappable_effects The method for modeling unmappable effects: "none", "inf".
#' @param estimate_residual_method The method used for estimating residual variance. For the original SuSiE model, "MLE" and "MoM" estimation is equivalent, but for the infinitesimal model, "MoM" is more stable. We recommend using "Servin_Stephens" when n < 80 for improved coverage.
#' @param xQTL.Nvec The vector of sample sizes of exposures.
#' @param xQTL.weight The vector of weights used in specifying the prior weights of SuSiE. Defaults to \code{NULL}.
#'
#' @importFrom MASS rlm ginv
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR coef.susie susie susie_rss susie_get_cs
#'
#' @return A list containing the results of the MRBEEX analysis using different methods:
#' \describe{
#'   \item{\code{bXest}}{A matrix where each column represents \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXestse}}{A matrix where each column contains the standard errors of \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXest0}}{A matrix where each column represents \code{beta_PLS} for each exposure.}
#'   \item{\code{bXestse0}}{A matrix where each column contains the standard errors of \code{beta_PLS} for each exposure.}
#'   \item{\code{xQTLfitList}}{A list containing the xQTL selection results, which can be used in \code{CisMRBEEX} by setting \code{xQTLfitList = xQTLfitList}.}
#' }
#' @export

Sparse_Prediction=function(bX,bXse,LD,xQTL.Nvec,
                           xQTL.max.L=10,xQTL.weight=NULL,xQTL.cred.thres=0.95,
                           estimate_residual_method="MoM",unmappable_effects="inf"){
p=ncol(bX)
m=nrow(bX)

xQTLfitList=list()
if(is.null(xQTL.weight)==T){
xQTL.weight=rep(1,m)
}
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000,prior_weights=xQTL.weight,estimate_residual_method=estimate_residual_method,unmappable_effects=unmappable_effects,convergence_method=ifelse(unmappable_effects=="none","elbo","pip"))
if(length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs)==xQTL.max.L){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L+5,max_iter=1000,prior_weights=xQTL.weight,estimate_residual_method=estimate_residual_method,unmappable_effects=unmappable_effects,convergence_method=ifelse(unmappable_effects=="none","elbo","pip"))
}else{
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=max(1,length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs)),max_iter=1000,prior_weights=xQTL.weight,estimate_residual_method=estimate_residual_method,unmappable_effects=unmappable_effects,convergence_method=ifelse(unmappable_effects=="none","elbo","pip"))
}
xQTLfitList[[i]]=fit
}

return(xQTLfitList)
}
