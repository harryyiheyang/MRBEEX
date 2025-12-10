#' @title Univariable Mendelian Randomization using Bias-corrected Estimating Equations
#' @description
#' This function removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP) and correlated horizontal pleiotropy (CHP) through IPOD or greedy search.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A vector of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A vector of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix. Default is the identity matrix, assuming independent instrumental variables (IVs).
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last element corresponds to the outcome.
#' @param cluster.index A vector indicating the LD block indices each IV belongs to. The length is equal to the number of IVs, and values are the LD block indices.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param Method Method for handling horizontal pleiotropy. Options are \code{"IPOD"} and \code{"Greedy"}.
#' @param tauvec When choosing \code{"IPOD"}, the candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(4, 8, by=0.5)}.
#' @param Kvec When choosing \code{"Greedy"}, the candidate vector of number of pleiotropy in the greedy search algorithm. Default is \code{seq(1, length(bX)/2, by=2)}.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to \code{1.5}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy Default is \code{0}.
#' @param sampling.strategy "bootstrap" or "subsampling" (0.5 sample without replacement).
#' @param sampling.time Number of resampling times. Default is \code{100}.
#' @param sampling.iter Number of iterations per resampling. Default is \code{20}.
#' @param prob_shrinkage_coef Shrinkage coefficient (alpha) used to smooth block sampling probabilities within each group. 0 = no shrinkage; 1 = full shrinkage to group median.
#' @param prob_shrinkage_size Number of blocks per smoothing group (e.g., 3–5). Blocks are sorted by weight and grouped before applying shrinkage.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.
#' @param gcov A matrix (p+1 x p+1) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.

#' @importFrom MASS rlm
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom MRBEE MRBEE.IMRP.UV
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mixtools regmixEM
#' @importFrom FDRestimation p.fdr
#'
#' @return A list containing the results of the MRBEEX.UV analysis:
#' \describe{
#'   \item{\code{theta}}{Causal effect estimate.}
#'   \item{\code{theta.se}}{Standard error of the causal effect estimate.}
#'   \item{\code{theta.cov}}{Covariance matrix of the causal effect estimate.}
#'   \item{\code{theta.bootstrap}}{Resampled causal effect estimates in bootstrap.}
#'   \item{\code{gamma}}{Estimate of horizontal pleiotropy.}
#'   \item{\code{Bic}}{A vector or matrix recording the Bayesian Information Criterion (BIC) values.}
#'   \item{\code{theta.ini}}{Initial value of \code{theta} used in the estimation procedure.}
#'   \item{\code{gamma.ini}}{Initial value of \code{gamma} used in the estimation procedure.}
#'   \item{\code{reliability.adjust}}{Estimated reliability-adjusted values.}
#'   \item{\code{thetalist}}{List of \code{theta} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{gammalist}}{List of \code{gamma} estimates recorded during each iteration in the subsampling procedure.}
#' }
#'
#' @export

MRBEEX_UV=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),
        reliability.thres=0.8,
        Method="IPOD",
        tauvec=seq(4,8,by=2),admm.rho=2,ebic.gamma=0,
        Kvec=seq(1,length(bX)/2,by=2),
        max.iter=100,max.eps=0.001,maxdiff=3,sampling.strategy="subsampling",
        sampling.time=1000,sampling.iter=20,prob_shrinkage_coef=1,prob_shrinkage_size=4,
        theta.ini=F,gamma.ini=F,ldsc=NULL,gcov=NULL){
##########################################################################
if(Method[1]=="IPOD"){
A=MRBEE_IPOD_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,ebic.gamma=ebic.gamma,max.iter=max.iter,max.eps=max.eps,tauvec=tauvec,reliability.thres=reliability.thres,rho=admm.rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini,LDSC=ldsc,Omega=gcov,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,sampling.strategy=sampling.strategy)
}
##########################################################################
if(Method[1]=="Greedy"){
A=MRBEE_Greedy_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Kvec=Kvec,reliability.thres=reliability.thres,sampling.time=sampling.time,max.iter=max.iter,max.eps=max.eps,sampling.iter=sampling.iter,LDSC=ldsc,Omega=gcov,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,sampling.strategy=sampling.strategy)
}

return(A)
}


