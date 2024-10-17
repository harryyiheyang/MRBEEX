#' @title Cis Multivariable Mendelian Randomization using Bias-corrected Estimating Equations
#' @description
#' This function performs multivariable cis-Mendelian randomization that removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP). Additionally, it integrates SuSiE for exposure selection, enhancing interpretability (use.susie=T).
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param xQTL.max.L The maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.sampling The number of subsampling times in estimating the standard errors of xQTL effects. Defaults to 1000.
#' @param xQTL.pip.thres The minimum empirical posterior inclustion probability (PIP) used in calibrating causality for xQTL data.
#' @param xQTL.Nvec The vector of sample sizes of exposures.
#' @param use.susie An indicator of whether using SuSiE to select causal exposures. Defaults to \code{T}.
#' @param estimate_residual_variance An indicator of whether estimating the residual_variance in SuSiE. Defaults to \code{T}.
#' @param ridge A ridge.parameter on causal effect estimate. Defaults to \code{0.05}.
#' @param residual_variance When \code{estimate_residual_variance = T}, the initial value of the residual variance. When \code{estimate_residual_variance = F}, the fixed value of the residual variance. Defaults to \code{1}.
#' @param tauvec When choosing \code{"IPOD"}, the candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param Lvec When SuSiE is used, the candidate vector for the number of single effects. Default is \code{c(1:min(10, nrow(bX)))}.
#' @param pip.thres A threshold of minimum posterior inclusion probability. Default is \code{0.2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param susie.iter Number of iterations in SuSiE per iteration. Default is \code{100}.
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to \code{1.5}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{1}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy Default is \code{2}.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.

#'
#' @importFrom MASS rlm ginv
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom varbvs varbvs
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR susie_suff_stat coef.susie susie susie_rss
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mixtools regmixEM
#' @importFrom FDRestimation p.fdr
#'
#' @return A list that contains the results of the MRBEEX with respect to different methods applied:
#' \describe{
#'   \item{\code{theta}}{Causal effect estimate.}
#'   \item{\code{theta.se}}{Standard error of the causal effect estimate.}
#'   \item{\code{theta.cov}}{Covariance matrix of the causal effect estimate.}
#'   \item{\code{theta.pip}}{Empirical posterior inclusion probability (PIP) of the causal effect in the subsampling procedure.}
#'   \item{\code{theta.pratt}}{Pratt index estimate of exposure.}
#'   \item{\code{gamma}}{Estimate of horizontal pleiotropy.}
#'   \item{\code{gamma.pratt}}{Pratt index estimate of horizontal pleiotropy.}
#'   \item{\code{Bic}}{A vector or matrix recording the Bayesian Information Criterion (BIC) values.}
#'   \item{\code{theta.ini}}{Initial value of \code{theta} used in the estimation procedure.}
#'   \item{\code{gamma.ini}}{Initial value of \code{gamma} used in the estimation procedure.}
#'   \item{\code{reliability.adjust}}{Estimated reliability-adjusted values.}
#'   \item{\code{thetalist}}{List of \code{theta} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{gammalist}}{List of \code{gamma} estimates recorded during each iteration in the subsampling procedure.}
#' }
#' @export

CisMRBEEX=function(by,bX,byse,bXse,LD,Rxy,use.susie=T,
                    estimate_residual_variance=T,residual_variance=1,
                    reliability.thres=0.9,Lvec=c(1:5),pip.thres=0.2,
                    xQTL.max.L=10,xQTL.sampling=1000,
                    xQTL.pip.thres=0.2,xQTL.Nvec,
                    tauvec=seq(3,30,by=3),rho=2,ridge=0.05,
                    max.iter=100,max.eps=0.001,susie.iter=100,
                    ebic.theta=1,ebic.gamma=2,maxdiff=3,
                    theta.ini=F,gamma.ini=F){

cat("Note: susie_rss() will be used to estimate the xQTL effect sizes\n")
cat("Please standardize data such that BETA = Zscore/sqrt n and SE = 1/sqrt n\n")
######################### Estimate xQTL effect size ############################
p=ncol(bX)
m=nrow(bX)
bXest=bX
bXest0=bXestse0=bX*0
bXestse=bXestse0=matrix(1000,m,p)
nonzero=c(1:p)
fitList=list()
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000)
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=sum(fit$pip>xQTL.pip.thres)+1,max_iter=1000)
fit=susie_xQTL_resampling(LD=LD,alpha=fit$alpha,mu=fit$mu,mu2=fit$mu2,pip.thres=xQTL.pip.thres,sampling=xQTL.sampling)
if(sum(abs(fit$mean_margin))>0){
bXest[,i]=fit$mean_margin
bXestse[,i]=fit$sd_margin
bXest0[,i]=fit$mean_joint
bXestse0[,i]=fit$sd_joint
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
}
}
##########################################################################
if(use.susie==T){
A=Cis_MRBEE_IPOD_SuSiE(by=by,bX=bXest,byse=byse,bXse=bXestse,LD=LD,Rxy=Rxy,pip.thres=pip.thres,Lvec=Lvec,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,susie.iter=susie.iter,ebic.theta=ebic.theta,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=rho,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,estimate_residual_variance=estimate_residual_variance,residual_variance=residual_variance,ridge=ridge)
}
##########################################################################
if(use.susie==F){
A=Cis_MRBEE_IPOD(by=by,bX=bXest,byse=byse,bXse=bXestse,LD=LD,Rxy=Rxy,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=rho,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,ebic.theta=ebic.theta,ridge=ridge)
}
A$bXest=bXest
A$bXestse=bXestse
A$bXest0=bXest0
A$bXestse0=bXestse0
return(A)
}


