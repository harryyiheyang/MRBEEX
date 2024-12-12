#' @title Cis Multivariable Mendelian Randomization using Bias-corrected Estimating Equations
#' @description
#' This function performs multivariable cis-Mendelian randomization that removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP). Additionally, it integrates SuSiE for exposure selection, enhancing interpretability.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param xQTL.max.L The maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres The minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to 0.95
#' @param xQTL.pip.thres The minimum empirical PIP used in purifying variables in each credible set. Defaults to 0.3
#' @param xQTL.Nvec The vector of sample sizes of exposures.
#' @param model.infinitesimal An indicator of whether using REML to model infinitesimal effects. Defaults to \code{F}.
#' @param ridge.diff A ridge.parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param tauvec When choosing \code{"IPOD"}, the candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param Lvec When SuSiE is used, the candidate vector for the number of single effects. Default is \code{c(1:min(10, nrow(bX)))}.
#' @param pip.thres A threshold of minimum posterior inclusion probability. Default is \code{0.2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param susie.iter Number of iterations in SuSiE per iteration. Default is \code{500}.
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to \code{1.5}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{1}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy Default is \code{2}.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.
#' @param eQTLfitList Initial fits of xQTLs of exposures. Default is \code{NULL}.
#'
#' @importFrom MASS rlm ginv
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom varbvs varbvs
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR susie_suff_stat coef.susie susie susie_rss susie_get_cs
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
#'   \item{var.error}{The variance of residuals.}
#'   \item{var.error}{The variance of infinitesimal effect.}
#' }
#' @export

CisMRBEEX=function(by,bX,byse,bXse,LD,Rxy,model.infinitesimal=F,
                 reliability.thres=0.75,Lvec=c(1:5),pip.thres=0.2,
                 xQTL.max.L=10,xQTL.cred.thres=0.2,xQTL.pip.thres=0.3,
                 xQTL.Nvec,tauvec=seq(3,30,by=3),admm.rho=2,ridge.diff=1e3,
                 max.iter=100,max.eps=0.001,susie.iter=500,
                 ebic.theta=1,ebic.gamma=2,maxdiff=3,
                 theta.ini=F,gamma.ini=F,eQTLfitList=NULL){

cat("Note: susie_rss() will be used to estimate the xQTL effect sizes\n")
cat("Please standardize data such that BETA = Zscore/sqrt n and SE = 1/sqrt n\n")
######################### Estimate xQTL effect size ############################
p=ncol(bX)
m=nrow(bX)
bXest=bX
bXest0=bXestse0=bX*0
bXestse=bXestse0=matrix(1000,m,p)
if(is.null(eQTLfitList)==T){
eQTLfitList=list()
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000)
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs)+1,max_iter=1000)
eQTLfitList[[i]]=fit
causal.cs=group.pip.filter(pip.summary=summary(fit)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.thres)
indj=causal.cs$ind.keep
if(length(indj)>0){
betaj=coef.susie(fit)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
}
}
}else{
for(i in 1:p){
fit=eQTLfitList[[i]]
causal.cs=group.pip.filter(pip.summary=summary(fit)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.thres)
indj=causal.cs$ind.keep
if(length(indj)>0){
betaj=coef.susie(fit)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
}
}
}
pleiotropy.rm=findUniqueNonZeroRows(bXest0)
##########################################################################
if(model.infinitesimal==F){
A=Cis_MRBEE_IPOD_SuSiE(by=by,bX=bXest,byse=byse,bXse=bXestse,LD=LD,Rxy=Rxy,pip.thres=pip.thres,Lvec=Lvec,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,susie.iter=susie.iter,ebic.theta=ebic.theta,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=admm.rho,maxdiff=maxdiff,theta.ini=theta.ini,gamma.ini=gamma.ini,ridge=ridge.diff,pleiotropy.rm=pleiotropy.rm)
}
##########################################################################
if(model.infinitesimal==T){
cat("awaiting development\n")}
A$bXest=bXest
A$bXestse=bXestse
A$bXest0=bXest0
A$bXestse0=bXestse0
A$eQTLfitList=eQTLfitList
return(A)
}


