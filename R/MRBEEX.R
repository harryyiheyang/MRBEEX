#' @title Multivariable Mendelian Randomization using Bias-corrected Estimating Equations
#' @description
#' This function removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP) and correlated horizontal pleiotropy (CHP) through two distinct methods. UHP is detected using the IPOD algorithm, where outliers are interpreted as UHP. For CHP, a two-mixture regression model (Mixture) is implemented by the mixtools R package. Additionally, it integrates SuSiE for exposure selection, enhancing interpretability (use.susie=T). Both the IPOD algorithm and the Mixture method support the inclusion of correlated instrumental variables using an LD matrix and provide advanced options for exposure selection and horizontal pleiotropy correction.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix. Default is the identity matrix, assuming independent instrumental variables (IVs).
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param cluster.index A vector indicating the LD block indices each IV belongs to. The length is equal to the number of IVs, and values are the LD block indices.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param method Method for handling horizontal pleiotropy. Options are \code{"IPOD"} and \code{"Mixture"}.
#' @param use.susie An indicator of whether using SuSiE to select causal exposures. Defaults to \code{T}.
#' @param estimate_residual_method The method used for estimating residual variance. For the original SuSiE model, "MLE" and "MoM" estimation is equivalent, but for the infinitesimal model, "MoM" is more stable.
#' @param group.penalize An indicator of whether using difference penalty to penalize highly correlated exposures. Defaults to \code{F}.
#' @param group.index A vector of the group index of exposure. Defaults to \code{NULL}.
#' @param group.diff The tuning penalizing difference of highly correlated exposure prediction. Defaults to \code{100}.
#' @param main.cluster.thres When choosing \code{"Mixture"}, a threshold for weights belonging to the first category. To prevent instability caused by small-effect IVs falling into both categories, we slightly lower the voting threshold for the first category to below 0.5, ensuring it remains dominant. Default is \code{0.48}.
#' @param min.cluster.size When choosing \code{"Mixture"}, a minimum sample size of the second mixture.  Default is \code{5}.
#' @param tauvec When choosing \code{"IPOD"}, the candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(4, 8, by=0.5)}.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param Lvec When SuSiE is used, the candidate vector for the number of single effects. Default is \code{c(1:min(10, nrow(bX)))}.
#' @param pip.thres Posterior inclusion probability (PIP) threshold. Individual PIPs less than this value will be shrunk to zero. Default is \code{0.5}.
#' @param pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.1}.
#' @param estimate_residual_variance When SuSiE is used, an indicator of whether estimating the variance of residuals. If setting F, the variance of residual will be fixed as \code{Rxy[p+1,p+1]}. Default is \code{T}.
#' @param cred.pip.thres The threshold of PIP of each credible set. Defaults to \code{0.95}.
#' @param coverage.causal The coverage of defining a credible set in MRBEEX when \code{use.susie = T}. Defaults to \code{0.95}.
#' @param standardize If standardize = TRUE, standardize the columns of X to unit variance prior to fitting (or equivalently standardize XtX and Xty to have the same effect) in SuSiE. Note that scaled_prior_variance specifies the prior on the coefficients of X after standardization (if it is performed). If you do not standardize, you may need to think more carefully about specifying scaled_prior_variance. Whatever your choice, the coefficients returned by coef are given for X on the original input scale. Any column of X that has zero variance is not standardized.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param susie.iter Number of iterations in SuSiE per iteration. Default is \code{100}.
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to \code{3}.
#' @param ridge.diff A ridge.parameter on the differences of causal effect estimate in one credible set. Defaults to \code{1e3}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{0}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy. Default is \code{1}.
#' @param sampling.strategy "bootstrap" or "subsampling" (0.5 sample without replacement).
#' @param sampling.time Number of blockwise bootstrapping times. Default is \code{100}.
#' @param sampling.iter Number of iterations per blockwise bootstrapping procedure. Default is \code{10}.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param prob_shrinkage_coef Shrinkage coefficient (alpha) used to smooth block sampling probabilities within each group. 0 = no shrinkage; 1 = full shrinkage to group median.
#' @param prob_shrinkage_size Number of blocks per smoothing group (e.g., 3–5). Blocks are sorted by weight and grouped before applying shrinkage.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.
#' @param verbose A logical indicator of whether to display the execution time of the method. Default is \code{T}.
#' @param gcov A matrix (p+1 x p+1) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.
#'
#' @importFrom MASS rlm
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR susie_ss coef.susie susie
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
#'   \item{\code{gamma}}{Estimate of horizontal pleiotropy.}
#'   \item{\code{Bic}}{A vector or matrix recording the Bayesian Information Criterion (BIC) values.}
#'   \item{\code{theta.ini}}{Initial value of \code{theta} used in the estimation procedure.}
#'   \item{\code{gamma.ini}}{Initial value of \code{gamma} used in the estimation procedure.}
#'   \item{\code{reliability.adjust}}{Estimated reliability-adjusted values.}
#'   \item{\code{thetalist}}{List of \code{theta} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{gammalist}}{List of \code{gamma} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{theta1}}{Causal effect estimate for the first mixture component (when \code{method="Mixture"}).}
#'   \item{\code{theta2}}{Causal effect estimate for the second mixture component (when \code{method="Mixture"}).}
#'   \item{\code{theta.se1}}{Standard error of \code{theta1}.}
#'   \item{\code{theta.se2}}{Standard error of \code{theta2}.}
#'   \item{\code{theta.cov1}}{Covariance matrix of \code{theta1}.}
#'   \item{\code{theta.cov2}}{Covariance matrix of \code{theta2}.}
#'   \item{\code{theta.pip1}}{Empirical PIP of \code{theta1} in the subsampling procedure.}
#'   \item{\code{theta.pip2}}{Empirical PIP of \code{theta2} in the subsampling procedure.}
#'   \item{\code{thetalist1}}{List of \code{theta1} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{thetalist2}}{List of \code{theta2} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{Voting}}{A list that contains (1) the weights of two mixtures and (2) the voting results of two mixture based on main.cluster.thres.}
#' }
#' @export

MRBEEX=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),
               method=c("IPOD","Mixture"),use.susie=T,standardize=T,
               group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=100,
               main.cluster.thres=0.48,min.cluster.size=5,
               tauvec=seq(4,8,by=0.5),admm.rho=2,
               Lvec=c(1:min(10,ncol(bX))),pip.thres=0.5,estimate_residual_variance=T,
               pip.min=0.1,cred.pip.thres=0.95,
               max.iter=100,max.eps=0.001,susie.iter=100,
               ebic.theta=0,ebic.gamma=1,ridge.diff=1e3,
               estimate_residual_method="MoM",sampling.strategy="subsampling",
               sampling.time=300,sampling.iter=25,prob_shrinkage_coef=0.5,prob_shrinkage_size=4,
               maxdiff=3,reliability.thres=2/3,coverage.causal=0.95,
               theta.ini=F,gamma.ini=F,verbose=T,gcov=NULL,ldsc=NULL){

##########################################################################
cluster.index <- as.integer(factor(cluster.index))
if(method[1]=="IPOD"&use.susie==T){
A=MRBEE_IPOD_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Lvec=Lvec,pip.thres=pip.thres,pip.min=pip.min,cred.pip.thres=cred.pip.thres,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,susie.iter=susie.iter,ebic.theta=ebic.theta,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=admm.rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini,ridge.diff=ridge.diff,verbose=verbose,group.penalize=group.penalize,group.index=group.index,group.diff=group.diff,coverage.causal=coverage.causal,LDSC=ldsc,Omega=gcov,estimate_residual_variance=estimate_residual_variance,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,estimate_residual_method=estimate_residual_method,sampling.strategy=sampling.strategy,standardize=standardize)
}
##########################################################################
if(method[1]=="IPOD"&use.susie==F){
A=MRBEE_IPOD(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=admm.rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini,ebic.theta=ebic.theta,verbose=verbose,group.penalize=group.penalize,group.index=group.index,group.diff=group.diff,LDSC=ldsc,Omega=gcov,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,sampling.strategy=sampling.strategy)
}
##########################################################################
if(method[1]=="Mixture"&use.susie==F){
A=MRBEE_Mixture(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,reliability.thres=reliability.thres,sampling.time=sampling.time,ebic.theta=ebic.theta,max.iter=max.iter,max.eps=max.eps,sampling.iter=sampling.iter,main.cluster.thres=main.cluster.thres,min.cluster.size=min.cluster.size,group.penalize=group.penalize,group.index=group.index,group.diff=group.diff,LDSC=ldsc,Omega=gcov,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,sampling.strategy=sampling.strategy)
}
###########################################################################
if(method[1]=="Mixture"&use.susie==T){
A=MRBEE_Mixture_SuSiE(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Lvec=Lvec,pip.thres=pip.thres,pip.min=pip.min,cred.pip.thres=cred.pip.thres,ebic.theta=ebic.theta,reliability.thres=reliability.thres,sampling.time=sampling.time,max.iter=max.iter,max.eps=max.eps,sampling.iter=sampling.iter,susie.iter=susie.iter,main.cluster.thres=main.cluster.thres,min.cluster.size=min.cluster.size,ridge.diff=ridge.diff,group.penalize=group.penalize,group.index=group.index,group.diff=group.diff,coverage.causal=coverage.causal,LDSC=ldsc,Omega=gcov,estimate_residual_variance=estimate_residual_variance,prob_shrinkage_size=prob_shrinkage_size,prob_shrinkage_coef=prob_shrinkage_coef,estimate_residual_method=estimate_residual_method,sampling.strategy=sampling.strategy,standardize=standardize)
}

return(A)
}


