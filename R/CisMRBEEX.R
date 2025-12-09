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
#' @param xQTL.method The method used in purifying the xQTLs. SuSiE or CARMA can be used here, where the latter can be more accurate but much most computationally costly. Defaults is SuSiE.
#' @param xQTL.selection.rule The method for purifying informative xQTLs within each credible set. Options include "minimum_pip", which selects all variables with PIPs exceeding a specified threshold, and "top_K", which ensures at least K variables are selected based on their PIP ranking. Defaults to "top_K".
#' @param top_K The maximum number of variables selected in each credible sets. Defaults to 1.
#' @param xQTL.pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.2}.
#' @param xQTL.pip.thres The threshold of individual PIP when selecting xQTL. Defaults to \code{0.5}.
#' @param xQTL.max.L When choosing \code{"SuSiE"}, the maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres When choosing \code{"SuSiE"}, the minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to \code{0.95}.
#' @param xQTL.Nvec When choosing \code{"SuSiE"}, the vector of sample sizes of exposures.
#' @param xQTL.weight When choosing \code{"SuSiE"}, the vector of weights used in specifying the prior weights of SuSiE. Defaults to \code{NULL}.
#' @param coverage.xQTL The coverage of defining a credible set in xQTL selection. Defaults to \code{0.95}.
#' @param coverage.causal The coverage of defining a credible set in cis-MRBEE. Defaults to \code{0.95}.
#' @param outlier.switch When choosing \code{"CARMA"}, an indicator of whether turning on outlier detection. Defaults to \code{F}.
#' @param Annotation When choosing \code{"CARMA"}, the annotation matrix of SNP. Default is NULL.
#' @param output.labels When choosing \code{"CARMA"}, output directory where output will be written while CARMA is running. Defaults to \code{NULL}, meaning that a temporary folder will be created and automatically deleted upon completion of the computation.
#' @param carma.iter When choosing \code{"CARMA"}, the maximum iterations for EM algorithm to run. Defaults to 5.
#' @param carma.inner.iter When choosing \code{"CARMA"}, the maximum iterations for Shotgun algorithm to run per iteration within EM algorithm. Defaults to 5.
#' @param xQTL.max.num When choosing \code{"CARMA"}, the maximum number of causal variants assumed per locus, which is similar to the number of single effects in SuSiE. Defaults to 10.
#' @param carma.epsilon.threshold When choosing \code{"CARMA"}, the convergence threshold measured by average of Bayes factors. Defaults to \code{1e-3}.
#' @param model.infinitesimal An indicator of whether using REML to model infinitesimal effects. Defaults to \code{F}.
#' @param ridge.diff A ridge.parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param tauvec The candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param admm.rho The tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param Lvec When SuSiE is used, the candidate vector for the number of single effects. Default is \code{c(1:min(10, nrow(bX)))}.
#' @param causal.pip.thres A threshold of minimum posterior inclusion probability. Default is \code{0.2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param susie.iter Number of iterations in SuSiE per iteration. Default is \code{500}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{1}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy Default is \code{2}.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.
#' @param xQTLfitList Initial fits of xQTLs for exposures. This should be a list. Each component corresponds to the susie.fit of each exposure when xQTL.method = "SuSiE". When xQTL.method = "CARMA", this should be the list of results from a CARMA analysis. Users can customize additional SuSiE or CARMA parameters to improve performance. Default is \code{NULL}.
#' @param verbose A logical indicator of whether to display the execution time of the method. Default is \code{T}.
#'
#' @importFrom MASS rlm ginv
#' @importFrom varbvs varbvs
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct matrixGeneralizedInverse
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR susie_ss coef.susie susie susie_rss susie_get_cs
#'
#' @return A list that contains the results of the MRBEEX with respect to different methods applied:
#' \describe{
#'   \item{\code{theta}}{Causal effect estimate.}
#'   \item{\code{theta.se}}{Standard error of the causal effect estimate.}
#'   \item{\code{theta.cov}}{Covariance matrix of the causal effect estimate.}
#'   \item{\code{theta.pip}}{Empirical posterior inclusion probability (PIP) of the causal effect in the subsampling procedure.}
#'   \item{\code{theta.pratt}}{Pratt index estimate of exposure.}
#'   \item{\code{susie.theta}}{The fit of causal effect resulted from susie.}
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
#'   \item{estimated.reliability.ratio}{The estimated reliability ratios of exposures.}
#'   \item{xQTLfitList}{The results of sparse predictions of exposures, yielded by SuSiE or CARMA.}
#' }
#' @export

CisMRBEEX=function(by,bX,byse,bXse,LD,Rxy,model.infinitesimal=F,
                    reliability.thres=0.75,Lvec=c(1:5),causal.pip.thres=0.2,
                    xQTL.method="SuSiE",xQTL.selection.rule="top_K",
                    top_K=1,xQTL.pip.min=0.2,
                    xQTL.max.L=10,xQTL.cred.thres=0.95,xQTL.pip.thres=0.5,
                    xQTL.Nvec,tauvec=seq(3,30,by=3),xQTL.weight=NULL,
                    outlier.switch=T,Annotation=NULL,output.labels=NULL,
                    carma.iter=5,carma.inner.iter=5,xQTL.max.num=10,
                    carma.epsilon.threshold=1e-3,
                    admm.rho=2,ridge.diff=1e3,
                    max.iter=100,max.eps=0.001,susie.iter=500,
                    coverage.xQTL=0.95,coverage.causal=0.95,
                    ebic.theta=0,ebic.gamma=1,
                    theta.ini=F,gamma.ini=F,xQTLfitList=NULL,
                    verbose=T){

cat("Please standardize data such that BETA = Zscore/sqrt n and SE = 1/sqrt n\n")
######################### Estimate xQTL effect size ############################
p=ncol(bX)
m=nrow(bX)
bXest=bX
bXest0=bXestse0=bX*0
bXestse=bXestse0=matrix(1000,m,p)

t1=Sys.time()
A=list()
estimated_reliability_ratio=c()
if(xQTL.method=="SuSiE"){
if(is.null(xQTLfitList)==T){
xQTLfitList=list()
if(is.null(xQTL.weight)==T){
xQTL.weight=rep(1,m)
}
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000,prior_weights=xQTL.weight,coverage=coverage.xQTL)
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs)+1,max_iter=1000,prior_weights=xQTL.weight,coverage=coverage.xQTL)
xQTLfitList[[i]]=fit
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(summary(fit)$vars,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
causal.cs=group.pip.filter(pip.summary=summary(fit)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
indj=union(causal.cs$ind.keep,which(fit$pip>xQTL.pip.thres))
}
if(length(indj)>0){
betaj=coef.susie(fit)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
ii=which.max(abs(betaj))
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
estimated_reliability_ratio[i]=(sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))-sum(diag(Thetaj)*bXse[indj,i]^2))/sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
estimated_reliability_ratio[i]=c(sum(bX[,i]^2)-sum(bXse[,i]^2))/sum(bX[,i]^2)
}
}
}else{
for(i in 1:p){
fit=xQTLfitList[[i]]
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(summary(fit)$vars,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
causal.cs=group.pip.filter(pip.summary=summary(fit)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
indj=union(causal.cs$ind.keep,which(fit$pip>xQTL.pip.thres))
}
if(length(indj)>0){
betaj=coef.susie(fit)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
ii=which.max(abs(betaj))
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
estimated_reliability_ratio[i]=(sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))-sum(diag(Thetaj)*bXse[indj,i]^2))/sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
estimated_reliability_ratio[i]=c(sum(bX[,i]^2)-sum(bXse[,i]^2))/sum(bX[,i]^2)
}
}
}
}
if(xQTL.method=="CARMA"){
if(is.null(xQTLfitList)==T){
if (!requireNamespace("CARMA", quietly = TRUE)) {
stop("Package 'CARMA' is required for xQTL.method='CARMA', but is not installed. ",
"Please install it with install.packages('CARMA').")
}
z.list=ld.list=w.list=lambda.list=list()
for(i in 1:p){
z.list[[i]]=bX[,i]/bXse[,i]
ld.list[[i]]=LD
w.list[[i]]=Annotation
lambda.list[[i]]=1
}
is.detect=F
if(is.null(output.labels)==T){
output.labels=tempfile("my_tmpdir_")
dir.create(output.labels)
is.delect=T
}
if(is.null(Annotation)==F){
fitxQTL=CARMA::CARMA(z.list,ld.list,w.list,lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}else{
fitxQTL=CARMA::CARMA(z.list,ld.list,lambda.list = lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}
if(is.delect==T){
unlink(output.labels, recursive = TRUE, force = TRUE)
}
for(i in 1:p){
sumstat.result = data.frame(variable=c(1:nrow(bX)),variable_prob = fitxQTL[[i]]$PIPs, cs = rep(0,nrow(bX)))
if(length(fitxQTL[[i]]$`Credible set`[[2]])!=0){
for(l in 1:length(fitxQTL[[i]]$`Credible set`[[2]])){
sumstat.result$cs[fitxQTL[[i]]$`Credible set`[[2]][[l]]]=l
}
}
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(sumstat.result,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
indj=which(sumstat.result$cs!=0&sumstat.result$pip>xQTL.pip.min)
}
if(length(indj)>0){
fitjj=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=1+sum(sumstat.result$cs>0),max_iter=100)
betaj=coef.susie(fitjj)[-1]
Diff=generate_block_matrix_CARMA(sumstat.result,rep(1,m),betaj)
Diff=Diff[indj,indj]
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
estimated_reliability_ratio[i]=(sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))-sum(diag(Thetaj)*bXse[indj,i]^2))/sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
estimated_reliability_ratio[i]=c(sum(bX[,i]^2)-sum(bXse[,i]^2))/sum(bX[,i]^2)
}
fitxQTL[[i]]$Summary=sumstat.result
}
}else{
fitxQTL=xQTLfitList
for(i in 1:p){
sumstat.result = data.frame(variable=c(1:nrow(bX)),variable_prob = fitxQTL[[i]]$PIPs, cs = rep(0,nrow(bX)))
if(length(fitxQTL[[i]]$`Credible set`[[2]])!=0){
for(l in 1:length(fitxQTL[[i]]$`Credible set`[[2]])){
sumstat.result$cs[fitxQTL[[i]]$`Credible set`[[2]][[l]]]=l
}
}
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(sumstat.result,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
indj=which(sumstat.result$cs!=0&sumstat.result$pip>xQTL.pip.min)
}
if(length(indj)>0){
fitjj=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=1+sum(sumstat.result$cs>0),max_iter=100,coverage=coverage.xQTL)
betaj=coef.susie(fitjj)[-1]
Diff=generate_block_matrix_CARMA(sumstat.result,rep(1,m),betaj)
Diff=Diff[indj,indj]
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff*1e3)
bXest0[indj,i]=as.vector(Thetaj%*%(bX[indj,i]/bXse[indj,i]))*bXse[indj,i]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0[,i]=sqrt(diag(Thetajj))*bXse[,i]
LDjj=LD%*%Thetajj%*%LD
bXest[,i]=as.vector(LD%*%(bXest0[,i]/bXse[,i]))*bXse[,i]
bXestse[,i]=sqrt(diag(LDjj))*bXse[,i]
estimated_reliability_ratio[i]=(sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))-sum(diag(Thetaj)*bXse[indj,i]^2))/sum(bXest0[indj,i]*(LDj%*%bXest0[indj,i]))
}else{
bXest[,i]=bX[,i]
bXestse[,i]=bXse[,i]
estimated_reliability_ratio[i]=c(sum(bX[,i]^2)-sum(bXse[,i]^2))/sum(bX[,i]^2)
}
}
}
}
t2=Sys.time()
sparse_prediction_time=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Sparse prediction ends: ",sparse_prediction_time," secs\n"))
}
pleiotropy.rm=findUniqueNonZeroRows(bXest0)
##########################################################################
if(model.infinitesimal==F){
t1=Sys.time()
A=Cis_MRBEE_IPOD_SuSiE(by=by,bX=bXest,byse=byse,bXse=bXestse,LD=LD,Rxy=Rxy,pip.thres=causal.pip.thres,Lvec=Lvec,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,susie.iter=susie.iter,ebic.theta=ebic.theta,ebic.gamma=ebic.gamma,reliability.thres=reliability.thres,rho=admm.rho,theta.ini=theta.ini,gamma.ini=gamma.ini,ridge=ridge.diff,pleiotropy.rm=pleiotropy.rm,coverage.causal=coverage.causal,xQTLfitList=xQTLfitList)
t2=Sys.time()
causal_estimation_time=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Causal effect estimation ends: ",causal_estimation_time," secs\n"))
}
}
##########################################################################
if(model.infinitesimal==T){
cat("Awaiting development\n")
}
A$bXest=bXest
A$bXestse=bXestse
A$bXest0=bXest0
A$bXestse0=bXestse0
A$estimated_reliability_ratio=estimated_reliability_ratio
if(xQTL.method=="SuSiE"){
A$xQTLfitList=xQTLfitList
}else{
A$xQTLfitList=fitxQTL
}
return(A)
}


