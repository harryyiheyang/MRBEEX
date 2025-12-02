#' CisMRBEE_UV: Cis Univariable Mendelian Randomization using Bias-corrected Estimating Equations
#'
#' This function performs univariable cis-Mendelian randomization that removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP).
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A vector of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A vector of standard errors of effect estimates from the exposure GWAS.
#' @param LD The LD matrix of variants.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param xQTL.max.L The maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres The minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to \code{0.95}.
#' @param xQTL.pip.thres If SuSiE fails to find any credible set, the threshold of individual PIP when selecting xQTL. Defaults to \code{0.5}.
#' @param xQTL.pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.2}.
#' @param xQTL.N The sample sizes of exposure.
#' @param xQTL.selection.rule The method for purifying informative xQTLs within each credible set. Options include "minimum_pip", which selects all variables with PIPs exceeding a specified threshold, and "top_K", which ensures at least K variables are selected based on their PIP ranking. Defaults to "top_K".
#' @param top_K The maximum number of variables selected in each credible sets. Defaults to 1.
#' @param tauvec A vector of tuning parameters used in penalizing the direct causal effect. Default is `seq(3,10,by=1)`.
#' @param admm.rho A parameter set in the ADMM algorithm. Default is \code{2}.
#' @param max.iter The maximum number of iterations for the ADMM algorithm. Default is \code{15}.
#' @param max.eps The convergence tolerance for the ADMM algorithm. Default is \code{0.005}.
#' @param ebic.gamma The extended BIC factor for model selection. Default is \code{2}.
#' @param coverage.xQTL The coverage of defining a credible set in xQTL selection. Defaults to \code{0.95}.
#' @param coverage.causal The coverage of defining a credible set in cis-MRBEE. Defaults to \code{0.95}.
#' @param xQTLfit  Initial fits of xQTLs for exposures. This should only be yielded by SuSiE, as CARMA is not allowed for cis-UVMR analysis currently. Default is \code{NULL}.
#' @return A list containing:
#' \describe{
#' \item{\code{theta}}{The estimated effect size of the tissue-gene pair.}
#' \item{\code{gamma}}{The estimated effect sizes of the direct causal variants.}
#' \item{\code{theta.cov}}{The variance of the estimated effect size `theta`.}
#' \item{\code{theta.se}}{The standard error of the estimated effect size `theta`.}
#' \item{\code{theta.z}}{The z-score of the estimated effect size `theta`.}
#' \item{\code{Bic}}{The BIC values for each tuning parameter.}
#' \item{\code{eQTL.fit}}{The SuSiE result of xQLT selection of exposure.}
#' \item{\code{var.error}}{The variance of residuals.}
#' \item{\code{var.error}}{The variance of infinitesimal effect.}
#' \item{\code{causal.fit}}{The SuSiE result of causal effect calibration of exposure.}
#' \item{\code{reliability.adjust}}{Estimated reliability-adjusted values.}
#' }
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom susieR susie_ss coef.susie susie susie_rss susie_get_cs
#' @export
#'
CisMRBEE_UV=function(by,bX,byse,bXse,LD,Rxy,xQTL.N,xQTL.selection.rule="top_K",
     top_K=1,xQTL.pip.min=0.2,
     xQTL.max.L=10,xQTL.cred.thres=0.95,
     xQTL.pip.thres=0.5,reliability.thres=0.75,
     tauvec=seq(3,30,by=1.5),admm.rho=2,
     coverage.xQTL=0.95,coverage.causal=0.95,
     max.iter=100,max.eps=0.001,ebic.gamma=2,
     xQTLfit=NULL){
m=length(by)
Theta=solve(LD)
bXest=bX
bXest0=bXestse0=bX*0
bXestse=bXestse0=c(1000,m)
if(is.null(xQTLfit)==T){
fit.susie=susie_rss(z=bX/bXse,R=LD,n=xQTL.N,L=xQTL.max.L,max_iter=1000,coverage=coverage.xQTL)
fit.susie=susie_rss(z=bX/bXse,R=LD,n=xQTL.N,L=length(susie_get_cs(fit.susie,coverage=xQTL.cred.thres)$cs)+1,max_iter=1000,coverage=coverage.xQTL)
xQTLfit=fit.susie
causal.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(summary(fit.susie)$vars,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
causal.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
indj=union(causal.cs$ind.keep,which(fit.susie$pip>xQTL.pip.thres))
}
if(length(indj)>0){
betaj=coef.susie(fit.susie)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit.susie)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj]=as.vector(Thetaj%*%(bX[indj]/bXse[indj]))*bXse[indj]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0=sqrt(diag(Thetajj))*bXse
LDjj=LD%*%Thetajj%*%LD
bXest=as.vector(LD%*%(bXest0/bXse))*bXse
bXestse=sqrt(diag(LDjj))*bXse
}else{
# If SuSiE cannot find a credible set, use ridge regression with tuning parameter = 0.5 instead.
Theta_ridge=matrixInverse(LD+0.5*diag(diag(LD)))
bXest0=c(Theta_ridge%*%(bX/bXse))*bXse
bXestse0=sqrt(diag(Theta_ridge))*bXse
bXest=bX
bXestse=bXse
}
}else{
fit.susie=xQTLfit
causal.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
if(xQTL.selection.rule=="top_K"){
indj=top_K_pip(summary(fit.susie)$vars,top_K=top_K,pip.min.thres=xQTL.pip.min,xQTL.pip.thres=xQTL.pip.thres)
}else{
causal.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=xQTL.cred.thres,xQTL.pip.thres=xQTL.pip.min)
indj=union(causal.cs$ind.keep,which(fit.susie$pip>xQTL.pip.thres))
}
if(length(indj)>0){
betaj=coef.susie(fit.susie)[-1]
betaj[-indj]=0
Diff=generate_block_matrix(summary(fit.susie)$vars,rep(1,m),betaj)
LDj=LD[indj,indj]
Thetaj=solve(LDj+Diff[indj,indj]*1e3)
bXest0[indj]=as.vector(Thetaj%*%(bX[indj]/bXse[indj]))*bXse[indj]
Thetajj=LD*0
Thetajj[indj,indj]=Thetaj
bXestse0=sqrt(diag(Thetajj))*bXse
LDjj=LD%*%Thetajj%*%LD
bXest=as.vector(LD%*%(bXest0/bXse))*bXse
bXestse=sqrt(diag(LDjj))*bXse
}else{
# If SuSiE cannot find a credible set, use ridge regression with tuning parameter = 1 instead.
bXest0=c(Theta%*%(bX/bXse))*bXse
bXestse0=sqrt(diag(Theta))*bXse
bXest=bX
bXestse=bXse
}
}
r=reliability.adj.uv(bXest,bXestse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r

by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
bXest=bXest*byseinv
bXestse=bXestse*byseinv
bXest0=bXest0*byseinv
bXestse0=bXestse0*byseinv
byse1=byse
byse=byse/byse

if(sum(bXest0!=0)==1){
pleiotropy.rm=which(bXest0!=0)
}else{
pleiotropy.rm=NULL
}
pleiotropy.keep=setdiff(c(1:m),pleiotropy.rm)
yinv=c(matrixVectorMultiply(Theta,by))
xtx=sum(bXest*bXest0)
xty=sum(by*bXest0)
theta.ini=xty/xtx
Thetarho=matrixInverse(LD[pleiotropy.keep,pleiotropy.keep]+admm.rho*diag(length(pleiotropy.keep)))
tauvec=sort(tauvec,decreasing=F)
w=length(tauvec)
Btheta=c(1:w)
Bgamma=matrix(0,m,w)
Bbic=c(1:w)
for(sss in c(w:1)){
theta=theta.ini
gamma=by*0
gamma1=gamma
u=admm.rho*(gamma-gamma1)
theta1=theta*0
error=1
iter=1
fit.theta=NULL
while(error>max.eps&iter<max.iter){
indvalid=which(gamma1==0)
theta1=theta
Hinv=1/(xtx-sum(bXestse[indvalid]^2)*Rxy[1,1])
g=sum(bXest0*(by-as.vector(LD%*%gamma)))-sum(bXestse[indvalid])*Rxy[2,1]
theta=g*Hinv
res=c(by-bXest*theta-u+admm.rho*gamma1)
gamma[pleiotropy.keep]=c(matrixVectorMultiply(Thetarho,res[pleiotropy.keep]))
gamma1=mcp(gamma+u/admm.rho,tauvec[sss])
u=u+admm.rho*(gamma-gamma1)
u[pleiotropy.rm]=0
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>3){
error=abs(theta-theta1)
}
Btheta[sss]=theta
Bgamma[,sss]=gamma
res=as.vector(by-bXest*theta-matrixVectorMultiply(LD,gamma))
df=sum(gamma!=0)
rss=sum(res*(matrixVectorMultiply(Theta,res)))
Bbic[sss]=m*log(rss)+log(m)*(1+ebic.gamma)*df
}
}

star=which.min(Bbic)
theta=Btheta[star]
gamma=Bgamma[,star]
indvalid=which(gamma==0)
indgamma=which(gamma!=0)
effn=m-length(indgamma)

res=as.vector(by-bXest*theta-LD%*%gamma)
upsilon=by*0
var_inf=1
var_error=1
for(vv in 1:30){
Hupsilon=solve(diag(m)/var_inf+LD/var_error)
upsilon=as.vector(Hupsilon%*%res)/var_error
df=sum(diag(Hupsilon))
var_inf=min((sum(upsilon^2)+df)/m,10)
res_inf=res-matrixVectorMultiply(LD,upsilon)
df=sum(diag(Hupsilon%*%LD))/var_error
var_error=sum(res_inf*(Theta%*%res_inf))/(m-df-1-length(indgamma))
var_error=max(1,var_error)
}

if(sum(indgamma)>0){
Z=cbind(bXest,LD[,indgamma])
Hinv=matrixMultiply(t(Z),matrixMultiply(Theta,Z))
Hinv[1,1]=Hinv[1,1]-sum(bXestse[indvalid])*Rxy[1,1]
Hinv=positiveinv(Hinv)
Hinv1=matrixListProduct(list(t(Z),Theta,var_error*LD+var_inf*LD%*%LD,Theta,Z))
COV=Hinv%*%Hinv1%*%Hinv
covg=COV[1,1]
covtheta=covg
}

if(sum(indgamma)==0){
h0=sum(bXest0*((var_error*LD+var_inf*LD%*%LD)%*%bXest0))
h1=(xtx-sum(bXestse[indvalid])*Rxy[1,1])
covtheta=h0/h1/h1
}

A=list()
A$theta=theta
A$gamma=gamma
A$theta.cov=as.numeric(covtheta)
A$theta.se=sqrt(A$theta.cov)
A$theta.z=A$theta/A$theta.se
A$Bic=Bbic
A$eQTL.fit=fit.susie
A$causal.fit=fit.theta
A$reliability.adjust=r
A$direct_exposure_effect=bXest0
return(A)
}
