#' Estimate Non-Transferable Causal Effect with MRBEE and SuSiE in UVMR
#'
#' This function estimates the non-transferable causal effect using a bias-correction estimating equation in a univariable MR model, considering potential pleiotropy and measurement errors, and using SuSiE to select the non-transferable causal effect.
#'
#' @param by A vector (n x 1) of the GWAS effect size of outcome.
#' @param bX A matrix (n x p) of the GWAS effect sizes of p exposures.
#' @param byse A vector (n x 1) of the GWAS effect size SE of outcome.
#' @param bXse A matrix (n x p) of the GWAS effect size SEs of p exposures.
#' @param theta.source A vector (p x 1) of the causal effect estimate learning from the source data.
#' @param theta.source.cov A matrix (p x p) of the covariance matrix of the causal effect estimate learning from the source data.
#' @param LD The linkage disequilibrium (LD) matrix. Default is the identity matrix, assuming independent instrumental variables (IVs).
#' @param cluster.index A vector indicating the LD block indices each IV belongs to. The length is equal to the number of IVs, and values are the LD block indices.
#' @param Rxy A matrix (p+1 x p+1) of the correlation matrix of the p exposures and outcome. The first one should be the transferred linear predictor and last one should be the outcome.
#' @param transfer.coef A scale of transfer.coef of theta.source to theta.target. Default is \code{1}.
#' @param tauvec The candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE. Default is \code{200}.
#' @param pip.thres A scale of PIP theshold for calibyating causality used in SuSiE. Default is \code{0.3}.
#' @param ebic.delta A scale of tuning parameter of causal effect estimate in extended BIC. Default is \code{1}.
#' @param ebic.gamma A scale of tuning parameter of horizontal pleiotropy in extended BIC. Default is \code{2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Default is \code{50}.
#' @param max.eps Tolerance for stopping criteria. Default is \code{1e-4}.
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold. Default is \code{0.8}.
#' @param sampling.time A scale of number of subsampling in estimating the standard error. Default is \code{100}.
#' @param sampling.iter A scale of iteration in subsampling in estimating the standard error. Default is \code{10}.
#' @param gcov A matrix (2 x 2) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.
#' @param prob.shrinkage Exponent for power-law scaling of selection probabilities based on Effective Sample Size (ESS). Controls the balance between favoring information-rich blocks and ensuring diversity. A value of 1 implies probability proportional to ESS; 0 implies uniform probability; 0.5 (default) uses square-root weighting to dampen the dominance of large blocks.

#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom MRBEE MRBEE.IMRP.UV
#' @importFrom MASS rlm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRBEE_TL_UV=function(by,bX,byse,bXse,Rxy,LD="identity",cluster.index=c(1:length(by)),theta.source,theta.source.cov,tauvec=seq(3,30,3),admm.rho=3,ebic.delta=1,ebic.gamma=2,transfer.coef=1,susie.iter=200,pip.thres=0.3,max.iter=50,max.eps=1e-4,reliability.thres=0.8,sampling.time=100,sampling.iter=10,ldsc=NULL,gcov=NULL,prob.shrinkage=0.5){
if(LD[1]=="identity"){
A=MRBEE_TL_UV_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,theta.source=theta.source,theta.source.cov=theta.source.cov,tauvec=tauvec,ebic.delta=ebic.delta,ebic.gamma=ebic.gamma,transfer.coef=transfer.coef,susie.iter=susie.iter,pip.thres=pip.thres,max.iter=max.iter,max.eps=max.eps,reliability.thres=reliability.thres,sampling.time=sampling.time,sampling.iter=sampling.iter,LDSC=ldsc,Omega=gcov)
return(A)
}else{
######### Basic Processing  ##############
fit.no.tran=MRBEE.IMRP.UV(by=by,bx=bX,byse=byse,bxse=bXse,Rxy=Rxy)
theta.source=transfer.coef*theta.source
theta.source.cov=transfer.coef^2*theta.source.cov
theta.ini=fit.no.tran$theta
gamma.ini=fit.no.tran$delta/byse
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
p=1
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
RC=as.matrix(TC%*%LD)
bXinv=as.vector(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.vector(TC%*%bX)
BtB=sum(bXinv*bX)
Thetarho=solve(LD+admm.rho*diag(n))
r=reliability.adj.uv(bX,bXse*sqrt(Rxy[1,1]),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=ldsc,Omega=gcov)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
########## Iteration ###################
Bic=Bic_direct=Btheta=tauvec*0
Bgamma=Bgamma_direct=array(0,c(length(tauvec),n))
theta=theta.ini
for(v in length(tauvec):1){
error=2
iter=0
gamma=gamma.ini
u=gamma1=gamma*0
fit.susie=NULL
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
theta.complement=ifelse(which.min(c(abs(theta),abs(theta-theta.source)))==1,0,theta.source)
by.complement=as.vector(by-bX*theta.complement-LD%*%gamma)
XtX=BtB-sum(bXse[indvalid]^2*Rxy[1,1])
Xty=sum(bXinv*by.complement)-Rxy[1,2]*sum(bXse[indvalid]*byse[indvalid])+sum(bXse[indvalid]^2*theta.complement)
yty=sum(by.complement*(Theta%*%by.complement))
tryCatch({
fit.susie=susie_ss(XtX=as.matrix(XtX),Xty=Xty,yty=yty,L=1,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_ss(XtX=as.matrix(XtX),Xty=Xty,yty=yty,L=1,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F)
})
if(fit.susie$pip>pip.thres){
theta=theta.complement+Xty/XtX
}
if(fit.susie$pip<=pip.thres&theta.source!=0){
theta=theta.source
}
if(fit.susie$pip<=pip.thres&theta.source==0){
theta=0
}
gamma=as.vector(Thetarho%*%(by-bX*theta-u+admm.rho*gamma1))
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
dftheta=as.numeric(fit.susie$pip>pip.thres)
e=as.vector(by-bX*theta-as.vector(LD%*%gamma))
vare=sum(e*(Theta%*%e))/(length(indvalid)-dftheta)
Bic[v]=log(vare)+log(n)/n*dftheta+(1+ebic.gamma)*log(n)*(n-length(indvalid))/n
Btheta[v]=theta
Bgamma[v,]=gamma1

gamma_direct=gamma;gamma_direct1=gamma1;u_direct=u;
for(iter_direct in 1:5){
gamma_direct=as.vector(Thetarho%*%(by-bX*theta.source-u_direct+admm.rho*gamma1_direct))
gamma1_direct=mcp(gamma_direct+u_direct/admm.rho,tauvec[v]/admm.rho)
u_direct=u_direct+admm.rho*(gamma_direct-gamma1_direct)
}
e=as.vector(by-bX*theta.source-as.vector(LD%*%gamma_direct))
vare=sum(e*(Theta%*%e))/sum(gamma1_direct==0)
Bic_direct[v]=log(vare)+(1+ebic.gamma)*log(n)*sum(gamma1_direct!=0)/n
Bgamma_direct[v,]=gamma1_direct
}
s1=min(Bic)
s2=min(Bic_direct)
if(s1<=s2){
############################### final estimate ##########################
vstar=which.min(Bic)
theta=theta.ini
gamma=gamma.ini
gamma1=u=gamma*0
fit.susie=NULL
error=2
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
theta.complement=ifelse(which.min(c(abs(theta),abs(theta-theta.source)))==1,0,theta.source)
by.complement=as.vector(by-bX*theta.complement-LD%*%gamma)
XtX=BtB-sum(bXse[indvalid]^2*Rxy[1,1])
Xty=sum(bXinv*by.complement)-Rxy[1,2]*sum(bXse[indvalid]*byse[indvalid])+sum(bXse[indvalid]^2*theta.complement)
yty=sum(by.complement*(Theta%*%by.complement))
tryCatch({
fit.susie=susie_ss(XtX=as.matrix(XtX),Xty=Xty,yty=yty,L=1,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_ss(XtX=as.matrix(XtX),Xty=Xty,yty=yty,L=1,n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F)
})
if(fit.susie$pip>pip.thres){
theta=theta.complement+Xty/XtX
}
if(fit.susie$pip<=pip.thres&theta.source!=0){
theta=theta.source
}
if(fit.susie$pip<=pip.thres&theta.source==0){
theta=0
}
gamma=as.vector(Thetarho%*%(by-bX*theta-u+admm.rho*gamma1))
gamma1=mcp(gamma+u/admm.rho,tauvec[vstar]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
############################### inference #########################
res=gamma1*byse1
names(res)=rownames(bX)
ThetaList=c(1:sampling.time)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
cluster.index <- as.integer(factor(cluster.index))
cluster_prob <- cluster_prob(cluster.index,LD,alpha=prob.shrinkage)
k <- floor(length(cluster_prob) * 0.5)
while(j<=sampling.time) {
setTxtProgressBar(pb, j)
indicator <- FALSE
tryCatch({
cluster.sampling <- sample(1:max(cluster.index), k, replace = F,prob = cluster_prob)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
nj=length(indj)
bXj=bX[indj]
byj=by[indj]
bXsej=bXse[indj]
bysej=byse[indj]
thetaj=theta*runif(length(theta),0.95,1.05)
RxyListj=RxyList[indj,,]
Rxyallj=biasterm(RxyList=RxyListj,c(1:nj))
gammaj=gamma[indj]*runif(1,0.975,1.025)
uj=gamma1j=gammaj*0
indvalidj=which(gammaj==0)
fit.susiej=fit.susie
deltaj=theta.source-thetaj
LDj=LD[indj,indj]
Thetaj <- Theta[indj,indj]
Thetarhoj <- Thetarho[indj,indj]
bXinvj <- as.vector(Thetaj %*% bXj)
BtBj <- sum(bXinvj*bXj)
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
if(length(indvalidj)==nj){
Rxysumj=Rxyallj
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:nj,indvalidj))
}
theta.complementj=ifelse(which.min(c(abs(thetaj),abs(thetaj-theta.source)))==1,0,theta.source)
by.complementj=as.vector(byj-bXj*theta.complementj-LDj%*%gammaj)
XtXj=BtBj-sum(bXsej[indvalidj]^2*Rxy[1,1])
Xtyj=sum(bXinvj*by.complementj)-Rxy[1,2]*sum(bXsej[indvalidj]*bysej[indvalidj])+sum(bXsej[indvalidj]^2*theta.complementj)
ytyj=sum(by.complementj*(Thetaj%*%by.complementj))
tryCatch({
fit.susiej=susie_ss(XtX=as.matrix(XtXj),Xty=Xtyj,yty=ytyj,L=1,n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=susie.iter,residual_variance_lowerbound=1)
},error = function(e) {
fit.susiej=susie_ss(XtX=as.matrix(XtXj),Xty=Xtyj,yty=ytyj,L=1,n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=susie.iter,estimate_residual_variance=F)
})
if(fit.susiej$pip>pip.thres){
thetaj=theta.complementj+Xtyj/XtXj
}
if(fit.susiej$pip<=pip.thres&theta.source!=0){
thetaj=theta.source
}
if(fit.susiej$pip<=pip.thres&theta.source==0){
thetaj=0
}
gammaj=as.vector(Thetarhoj%*%(byj-bXj*thetaj-uj+admm.rho*gamma1j))
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
}
ThetaList[j]=thetaj
j=j+1
}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, "\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
j <<- j - 1  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
}
close(pb)
theta.cov=var(ThetaList)
if(theta!=0){
theta.cov=theta.cov+theta.source.cov
}
theta.se=sqrt(theta.cov)

A=list()
A$theta=theta
A$gamma=res
A$theta.se=theta.se
A$theta.cov=theta.cov
A$reliability.adjust=r
A$susie.delta=fit.susie
A$theta.list=ThetaList
A$Bic=Bic
A$tau.optimal=tauvec[vstar]
return(A)
}else{
A=list()
A$theta=theta.source
A$theta.se=sqrt(theta.source.cov)
A$gamma=Bgamma_direct[which.min(Bic_direct),]/byse1
cat("BIC shows that theta.source can be directly applied\n")
return(A)
}
}
}
