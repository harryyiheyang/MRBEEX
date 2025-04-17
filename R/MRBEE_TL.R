#' Estimate Non-Transferable Causal Effect with MRBEE and SuSiE
#'
#' This function estimates the non-transferable causal effect using a bias-correction estimating equation, considering potential pleiotropy and measurement errors, and using SuSiE to select the non-transferable causal effect.
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
#' @param Lvec A vector of the number of single effects used in SuSiE. Default is \code{c(1:6)}.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE. Default is \code{200}.
#' @param pip.thres Posterior inclusion probability (PIP) threshold. Individual PIPs less than this value will be shrunk to zero. Default is \code{0.5}.
#' @param pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.1}.
#' @param cred.pip.thres The threshold of PIP of each credible set. Defaults to \code{0.95}.
#' @param ebic.delta A scale of tuning parameter of causal effect estimate in extended BIC. Default is \code{1}.
#' @param ebic.gamma A scale of tuning parameter of horizontal pleiotropy in extended BIC. Default is \code{2}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Default is \code{50}.
#' @param max.eps Tolerance for stopping criteria. Default is \code{1e-4}.
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold. Default is \code{0.8}.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param sampling.time A scale of number of subsampling in estimating the standard error. Default is \code{100}.
#' @param sampling.iter A scale of iteration in subsampling in estimating the standard error. Default is \code{10}.
#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_suff_stat coef.susie susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom MASS rlm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRBEE_TL=function(by,bX,byse,bXse,Rxy,LD="identity",cluster.index=c(1:length(by)),
                  theta.source,theta.source.cov,tauvec=seq(3,30,3),Lvec=c(1:6),
                  admm.rho=3,ebic.delta=1,ebic.gamma=2,transfer.coef=1,susie.iter=200,
                  pip.thres=0.5, pip.min=0.1,cred.pip.thres=0.95,max.iter=50,
                  max.eps=1e-4,reliability.thres=0.8,ridge.diff=100,
                  sampling.time=100,sampling.iter=10){
if(LD[1]=="identity"){
A=MRBEE_TL_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,
                       theta.source=theta.source,theta.source.cov=theta.source.cov,
                       tauvec=tauvec,Lvec=Lvec,ebic.delta=ebic.delta,ebic.gamma=ebic.gamma,
                       transfer.coef=transfer.coef,susie.iter=susie.iter,pip.thres=0.5,
                       pip.min=0.1,cred.pip.thres=0.95,max.iter=max.iter,max.eps=max.eps,
                       reliability.thres=reliability.thres,ridge.diff=ridge.diff,
                       sampling.time=sampling.time,sampling.iter=sampling.iter)
return(A)
}else{
######### Basic Processing  ##############
fit.no.tran=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy)
theta.source=transfer.coef*theta.source
theta.source.cov=transfer.coef^2*theta.source.cov
theta.ini=fit.no.tran$theta
gamma.ini=fit.no.tran$gamma/byse
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
p=ncol(bX)
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
RC=as.matrix(TC%*%LD)
bXinv=as.matrix(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
Bt=t(bXinv)
BtB=matrixMultiply(Bt,bX)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
Thetarho=solve(LD+admm.rho*diag(n))
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
br=as.vector(by-bX%*%theta.source)
########## Iteration ###################
Bic=matrix(0,length(Lvec),length(tauvec))
Btheta=array(0,c(length(Lvec),length(tauvec),p))
Bgamma=array(0,c(length(Lvec),length(tauvec),n))
for(i in 1:length(Lvec)){
fit.susie=NULL
delta=theta.source-theta.ini
for(v in length(tauvec):1){
error=2
iter=0
gamma=gamma.ini
u=gamma1=gamma*0
while(error>max.eps&iter<max.iter){
delta1=delta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
fit.cluster=center.classifying(delta,-theta.source)
delta.complement=fit.cluster$complement
delta.cluster=fit.cluster$cluster
br.complement=as.vector(br-bX%*%delta.complement-as.vector(LD%*%gamma))
addbias=matrixVectorMultiply(Rxysum[1:p,1:p],theta.source+delta.complement)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2
Xty=matrixVectorMultiply(Bt,br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum(br.complement*(Theta%*%br.complement))
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>pip.min)
delta.latent.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=delta.latent.cs$ind.keep
delta.latent[-pip.alive]=0
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1))
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
e=as.vector(br-bX%*%theta-as.vector(LD%*%gamma))
vare=sum(e*(Theta%*%e))/(length(indvalid)-length(inddelta))
Bic[i,v]=log(vare)+(log(n)+log(p)*ebic.delta)/n*length(inddelta)+(1+ebic.gamma)*log(n)*(n-length(indvalid))/n
Btheta[i,v,]=theta
Bgamma[i,v,]=gamma
}
}
############################### final estimate ##########################
istar=bimin(Bic)[1]
vstar=bimin(Bic)[2]
theta=Btheta[istar,vstar,]
gamma=Bgamma[istar,vstar,]
u=gamma1=gamma*0
delta=theta.source-theta
fit.susie=NULL
error=2
iter=0
while(error>max.eps&iter<max.iter){
delta1=delta
indvalid=which(gamma1==0)
if(length(indvalid)==n){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:n,indvalid))
}
fit.cluster=center.classifying(delta,-theta.source)
delta.complement=fit.cluster$complement
delta.cluster=fit.cluster$cluster
br.complement=as.vector(br-bX%*%delta.complement-as.vector(LD%*%gamma))
addbias=matrixVectorMultiply(Rxysum[1:p,1:p],theta.source+delta.complement)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2
Xty=matrixVectorMultiply(Bt,br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum(br.complement*(Theta%*%br.complement))
tryCatch({
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susie=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,s_init=fit.susie,standardize=F,max_iter=susie.iter,intercept=F,estimate_residual_variance=F)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>pip.min)
delta.latent.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=delta.latent.cs$ind.keep
delta.latent[-pip.alive]=0
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1))
gamma1=mcp(gamma+u/admm.rho,tauvec[vstar]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
############################### inference #########################
names(delta)=colnames(bX)
theta=theta.source+delta
res=gamma1*byse1
names(res)=rownames(bX)
ThetaList=DeltaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(DeltaList)=colnames(bX)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time) {
setTxtProgressBar(pb, j)
indicator <- FALSE
tryCatch({
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
nj=length(indj)
bXj=bX[indj,]
byj=by[indj]
bXsej=bXse[indj,]
bysej=byse[indj]
brj=br[indj]
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
Btj <- as.matrix(t(bX[indj, ]) %*% Thetaj)
BtBj <- Btj%*%bX[indj, ]
BtBj=(t(BtBj)+BtBj)/2
dBtBj=diag(BtBj)
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
if(length(indvalidj)==nj){
Rxysumj=Rxyallj
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:nj,indvalidj))
}
fit.clusterj=center.classifying(deltaj,-theta.source)
delta.complementj=fit.clusterj$complement
delta.clusterj=fit.clusterj$clusterj
br.complementj=c(brj-bXj%*%delta.complementj-gammaj)
addbiasj=matrixVectorMultiply(Rxysumj[1:p,1:p],theta.source+delta.complementj)
XtXj=BtBj-Rxysumj[1:p,1:p]
XtXj=t(XtXj)/2+XtXj/2
Xtyj=matrixVectorMultiply(Btj,br.complementj)-Rxysumj[1+p,1:p]+addbiasj
ytyj=sum(br.complementj*(Thetaj%*%br.complementj))
tryCatch({
fit.susiej=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=15,intercept=F,residual_variance_lowerbound=1)
},error = function(e) {
fit.susiej=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,s_init=fit.susiej,standardize=F,max_iter=15,intercept=F,estimate_residual_variance=F)
})
delta.latentj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>pip.min)
delta.latent.csj=group.pip.filter(pip.summary=summary(fit.susiej)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alivej=delta.latent.csj$ind.keep
delta.latentj[-pip.alivej]=0
inddeltaj=which(delta.latentj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(BtBj),delta.latentj)
deltaj=deltaj*0
if(length(inddeltaj)==1){
xtxj=XtXj[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=xtyj/xtxj
}
if(length(inddeltaj)>1){
xtxj=XtXj[inddeltaj,inddeltaj]+ridge.diff*Diffj[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=c(solve(xtxj)%*%xtyj)
}
deltaj=delta.latentj+delta.complementj
thetaj=deltaj+theta.source
gammaj=as.vector(Thetarhoj%*%(byj-matrixVectorMultiply(bXj,thetaj)-uj+admm.rho*gamma1j))
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
}
ThetaList[j,]=theta.source+deltaj
DeltaList[j,]=deltaj
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
theta.cov=cov(ThetaList)*n/length(indvalid)
theta.cov[which(delta==0),which(delta==0)]=theta.cov[which(delta==0),which(delta==0)]+theta.source.cov[which(delta==0),which(delta==0)]
theta.se=sqrt(diag(theta.cov))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
delta.se=colSD(DeltaList)*sqrt(n/length(indvalid))
delta.cov=cov(DeltaList)*n/length(indvalid)
colnames(delta.cov)=rownames(delta.cov)=names(delta.se)=colnames(bX)

A=list()
A$delta=delta
A$theta=theta.source+delta
A$gamma=res
A$delta.se=delta.se
A$theta.se=theta.se
A$delta.cov=delta.cov
A$theta.cov=theta.cov
A$reliability.adjust=r
A$susie.delta=fit.susie
A$delta.latent=delta.latent
A$delta.list=DeltaList
A$theta.list=ThetaList
A$Bic=Bic
A$L.optimal=Lvec[istar]
A$tau.optimal=tauvec[vstar]
return(A)
}
}
