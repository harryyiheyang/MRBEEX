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
#' @param tauvec The candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(4, 8, by=0.5)}.
#' @param Lvec A vector of the number of single effects used in SuSiE. Default is \code{c(1:6)}.
#' @param standardize If standardize = TRUE, standardize the columns of X to unit variance prior to fitting (or equivalently standardize XtX and Xty to have the same effect) in SuSiE. Note that scaled_prior_variance specifies the prior on the coefficients of X after standardization (if it is performed). If you do not standardize, you may need to think more carefully about specifying scaled_prior_variance. Whatever your choice, the coefficients returned by coef are given for X on the original input scale. Any column of X that has zero variance is not standardized.
#' @param admm.rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param group.penalize An indicator of whether using difference penalty to penalize highly correlated exposures. Defaults to \code{F}.
#' @param group.index A vector of the group index of exposure. Defaults to \code{NULL}.
#' @param group.diff The tuning penalizing difference of highly correlated exposure prediction. Defaults to \code{10}.
#' @param susie.iter A scale of the maximum number of iterations used in SuSiE. Default is \code{200}.
#' @param pip.thres Posterior inclusion probability (PIP) threshold. Individual PIPs less than this value will be shrunk to zero. Default is \code{0.25}.
#' @param pip.min The minimum empirical PIP used in purifying variables in each credible set. Defaults to \code{0.1}.
#' @param coverage.causal The coverage of defining a credible set in MRBEEX when \code{use.susie = T}. Defaults to \code{0.95}.
#' @param estimate_residual_method The method used for estimating residual variance. For the original SuSiE model, "MLE" and "MoM" estimation is equivalent, but for the infinitesimal model, "MoM" is more stable.
#' @param cred.pip.thres The threshold of PIP of each credible set. Defaults to \code{0.95}.
#' @param projection.eigen.floor The minimum eigenvalue used when projecting SuSiE and selected refit cross-product matrices. The full-data floor is this value; resampled matrices are scaled by their current row count divided by the full row count. Defaults to \code{1}.
#' @param ebic.delta A scale of tuning parameter of causal effect estimate in extended BIC. Default is \code{0}.
#' @param ebic.gamma A scale of tuning parameter of horizontal pleiotropy in extended BIC. Default is \code{1}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Default is \code{50}.
#' @param max.eps Tolerance for stopping criteria. Default is \code{1e-4}.
#' @param reliability.thres A scale of threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold. Default is \code{0.6}.
#' @param ridge.diff A scale of parameter on the differences of causal effect estimate in one credible set. Defaults to \code{10}.
#' @param sampling.strategy Resampling scheme used only by the independent branch where \code{LD="identity"}.
#' @param sampling.time Number of fixed half-SNP sampling repeats for standard-error estimation. Default is \code{300}.
#' @param sampling.iter Number of estimation iterations per sampling repeat. Default is \code{25}.
#' @param gcov A matrix (p+1 x p+1) of the per-snp genetic covariance matrix of the p exposures and outcome. The last one should be the outcome.
#' @param ldsc A vector (n x 1) of the LDSCs of the IVs.

#'
#' @return A list containing the estimated causal effect, its covariance, and pleiotropy.
#' @importFrom susieR susie_ss coef.susie susie
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom Matrix Matrix solve chol
#' @importFrom MASS rlm
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
MRBEE_TL=function(by,bX,byse,bXse,Rxy,LD="identity",cluster.index=c(1:length(by)),
  group.penalize=F,group.index=NULL,group.diff=100,
  theta.source,theta.source.cov,tauvec=seq(4,8,0.5),Lvec=c(1:6),standardize=F,
  admm.rho=2,ebic.delta=0,ebic.gamma=1,transfer.coef=1,susie.iter=200,
  pip.thres=0.25,pip.min=0.1,cred.pip.thres=0.95,max.iter=50,coverage.causal=0.95,
  max.eps=1e-6,reliability.thres=0.5,ridge.diff=100,
  estimate_residual_method="MoM",sampling.strategy="bootstrap",
  projection.eigen.floor=1,sampling.time=300,sampling.iter=25,ldsc=NULL,gcov=NULL){
if(LD[1]=="identity"){
A=MRBEE_TL_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,
       theta.source=theta.source,theta.source.cov=theta.source.cov,
       group.penalize=group.penalize,group.index=group.index,group.diff=group.diff,
       tauvec=tauvec,Lvec=Lvec,ebic.delta=ebic.delta,ebic.gamma=ebic.gamma,standardize=standardize,
       transfer.coef=transfer.coef,susie.iter=susie.iter,pip.thres=pip.thres,
       pip.min=pip.min,cred.pip.thres=cred.pip.thres,max.iter=max.iter,max.eps=max.eps,
       estimate_residual_method=estimate_residual_method,sampling.strategy=sampling.strategy,
       reliability.thres=reliability.thres,ridge.diff=ridge.diff,coverage.causal=coverage.causal,
       projection.eigen.floor=projection.eigen.floor,
       sampling.time=sampling.time,sampling.iter=sampling.iter,LDSC=ldsc,Omega=gcov)
return(A)
}else{
######### Basic Processing  ##############
fit.no.tran=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy)
theta.source=transfer.coef*theta.source
theta.source.cov=transfer.coef^2*theta.source.cov
theta.ini=theta.source
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
r=reliability.adj(bX,bXse%*%diag(sqrt(diag(Rxy[1:p,1:p]))),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=ldsc,Omega=gcov)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
br=as.vector(by-bX%*%theta.source)
Diff_matrix=diag(p)*0
if(group.penalize==T){
Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
Veigen=FProject_basis(BtB+Diff_matrix)
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
project_XtX <- new_FProjector(Veigen, eigen.floor=projection.eigen.floor)
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
Cmat=Rxysum[1:p,1:p]
XtX.raw=BtB-Cmat
XtX=project_XtX(XtX.raw, Diff_matrix, indvalid)
Xty=matrixVectorMultiply(Bt,br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum(br.complement*(Theta%*%br.complement))
fit.susie=tryCatch({
susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=0.9,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>pip.min)
delta.latent.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=delta.latent.cs$ind.keep
if(length(pip.alive)>0){
delta.latent[-pip.alive]=0
}else{
delta.latent=delta.latent*0
}
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=project_select_xtx(XtX.raw[inddelta,inddelta,drop=FALSE]+Diff_matrix[inddelta,inddelta,drop=FALSE],eigen.floor=projection.eigen.floor)
xtx=xtx[1,1]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=project_select_xtx(XtX.raw[inddelta,inddelta,drop=FALSE]+Diff_matrix[inddelta,inddelta,drop=FALSE]+ridge.diff*Diff[inddelta,inddelta,drop=FALSE],eigen.floor=projection.eigen.floor)
xty=Xty[inddelta]
delta.latent[inddelta]=c(CppMatrix::matrixSolve(xtx,xty))
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
project_XtX <- new_FProjector(Veigen, eigen.floor=projection.eigen.floor)
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
Cmat=Rxysum[1:p,1:p]
XtX.raw=BtB-Cmat
XtX=project_XtX(XtX.raw, Diff_matrix, indvalid)
Xty=matrixVectorMultiply(Bt,br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum(br.complement*(Theta%*%br.complement))
fit.susie=tryCatch({
susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=0.9,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
delta.latent=coef.susie(fit.susie)[-1]*(fit.susie$pip>pip.min)
delta.latent.cs=group.pip.filter(pip.summary=summary(fit.susie)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=delta.latent.cs$ind.keep
if(length(pip.alive)>0){
delta.latent[-pip.alive]=0
}else{
delta.latent=delta.latent*0
}
inddelta=which(delta.latent!=0)
Diff=generate_block_matrix(summary(fit.susie)$vars,n/diag(BtB),delta.latent)
delta=delta*0
if(length(inddelta)==1){
xtx=project_select_xtx(XtX.raw[inddelta,inddelta,drop=FALSE]+Diff_matrix[inddelta,inddelta,drop=FALSE],eigen.floor=projection.eigen.floor)
xtx=xtx[1,1]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=project_select_xtx(XtX.raw[inddelta,inddelta,drop=FALSE]+Diff_matrix[inddelta,inddelta,drop=FALSE]+ridge.diff*Diff[inddelta,inddelta,drop=FALSE],eigen.floor=projection.eigen.floor)
xty=Xty[inddelta]
delta.latent[inddelta]=c(CppMatrix::matrixSolve(xtx,xty))
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
cat("Resampling process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
resampling.retries=0
while(j<=sampling.time) {
setTxtProgressBar(pb, j)
indicator <- FALSE
tryCatch({
indj <- sort(sample.int(n, size = max(2L, floor(0.5 * n)), replace = FALSE))
nj <- length(indj)
LDj <- Matrix(LD[indj, indj, drop = FALSE], sparse = TRUE)
Thetaj <- solve(LDj)
Cinvj <- LDj
A_gammaj <- Matrix(LD[indj, , drop = FALSE], sparse = TRUE)
bXj <- bX[indj,,drop=FALSE]
byj <- by[indj]
bXsej <- bXse[indj,,drop=FALSE]
bysej <- byse[indj]
brj <- br[indj]
Btj <- as.matrix(t(bXj)%*%Thetaj)
BtBj <- matrixMultiply(Btj,bXj)
BtBj <- (t(BtBj) + BtBj) / 2
dBtBj <- diag(BtBj)
thetaj=theta*runif(length(theta),0.95,1.05)
Rxyallj <- biasterm(RxyList = RxyList, indj)
gammaj=gamma1j=as.vector(gamma*runif(1,0.975,1.025))
uj=gammaj*0
deltaj=theta.source-thetaj
errorj=1
fit.susiej=NULL
projection.eigen.floorj <- projection.eigen.floor*nj/n
project_XtXj <- new_FProjector(Veigen, eigen.floor=projection.eigen.floorj)
for(jiter in 1:sampling.iter){
theta_prevj=thetaj
indvalidj <- which(gamma1j[indj]==0)
if(length(indvalidj)<(0.55*length(indj))){
indvalidj=sample(seq_along(indj), max(1L, floor(0.6*length(indj))))
gamma1j[indj[indvalidj]]=gammaj[indj[indvalidj]]=0
}
invalidj <- which(gamma1j[indj]!=0)
if(length(invalidj)<length(indvalidj)){
Rxysumj <- Rxyallj-biasterm(RxyList = RxyList, indj[invalidj])
}else{
Rxysumj <- biasterm(RxyList = RxyList, indj[indvalidj])
}
fit.clusterj=center.classifying(deltaj,-theta.source)
delta.complementj=fit.clusterj$complement
delta.clusterj=fit.clusterj$clusterj
br.complementj=c(brj-bXj%*%delta.complementj-as.vector(A_gammaj%*%gammaj))
addbiasj=matrixVectorMultiply(Rxysumj[1:p,1:p],theta.source+delta.complementj)
Cmatj=Rxysumj[1:p,1:p]
XtXj.raw=BtBj-Cmatj
XtXj=project_XtXj(XtXj.raw, Diff_matrix, indvalidj)
Xtyj=matrixVectorMultiply(Btj,br.complementj)-Rxysumj[1+p,1:p]+addbiasj
ytyj=sum(br.complementj*(Thetaj%*%br.complementj))
fit.susiej=tryCatch({
susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=ifelse(jiter==1,susie.iter,min(susie.iter,30)),residual_variance_lowerbound=0.9,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=ifelse(jiter==1,susie.iter,min(susie.iter,30)),estimate_residual_variance=F,coverage=coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
delta.latentj=coef.susie(fit.susiej)[-1]*(fit.susiej$pip>pip.min)
delta.latent.csj=group.pip.filter(pip.summary=summary(fit.susiej)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alivej=delta.latent.csj$ind.keep
if(length(pip.alivej)>0){
  delta.latentj[-pip.alivej]=0
}else{
  delta.latentj=delta.latentj*0
}
inddeltaj=which(delta.latentj!=0)
Diffj=generate_block_matrix(summary(fit.susiej)$vars,nj/diag(BtBj),delta.latentj)
deltaj=deltaj*0
if(length(inddeltaj)==1){
xtxj=project_select_xtx(XtXj.raw[inddeltaj,inddeltaj,drop=FALSE]+Diff_matrix[inddeltaj,inddeltaj,drop=FALSE],eigen.floor=projection.eigen.floorj)
xtxj=xtxj[1,1]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=xtyj/xtxj
}
if(length(inddeltaj)>1){
xtxj=project_select_xtx(XtXj.raw[inddeltaj,inddeltaj,drop=FALSE]+Diff_matrix[inddeltaj,inddeltaj,drop=FALSE]+ridge.diff*Diffj[inddeltaj,inddeltaj,drop=FALSE],eigen.floor=projection.eigen.floorj)
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=c(CppMatrix::matrixSolve(xtxj,xtyj))
}
deltaj=delta.latentj+delta.complementj
thetaj=deltaj+theta.source
gamma_centerj <- as.vector(gamma1j - uj/admm.rho)
gamma_residj <- as.vector(byj - matrixVectorMultiply(bXj,thetaj) - A_gammaj%*%gamma_centerj)
gamma_hessianj <- Matrix::forceSymmetric(admm.rho*Cinvj + tcrossprod(A_gammaj))
gamma_middlej <- Matrix::solve(gamma_hessianj, gamma_residj)
gammaj=as.vector(gamma_centerj + crossprod(A_gammaj,gamma_middlej))
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
if(jiter>4) errorj=norm(thetaj-theta_prevj,"2")
if(errorj<max.eps) break
}
ThetaList[j,]=theta.source+deltaj
DeltaList[j,]=deltaj
resampling.retries=0
j=j+1
}, error = function(e) {
resampling.retries <<- resampling.retries+1
if(resampling.retries>20){
stop("MRBEE_TL resampling failed repeatedly: ", conditionMessage(e))
}
indicator <<- TRUE
})
if (indicator) {
next
}
}
close(pb)
theta.cov=covmad(ThetaList)
ind_source=which(delta==0&theta!=0)
if(length(ind_source)>0){
theta.cov[ind_source,ind_source]=theta.cov[ind_source,ind_source]+theta.source.cov[ind_source,ind_source]
}
theta.se=sqrt(diag(theta.cov))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
delta.se=colSDMAD(DeltaList)
delta.cov=covmad(DeltaList)
colnames(delta.cov)=rownames(delta.cov)=names(delta.se)=colnames(bX)

A=list()
A$delta=delta
A$theta=theta.source+delta
A$gamma=res
A$delta.se=delta.se
A$theta.se=theta.se
A$theta.pip=colMeans(ThetaList!=0)
A$delta.pip=colMeans(DeltaList!=0)
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
