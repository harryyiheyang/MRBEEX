#' @title MRBEE_LDA with Data Fission Inference
#' @description
#' A resampling-free inference method for LD-aware multivariable MR. The data is
#' split into two independent views via Gaussian data fission. The selection view
#' determines both theta and gamma supports. The inference view performs an
#' unpenalized fixed-support refit with a joint sandwich standard error.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param cluster.index A vector indicating the LD block indices each IV belongs to.
#' @param epsilon The thinning fraction allocated to fold 1 (selection). Must be between \code{0.1} and \code{0.9}, inclusive. Default is \code{0.5}.
#' @param reliability.thres Minimum reliability ratio threshold. Default is \code{0.6}.
#' @param estimate_residual_method Method for estimating residual variance in SuSiE. Default is \code{"MoM"}.
#' @param group.penalize Whether to penalize differences of correlated exposures. Default is \code{FALSE}.
#' @param group.index Group index of exposures. Default is \code{c(1:ncol(bX))}.
#' @param group.diff Penalty on correlated exposure differences. Default is \code{100}.
#' @param tauvec Candidate tuning parameters for MCP penalty. Default is \code{seq(4,8,by=0.5)}.
#' @param admm.rho ADMM tuning parameter. Default is \code{2}.
#' @param Lvec Candidate number of single effects in SuSiE. Default is \code{c(1:min(10, ncol(bX)))}.
#' @param pip.thres PIP threshold. Default is \code{0.25}.
#' @param pip.min Minimum PIP for credible set purification. Default is \code{0.1}.
#' @param estimate_residual_variance Whether to estimate residual variance in SuSiE. Default is \code{TRUE}.
#' @param cred.pip.thres Credible set PIP threshold. Default is \code{0.95}.
#' @param coverage.causal Coverage for credible sets. Default is \code{0.95}.
#' @param standardize Whether to standardize in SuSiE. Default is \code{FALSE}.
#' @param projection.eigen.floor Minimum eigenvalue for projections. Default is \code{1}.
#' @param max.iter Maximum iterations. Default is \code{50}.
#' @param max.eps Convergence tolerance. Default is \code{1e-5}.
#' @param susie.iter SuSiE iterations per IPOD iteration. Default is \code{100}.
#' @param maxdiff Maximum norm ratio for theta. Default is \code{3}.
#' @param ridge.diff Ridge on credible set differences. Default is \code{1e3}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{0}.
#' @param ebic.gamma EBIC factor on pleiotropy. Default is \code{1}.
#' @param theta.ini Initial theta. For valid data-fission inference, supplied
#'   values must be independent of the input summary data. Default is \code{FALSE}.
#' @param gamma.ini Initial gamma. For valid data-fission inference, supplied
#'   values must be independent of the input summary data. Default is \code{FALSE}.
#' @param verbose Whether to print progress. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{theta}}{Causal effect estimate.}
#'   \item{\code{theta.se}}{Standard error of the causal effect estimate.}
#'   \item{\code{theta.cov}}{Covariance matrix of the causal effect estimate.}
#'   \item{\code{gamma}}{Estimate of horizontal pleiotropy (on original scale).}
#'   \item{\code{Bic}}{BIC values from fold 1 selection.}
#'   \item{\code{reliability.adjust}}{The full-data global vector of covariance scaling factors.}
#'   \item{\code{epsilon}}{Thinning fraction used.}
#' }
#'
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixListProduct
#' @importFrom Matrix Matrix solve chol
#' @importFrom susieR susie_ss coef.susie
#' @importFrom stats rnorm
#' @noRd

MRBEE_LDA_Datafission=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),
                            epsilon=0.5,standardize=F,
                            group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=100,
                            tauvec=seq(4,8,by=0.5),admm.rho=2,
                            Lvec=c(1:min(10,ncol(bX))),pip.thres=0.25,
                            estimate_residual_variance=T,estimate_residual_method="MoM",
                            pip.min=0.1,cred.pip.thres=0.95,
                            max.iter=50,max.eps=1e-5,susie.iter=100,
                            ebic.theta=0,ebic.gamma=1,ridge.diff=1e3,
                            maxdiff=3,reliability.thres=0.6,coverage.causal=0.95,
                            theta.ini=F,gamma.ini=F,verbose=T,
                            projection.eigen.floor=1){

t1=Sys.time()
m=length(by)
p=ncol(bX)
if(length(epsilon)!=1||!is.numeric(epsilon)||!is.finite(epsilon)||epsilon<0.1||epsilon>0.9){
stop("epsilon must be between 0.1 and 0.9, inclusive.")
}
if(length(cluster.index)!=m||any(is.na(cluster.index))){
stop("cluster.index must have one non-missing value for each IV.")
}
cluster.index=as.integer(factor(cluster.index))
isLD=!(is.character(LD)&&length(LD)==1&&identical(LD,"identity"))

########################################################################
## Step 0: Data Fission
########################################################################
if(verbose) cat("Data fission (epsilon =",epsilon,")...\n")

L_Rxy=t(chol(Rxy))
if(isLD){
if(inherits(LD,"sparseMatrix")){
LD_sp=Matrix::forceSymmetric(LD)
RC_LD=chol(LD_sp)
}else{
LD_dense=as.matrix(LD)
RC_LD=chol(LD_dense)
LD_sp=Matrix(LD_dense,sparse=T)
}
Theta=solve(LD_sp)
}else{
LD_sp=Matrix::Diagonal(m)
Theta=LD_sp
}

byseinv=1/byse
bX0=bX*byseinv
bXse0=bXse*byseinv
r=reliability.adj(bX0,bXse0%*%diag(sqrt(diag(Rxy[1:p,1:p]))),
                  Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy_adj=t(t(Rxy)*r)*r

W=matrix(rnorm(m*(p+1)),m,p+1)
if(isLD){
Z_raw=as.matrix(Matrix::t(RC_LD)%*%W%*%t(L_Rxy))
}else{
Z_raw=W%*%t(L_Rxy)
}

se_matrix=cbind(bXse,byse)
Z=Z_raw*se_matrix
noise=sqrt(epsilon*(1-epsilon))*Z

bX1=(epsilon*bX+noise[,1:p,drop=FALSE])/epsilon
by1=(epsilon*by+noise[,p+1])/epsilon
bXse1=bXse/sqrt(epsilon)
byse1=byse/sqrt(epsilon)

bX2=((1-epsilon)*bX-noise[,1:p,drop=FALSE])/(1-epsilon)
by2=((1-epsilon)*by-noise[,p+1])/(1-epsilon)
bXse2=bXse/sqrt(1-epsilon)
byse2=byse/sqrt(1-epsilon)

t2=Sys.time()
if(verbose) cat(paste0("Fission done: ",round(difftime(t2,t1,units="secs"),3)," secs\n"))

########################################################################
## Step 1: Model Selection on Fold 1
########################################################################
t1=Sys.time()
if(verbose) cat("Model selection on fold 1...\n")

by_s=by1/byse1
byseinv1=1/byse1
bX_s=bX1*byseinv1
bXse_s=bXse1*byseinv1
byse_s=byse1/byse1

if(isLD){
bXinv=as.matrix(Theta%*%bX_s)
Bt=t(bXinv)
BtB=matrixMultiply(Bt,bX_s)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
Thetarho=solve(LD_sp+Matrix::Diagonal(m,admm.rho))
}else{
bXinv=bX_s
Bt=t(bX_s)
BtB=matrixMultiply(bX_s,bX_s,transA=TRUE)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
Thetarho=Matrix::Diagonal(m,1/(1+admm.rho))
}

RxyList=IVweight(byse_s,bXse_s,Rxy_adj)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
Diff_matrix=diag(p)*0
if(group.penalize==T){
Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
Veigen=FProject_basis(BtB+Diff_matrix)

theta.ini.missing=length(theta.ini)==1&&((is.logical(theta.ini)&&!theta.ini)||
                  (is.numeric(theta.ini)&&is.finite(theta.ini)&&theta.ini==0))
if(theta.ini.missing){
fit0=tryCatch({
susie_ss(XtX=BtB+Diff_matrix,Xty=c(matrixMultiply(bXinv,by_s,transA=TRUE)),
         yty=sum(by_s*(Theta%*%by_s)),L=10,n=m,coverage=coverage.causal,standardize=standardize)
},error=function(e){
susie_ss(XtX=BtB+Diff_matrix,Xty=c(matrixMultiply(bXinv,by_s,transA=TRUE)),
         yty=sum(by_s*(Theta%*%by_s)),L=10,n=m,coverage=coverage.causal,
         standardize=standardize,estimate_residual_variance=F)
})
theta.ini=coef.susie(fit0)[-1]*(fit0$pip>0.5)
if(!any(theta.ini!=0)){
XtX.ini.raw=BtB-Rxyall[1:p,1:p]
project_XtX.ini=new_FProjector(Veigen,eigen.floor=projection.eigen.floor)
XtX.ini=project_XtX.ini(XtX.ini.raw,Diff_matrix,c(1:m))
Xty.ini=matrixVectorMultiply(Bt,by_s)-Rxyall[1:p,p+1]
theta.ini=theta_mcp_initialization(
XtX=XtX.ini,Xty=Xty.ini,yty=sum(by_s*(Theta%*%by_s)),n=m,
dfmax=min(max(Lvec),p)
)
}
gamma.ini=by_s*0
}else{
if(length(theta.ini)!=p){
stop("theta.ini must have one value for each exposure.")
}
gamma.ini=gamma.ini/byse1
}

vary=Rxy[p+1,p+1]

## BIC grid search over (L, tau) with SuSiE
Btheta=array(0,c(p,length(tauvec),length(Lvec)))
Bgamma=array(0,c(m,length(tauvec),length(Lvec)))
Bbic=matrix(0,length(tauvec),length(Lvec))

for(v in length(Lvec):1){
fit.theta=NULL
for(j in length(tauvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
project_XtX=new_FProjector(Veigen,eigen.floor=projection.eigen.floor)
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)<(0.55*m)){
indvalid=sample(m,0.6*m)
gamma1[indvalid]=gamma[indvalid]=0
}
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by_s-as.vector(LD_sp%*%gamma)
Cmat=Rxysum[1:p,1:p]
XtX.raw=BtB-Cmat
XtX=project_XtX(XtX.raw,Diff_matrix,indvalid)
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=tryCatch({
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],
         estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),
         max_iter=susie.iter,model_init=fit.theta,coverage=coverage.causal,
         estimate_residual_variance=estimate_residual_variance,
         residual_variance=max(0.9,vary),
         estimate_residual_method=estimate_residual_method,standardize=standardize)
},error=function(e){
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],
         estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),
         estimate_residual_variance=F,residual_variance=max(0.9,vary),
         max_iter=susie.iter,model_init=fit.theta,coverage=coverage.causal,
         estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.min)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,
         xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
if(length(pip.alive)>0){
theta[-pip.alive]=0
}else{
theta=theta*0
}
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.theta)$vars,m/dBtB,theta)
if(length(indtheta)==1){
xtx=project_select_xtx(XtX.raw[indtheta,indtheta,drop=FALSE]+Diff_matrix[indtheta,indtheta,drop=FALSE],eigen.floor=projection.eigen.floor)
xtx=xtx[1,1]
xty=Xty[indtheta]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX_sel=project_select_xtx(XtX.raw[indtheta,indtheta,drop=FALSE]+Diff_matrix[indtheta,indtheta,drop=FALSE]+ridge.diff*Diff[indtheta,indtheta,drop=FALSE],eigen.floor=projection.eigen.floor)
xty=Xty[indtheta]
theta[indtheta]=c(CppMatrix::matrixSolve(XtX_sel,xty))
}
if((norm(theta,"2")/max(norm(theta.ini,"2"),1e-6))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}
if(isLD){
gamma=as.vector(Thetarho%*%(by_s-matrixVectorMultiply(bX_s,theta)-delta+admm.rho*gamma1))
gamma1=mcp(gamma+delta/admm.rho,tauvec[j]/admm.rho)
delta=delta+admm.rho*(gamma-gamma1)
}else{
gamma=(by_s-matrixVectorMultiply(bX_s,theta)-delta+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+delta/admm.rho,tauvec[j]/admm.rho)
delta=delta+admm.rho*(gamma-gamma1)
}
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta[,j,v]=theta
Bgamma[,j,v]=gamma
df1=sum(gamma1!=0)
df2=min(Lvec[v],sum(theta!=0))
res=c(by_s-matrixVectorMultiply(bX_s,theta)-as.vector(LD_sp%*%gamma))
rss=sum(res*(Theta%*%res))/(m-df1-df2)
Bbic[j,v]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+df2*(log(m)+ebic.theta*log(p))
}
}
Bbic=Bbic/m
jstar=bimin(Bbic)[1]
vstar=bimin(Bbic)[2]
theta_sel=Btheta[,jstar,vstar]
gamma_sel=Bgamma[,jstar,vstar]
Bic_out=Bbic

indtheta=which(theta_sel!=0)
indgamma_sel=which(gamma_sel!=0)

t2=Sys.time()
if(verbose) cat(paste0("Selection done: ",round(difftime(t2,t1,units="secs"),3)," secs\n"))

########################################################################
## Step 2: Fixed-Support Inference on Fold 2
########################################################################
t1=Sys.time()
if(verbose) cat("Fixed-support inference on fold 2...\n")

by2_s=by2/byse2
byseinv2=1/byse2
bX2_s=bX2*byseinv2
bXse2_s=bXse2*byseinv2
byse2_s=byse2/byse2

if(isLD){
bXinv2=as.matrix(Theta%*%bX2_s)
Bt2=t(bXinv2)
BtB2=matrixMultiply(Bt2,bX2_s)
BtB2=t(BtB2)/2+BtB2/2
}else{
bXinv2=bX2_s
Bt2=t(bX2_s)
BtB2=matrixMultiply(bX2_s,bX2_s,transA=TRUE)
BtB2=t(BtB2)/2+BtB2/2
}

RxyList2=IVweight(byse2_s,bXse2_s,Rxy_adj)
Rxyall2=biasterm(RxyList2,c(1:m))
indgamma=indgamma_sel
indvalid2=setdiff(seq_len(m),indgamma)
if(length(indvalid2)==m){
Rxysum2=Rxyall2
}else{
Rxysum2=Rxyall2-biasterm(RxyList2,setdiff(1:m,indvalid2))
}

theta=rep(0,p)
gamma=rep(0,m)
if(length(indtheta)>0){
nt=length(indtheta)
ng=length(indgamma)
Htt=BtB2[indtheta,indtheta,drop=FALSE]-
    Rxysum2[indtheta,indtheta,drop=FALSE]+Diff_matrix[indtheta,indtheta,drop=FALSE]
gtheta=matrixVectorMultiply(Bt2,by2_s)-Rxysum2[1:p,p+1]
if(length(indgamma)>0){
Htg=t(bX2_s[indgamma,indtheta,drop=FALSE])
Hgg=as.matrix(LD_sp[indgamma,indgamma,drop=FALSE])
H=rbind(cbind(Htt,Htg),cbind(t(Htg),Hgg))
g=c(gtheta[indtheta],by2_s[indgamma])
}else{
H=Htt
g=gtheta[indtheta]
}
H=(H+t(H))/2
if(rcond(H)<sqrt(.Machine$double.eps)){
stop("The fixed-support data-fission estimating equation is numerically singular.")
}
eta=c(CppMatrix::matrixSolve(H,g))
theta[indtheta]=eta[seq_len(nt)]
if(ng>0) gamma[indgamma]=eta[nt+seq_len(ng)]

res=as.vector(by2_s-matrixVectorMultiply(bX2_s,theta)-as.vector(LD_sp%*%gamma))
npar=nt+ng
E=matrix(0,m,npar)
E[,seq_len(nt)]=-bXinv2[,indtheta,drop=FALSE]*as.vector(res)
if(ng>0){
E[cbind(indgamma,nt+seq_len(ng))]=-res[indgamma]
}

for(i in 1:length(indvalid2)){
ii=indvalid2[i]
Rxy_ii=RxyList2[ii,,]
E[ii,seq_len(nt)]=E[ii,seq_len(nt)]+Rxy_ii[indtheta,p+1]-
    as.vector(Rxy_ii[indtheta,indtheta,drop=FALSE]%*%theta[indtheta])
}
V=crossprod(E)
dfden=length(indvalid2)-nt
if(dfden<=0){
stop("Too few valid IVs for fixed-support data-fission inference.")
}
adjf=m/dfden
Hinv=CppMatrix::matrixSolve(H,diag(npar))
covtheta_full=(Hinv%*%V%*%Hinv)*adjf

theta.cov=diag(p)*0
theta.cov[indtheta,indtheta]=covtheta_full[seq_len(nt),seq_len(nt),drop=FALSE]
theta.se=sqrt(pmax(diag(theta.cov),0))
}else{
if(length(indgamma)>0){
Hgg=as.matrix(LD_sp[indgamma,indgamma,drop=FALSE])
Hgg=(Hgg+t(Hgg))/2
if(rcond(Hgg)<sqrt(.Machine$double.eps)){
stop("The selected gamma design is numerically singular on the inference view.")
}
gamma[indgamma]=c(CppMatrix::matrixSolve(Hgg,by2_s[indgamma]))
}
theta.cov=diag(p)*0
theta.se=rep(0,p)
}

names(theta)=colnames(bX)
names(gamma)=rownames(bX)

colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)

t2=Sys.time()
if(verbose) cat(paste0("Inference done: ",round(difftime(t2,t1,units="secs"),3)," secs\n"))

########################################################################
## Output
########################################################################
A=list()
A$theta=theta
A$theta.se=theta.se
A$theta.cov=theta.cov
A$gamma=gamma*byse2
A$Bic=Bic_out
A$reliability.adjust=r
A$epsilon=epsilon
return(A)
}
