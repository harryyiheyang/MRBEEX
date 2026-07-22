MRBEE_IPOD_SuSiE=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),Lvec=c(1:min(10,nrow(bX))),pip.thres=0.5,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=1,ebic.gamma=2,reliability.thres=0.5,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=10,theta.ini=F,gamma.ini=F,ridge.diff=1e5,projection.eigen.floor=1,verbose=T,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,coverage.causal=0.95,estimate_residual_variance=T,estimate_residual_method="MoM",standardize=F,group.size=4){
########################### Basic information #######################
t1=Sys.time()
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
vary=Rxy[p+1,p+1]
if(LD[1]!="identity"){
isLD=T
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
bXinv=as.matrix(Theta%*%bX)
Bt=t(bXinv)
BtB=matrixMultiply(Bt,bX)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
Thetarho=solve(LD+rho*diag(m))
}else{
isLD=F
LD=Theta=TC=Matrix(diag(m),sparse=T)
bXinv=tilde.X=bX
Bt=t(bX)
BtB=matrixMultiply(bX,bX,transA=TRUE)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
tilde.y=by
Thetarho=diag(m)*1/(1+rho)
Thetarho=Matrix(Thetarho,sparse=T)
}
r=reliability.adj(bX,bXse%*%diag(sqrt(diag(Rxy[1:p,1:p]))),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
Diff_matrix=diag(p)*0
if(group.penalize==T){
Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
Veigen=FProject_basis(BtB+Diff_matrix)
############################ Initial Estimate #######################
theta.ini.missing=length(theta.ini)==1&&((is.logical(theta.ini)&&!theta.ini)||
                  (is.numeric(theta.ini)&&is.finite(theta.ini)&&theta.ini==0))
if(theta.ini.missing){
fit0=susie_ss(XtX=BtB+Diff_matrix,Xty=c(matrixMultiply(bXinv,by,transA=TRUE)),yty=sum(by*(Theta%*%by)),L=10,n=m,standardize=standardize)
theta.ini=coef(fit0)[-1]*(fit0$pip>0.5)
if(!any(theta.ini!=0)){
XtX.ini.raw=BtB-Rxyall[1:p,1:p]
project_XtX.ini=new_FProjector(Veigen,eigen.floor=projection.eigen.floor)
XtX.ini=project_XtX.ini(XtX.ini.raw,Diff_matrix,c(1:m))
Xty.ini=matrixVectorMultiply(Bt,by)-Rxyall[1:p,p+1]
theta.ini=theta_mcp_initialization(
XtX=XtX.ini,Xty=Xty.ini,yty=sum(by*(Theta%*%by)),n=m,
dfmax=min(max(Lvec),p)
)
}
theta.ini1=theta.ini
gamma.ini=gamma.ini1=by*0
}else{
if(length(theta.ini)!=p){
stop("theta.ini must have one value for each exposure.")
}
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
}
theta.ini.norm=norm(theta.ini1,"2")
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Initialization ends: ",time_to_print," secs\n"))
}
############################## Tuning Parameter ######################
t1=Sys.time()
w=length(tauvec)
q=length(Lvec)
Btheta=array(0,c(p,w,q))
Bgamma=array(0,c(m,w,q))
Bbic=matrix(0,w,q)
for(v in length(Lvec):1){
fit.theta=NULL
for(j in length(tauvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
project_XtX <- new_FProjector(Veigen, eigen.floor=projection.eigen.floor)

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
res.theta=by-as.vector(LD%*%gamma)
Cmat=Rxysum[1:p,1:p]
XtX.raw=BtB-Cmat
XtX=project_XtX(XtX.raw, Diff_matrix, indvalid)
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=tryCatch({
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta=coef.susie(fit.theta)[-1]
theta.summary=summary(fit.theta)$vars
pip.alive=theta.summary$variable[theta.summary$cs>0&theta.summary$variable_prob>=pip.thres]
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
XtX=project_select_xtx(XtX.raw[indtheta,indtheta,drop=FALSE]+Diff_matrix[indtheta,indtheta,drop=FALSE]+ridge.diff*Diff[indtheta,indtheta,drop=FALSE],eigen.floor=projection.eigen.floor)
Xty=Xty[indtheta]
theta[indtheta]=c(CppMatrix::matrixSolve(XtX,Xty))
}
theta.norm=norm(theta,"2")
if(!is.finite(theta.norm)){
stop("Theta estimation produced non-finite values.")
}
if(theta.ini.norm>0&&(theta.norm/theta.ini.norm)>maxdiff){
theta=theta/theta.norm*maxdiff*theta.ini.norm
}
if(isLD){
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
}else{
gamma=(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)/(1+rho)
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
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
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
rss=sum(res*(Theta%*%res))/(m-df1-df2)
Bbic[j,v]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+df2*(log(m)+ebic.theta*log(p))
}
}
Bbic=Bbic/m
jstar=bimin(Bbic)[1]
vstar=bimin(Bbic)[2]
theta.ini=Btheta[,jstar,vstar]
gamma.ini=Bgamma[,jstar,vstar]
indvalid=which(gamma.ini==0)
theta=theta.ini
gamma=gamma1=gamma.ini
delta=0*gamma
error=1
iter=1
fit.theta=NULL
project_XtX <- new_FProjector(Veigen, eigen.floor=projection.eigen.floor)
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
res.theta=by-as.vector(LD%*%gamma)
Cmat=Rxysum[1:p,1:p]
XtX.raw=BtB-Cmat
XtX=project_XtX(XtX.raw, Diff_matrix, indvalid)
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=tryCatch({
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],estimate_prior_method=ifelse(is.null(fit.theta),"optim","EM"),estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta=coef.susie(fit.theta)[-1]
theta.summary=summary(fit.theta)$vars
pip.alive=theta.summary$variable[theta.summary$cs>0&theta.summary$variable_prob>=pip.thres]
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
XtX=project_select_xtx(XtX.raw[indtheta,indtheta,drop=FALSE]+Diff_matrix[indtheta,indtheta,drop=FALSE]+ridge.diff*Diff[indtheta,indtheta,drop=FALSE],eigen.floor=projection.eigen.floor)
Xty=Xty[indtheta]
theta[indtheta]=c(CppMatrix::matrixSolve(XtX,Xty))
}
theta.norm=norm(theta,"2")
if(!is.finite(theta.norm)){
stop("Theta estimation produced non-finite values.")
}
if(theta.ini.norm>0&&(theta.norm/theta.ini.norm)>maxdiff){
theta=theta/theta.norm*maxdiff*theta.ini.norm
}
if(isLD){
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
}else{
gamma=(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)/(1+rho)
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
}
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Estimation ends: ",time_to_print," secs\n"))
}
############################### inference #########################
t1=Sys.time()
gamma=gamma
theta=theta
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma1!=0)
indvalid=which(gamma1==0)
res=by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma1)
ThetaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(bX)
GammaList=matrix(0,sampling.time,m)
colnames(GammaList)=rownames(bX)
if(verbose){
cat("Resampling starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
}
j=1
resampling.retries=0
while(j<=sampling.time){
indicator <- FALSE
resampling.stage <- "start"
if(verbose) setTxtProgressBar(pb,j)
tryCatch({
indj <- group_subsampling_indices(cluster.index, group.size = group.size, min.size = 2L)
mj <- length(indj)
LDj <- Matrix(LD[indj, indj, drop = FALSE], sparse = TRUE)
resampling.stage <- "LD solve"
Thetaj <- solve(LDj)
Cinvj <- LDj
A_gammaj <- Matrix(LD[indj, , drop = FALSE], sparse = TRUE)
bXj <- bX[indj,,drop=FALSE]
bXsej <- bXse[indj,,drop=FALSE]
byj <- by[indj]
bysej <- byse[indj]
Btj <- as.matrix(t(bXj) %*% Thetaj)
BtBj <- matrixMultiply(Btj,bXj)
BtBj <- (t(BtBj) + BtBj) / 2
dBtBj <- diag(BtBj)
Rxyallj <- biasterm(RxyList = RxyList, indj)
thetaj=theta*runif(p,0.95,1.05)
gammaj=gamma1j=as.vector(gamma*runif(1,0.975,1.025))
deltaj=gammaj*0
errorj=1
fit.thetaj=NULL
projection.eigen.floorj <- projection.eigen.floor*mj/m
project_XtXj <- new_FProjector(Veigen, eigen.floor=projection.eigen.floorj)

for(jiter in 1:sampling.iter){
theta_prevj=thetaj
resampling.stage <- "bias correction"
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
res.thetaj=byj-as.vector(A_gammaj%*%gammaj)
Cmatj=Rxysumj[1:p,1:p]
XtXj.raw=BtBj-Cmatj
XtXj=project_XtXj(XtXj.raw, Diff_matrix/2, indvalidj)
Xtyj=matrixVectorMultiply(Btj,res.thetaj)-Rxysumj[1:p,p+1]
ytyj=sum(res.thetaj*(Thetaj%*%res.thetaj))
resampling.stage <- "susie"
fit.thetaj=tryCatch({
susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,n=mj,L=Lvec[vstar],estimate_prior_method="EM",max_iter=ifelse(jiter==1,susie.iter,min(susie.iter,30)),model_init=fit.thetaj,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,n=mj,L=Lvec[vstar],estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=ifelse(jiter==1,susie.iter,min(susie.iter,30)),model_init=fit.thetaj,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
thetaj=coef.susie(fit.thetaj)[-1]
theta.summaryj=summary(fit.thetaj)$vars
pip.thresj=max(pip.thres/sqrt(2),0.1)
pip.alivej=theta.summaryj$variable[theta.summaryj$cs>0&theta.summaryj$variable_prob>=pip.thresj]
if(length(pip.alivej)>0){
  thetaj[-pip.alivej]=0
}else{
  thetaj=thetaj*0
}
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.thetaj)$vars,m/dBtBj,thetaj)
if(length(indthetaj)==1){
resampling.stage <- "single-effect refit"
xtxj=project_select_xtx(XtXj.raw[indthetaj,indthetaj,drop=FALSE]+Diff_matrix[indthetaj,indthetaj,drop=FALSE]/2,eigen.floor=projection.eigen.floorj)
xtxj=xtxj[1,1]
xtyj=Xtyj[indthetaj]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
resampling.stage <- "multi-effect refit"
XtXj=project_select_xtx(XtXj.raw[indthetaj,indthetaj,drop=FALSE]+Diff_matrix[indthetaj,indthetaj,drop=FALSE]/2+ridge.diff*Diffj[indthetaj,indthetaj,drop=FALSE],eigen.floor=projection.eigen.floorj)
Xtyj=Xtyj[indthetaj]
thetaj[indthetaj]=c(CppMatrix::matrixSolve(XtXj,Xtyj))
}
thetaj.norm=norm(thetaj,"2")
if(!is.finite(thetaj.norm)){
stop("Resampling produced non-finite theta values.")
}
if(theta.ini.norm>0&&(thetaj.norm/theta.ini.norm)>maxdiff){
thetaj=thetaj/thetaj.norm*maxdiff*theta.ini.norm
}
resampling.stage <- "gamma projection"
gamma_centerj <- as.vector(gamma1j - deltaj/rho)
gamma_residj <- as.vector(byj - matrixVectorMultiply(bXj,thetaj) - A_gammaj%*%gamma_centerj)
gamma_hessianj <- Matrix::forceSymmetric(rho*Cinvj + tcrossprod(A_gammaj))
gamma_middlej <- Matrix::solve(gamma_hessianj, gamma_residj)
gammaj=as.vector(gamma_centerj + crossprod(A_gammaj,gamma_middlej))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
if(jiter>4) errorj=norm(thetaj-theta_prevj,"2")
if(errorj<max.eps) break
}
ThetaList[j, ] <- thetaj
GammaList[j, ] <- gammaj
resampling.retries=0
j=j+1
}, error = function(e) {
resampling.retries <<- resampling.retries+1
if(resampling.retries>20){
stop("MRBEE_IPOD_SuSiE resampling failed repeatedly at ", resampling.stage, ": ", conditionMessage(e))
}
indicator <<- TRUE
})
if (indicator) {
next
}
}
if(verbose) close(pb)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Resampling ends: ",time_to_print," secs\n"))
}
theta.se=colSDMAD(ThetaList)
theta.cov=covmad(ThetaList)
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$theta.pip=colMeans(ThetaList!=0)
A$Bic=Bbic
A$reliability.adjust=r
A$susie.theta=fit.theta
A$thetalist=ThetaList
A$gammalist=GammaList
A$Diff=Diff
A$Group_Penalty=Diff_matrix
return(A)
}
