MRBEE_IPOD_SuSiE=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),Lvec=c(1:min(10,nrow(bX))),pip.thres=0.5,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=1,ebic.gamma=2,reliability.thres=0.8,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=10,theta.ini=F,gamma.ini=F,ridge.diff=1e5,verbose=T,pip.min=0.1,cred.pip.thres=0.95,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10){
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
if(LD[1]!="identity"){
isLD=T
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
Thetarho=solve(LD+rho*diag(m))
}else{
isLD=F
LD=Theta=TC=Matrix(diag(m),sparse=T)
RC=diag(m)
bXinv=tilde.X=bX
Bt=t(bX)
BtB=matrixMultiply(t(bX),bX)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
tilde.y=by
Thetarho=diag(m)*1/(1+rho)
Thetarho=Matrix(Thetarho,sparse=T)
}
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
Diff_matrix=diag(p)*0
if(group.penalize==T){
  Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
############################ Initial Estimate #######################
if(theta.ini[1]==F){
if(length(tilde.y)<2000){
fit0=varbvs(X=RC,Z=tilde.X,y=tilde.y,verbose=F,maxiter=100)
gamma.ini=gamma.ini1=fit0$beta*(fit0$pip>0.5)
theta.ini=theta.ini1=fit0$beta.cov[-1]
}else{
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak",pv.thres=0.01)
gamma.ini=gamma.ini1=fit0$gamma/byse1
theta.ini=theta.ini1=fit0$theta
}
}else{
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
}
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
empirical.variance=5
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB+Diff_matrix-Rxysum[1:p,1:p]
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[v],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=susie.iter,s_init=fit.theta)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.min)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.theta)$vars,m/dBtB,theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]+Diff_matrix[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
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
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta[,j,v]=theta
Bgamma[,j,v]=gamma1
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
empirical.variance=5
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-as.vector(LD%*%gamma)
XtX=BtB+Diff_matrix-Rxysum[1:p,1:p]
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,p+1]
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.min)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.theta)$vars,m/dBtB,theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]+Diff_matrix[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
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
gamma=gamma1
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
cat("Bootstrapping starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
LDj=LD[indj,indj]
Thetaj <- Theta[indj,indj]
Thetarhoj <- Thetarho[indj,indj]
Btj <- as.matrix(t(bX[indj, ]) %*% Thetaj)
BtBj <- Btj%*%bX[indj, ]
BtBj=(t(BtBj)+BtBj)/2
dBtBj=diag(BtBj)
thetaj=theta*runif(1,0.95,1.05)
gammaj=gamma1j=gamma*runif(1,0.975,1.025)
deltaj=gammaj*0
for(jiter in 1:sampling.iter){
indvalidj <- which(gamma1j==0)
indvalidj <- intersect(indvalidj, indj)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
res.thetaj=by[indj]-as.vector(LD[indj,indj]%*%gammaj[indj])
XtXj=BtBj+Diff_matrix-Rxysumj[1:p,1:p]
Xtyj=matrixVectorMultiply(Btj,res.thetaj)-Rxysumj[1:p,1+p]
ytyj=sum(res.thetaj*(Thetaj%*%res.thetaj))
fit.thetaj=susie_suff_stat(XtX=BtBj,Xty=Xtyj,yty=ytyj,n=length(indvalidj),L=Lvec[vstar],estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=sampling.iter,s_init=fit.theta)
thetaj=coef.susie(fit.thetaj)[-1]*(fit.thetaj$pip>max(pip.min/sqrt(2),0.1))
theta.csj=group.pip.filter(pip.summary=summary(fit.thetaj)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=max(pip.thres/sqrt(2),0.1))
pip.alivej=theta.csj$ind.keep
thetaj[-pip.alivej]=0
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.thetaj)$vars,m/dBtBj,thetaj)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
XtXj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]+ridge.diff*Diffj[indthetaj,indthetaj]+Diff_matrix[indthetaj,indthetaj]
Xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=c(solve(XtXj)%*%Xtyj)
}
if((norm(thetaj, "2") / norm(theta.ini1, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini1, "2")
}
if(isLD){
gammaj[indj]=as.vector(Thetarhoj%*%(by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj]))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}else{
gammaj[indj]=as.vector((by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj]))/(1+rho)
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
}
ThetaList[j, ] <- thetaj
GammaList[j, ] = gamma1j
}
close(pb)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Bootstrapping ends: ",time_to_print," secs\n"))
}
theta.se=colSD(ThetaList)*sqrt((m-length(indtheta))/(m-length(indtheta)-length(indgamma)))
theta.cov=cov(ThetaList)*(m-length(indtheta))/(m-length(indtheta)-length(indgamma))
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
A$theta.pratt=getPratt(bX=bX,by=by,bXse=bXse,byse=byse,Theta=Theta,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma,Theta=Theta,LD=LD)
A$Diff=Diff
A$Group_Penalty=Diff_matrix
return(A)
}
