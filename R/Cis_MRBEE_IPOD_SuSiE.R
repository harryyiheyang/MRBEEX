Cis_MRBEE_IPOD_SuSiE=function(by,bX,byse,bXse,LD,Rxy,Lvec=c(1:min(10,nrow(bX))),pip.thres=0.2,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=1,ebic.gamma=2,reliability.thres=0.8,rho=2,theta.ini=F,gamma.ini=F,ridge=ridge,pleiotropy.rm=NULL,sandwich=F,sampling.time=300,sampling.iter=15,coverage.causal=0.95,xQTLfitList=NULL){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
LD=as.matrix(LD)
if(m<=1000){
Theta=matrixInverse(LD)
}else{
Theta=positiveinv(LD,min.eps=0)
}
bXinv=matrixMultiply(matrixInverse(LD*0.975+diag(m)*0.025),bX)
Bt=matrixMultiply(t(bX),Theta)
BtB=matrixMultiply(Bt,bX)
BtB=(t(BtB)+BtB)/2
dBtB=sqrt(diag(BtB)/m)
ridge=ridge*m
Thetarho=solve(LD+rho*diag(m))
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=susie_suff_stat(XtX=BtB,Xty=c(t(bXinv)%*%by),yty=sum(by*(Theta%*%by)),n=m,intercept=F,estimate_prior_method="EM",coverage=coverage.causal)
theta.ini=coef.susie(fit0)[-1]
gamma.ini=gamma.ini1=by*0
theta.ini1=theta.ini
}else{
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
}
############################## Tuning Parameter ######################
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
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-matrixVectorMultiply(LD,gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*matrixVectorMultiply(Theta,res.theta))
fit.theta=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],residual_variance=1,estimate_residual_variance=T,estimate_prior_method="EM",intercept=F,max_iter=susie.iter,standardize=F,s_init=fit.theta,coverage=coverage.causal)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=0.95,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
Diff=generate_block_matrix(summary(fit.theta)$vars,1/dBtB,theta)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge*Diff[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma[pleiotropy.rm]=0
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
gamma1[pleiotropy.rm]=0
delta=delta+rho*(gamma-gamma1)
delta[pleiotropy.rm]=0
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
res=by-matrixVectorMultiply(bX,theta)-matrixVectorMultiply(LD,gamma)
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
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-matrixVectorMultiply(LD,gamma)
XtX=BtB
Xty=matrixVectorMultiply(Bt,res.theta)
yty=sum(res.theta*(Theta%*%res.theta))
fit.theta=susie_suff_stat(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],residual_variance=1,estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=susie.iter*10,standardize=F,s_init=fit.theta,coverage=coverage.causal)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=0.95,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
Diff=generate_block_matrix(summary(fit.theta)$vars,1/dBtB,theta)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]+ridge*Diff[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
gamma=matrixVectorMultiply(Thetarho,by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)
gamma[pleiotropy.rm]=0
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
gamma1[pleiotropy.rm]=0
delta=delta+rho*(gamma-gamma1)
delta[pleiotropy.rm]=0
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
############################### inference #########################
gamma=gamma
theta=theta
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma1!=0)
indvalid=which(gamma1==0)
res=by-matrixVectorMultiply(bX,theta)-matrixVectorMultiply(LD,gamma)
if(length(indvalid)==m){
  Rxysum=Rxyall
}else{
  Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
adjf=m/(length(indvalid)-length(indtheta))

ThetaList=NULL
var_inf=0
var_error=1
if(sandwich==T){
if(length(indtheta)>0){
upsilon=by*0
var_inf=1e-2
var_error=1
for(vv in 1:30){
Hupsilon=solve(diag(m)/var_inf+LD/var_error)
upsilon=as.vector(Hupsilon%*%res)/var_error
df=sum(diag(Hupsilon))
var_inf=min((sum(upsilon^2)+df)/m,10)
res_inf=res-matrixVectorMultiply(LD,upsilon)
df=sum(diag(Hupsilon%*%LD))/var_error
var_error=sum(res_inf*(Theta%*%res_inf))/(m-df-length(indtheta)-length(indgamma))
var_error=max(1,var_error)
}

bZ=cbind(bX[,indtheta],as.matrix(LD[,indgamma]))
Hinv=matrixListProduct(list(t(bZ),Theta,bZ))
Hinv[1:length(indtheta),1:length(indtheta)]=Hinv[1:length(indtheta),1:length(indtheta)]-Rxysum[indtheta,indtheta]+ridge*Diff[indtheta,indtheta]
Hinv=positiveinv(as.matrix(Hinv))
Hinv1=matrixListProduct(list(t(bZ),Theta,var_error*LD+var_inf*LD%*%LD,Theta,bZ))
COV=Hinv%*%Hinv1%*%Hinv
theta.cov=diag(p)*0
theta.cov[indtheta,indtheta]=COV[1:length(indtheta),1:length(indtheta)]
theta.se=sqrt(diag(theta.cov))
}

if(length(indtheta)==0){
theta.cov=diag(p)*0
theta.se=rep(0,p)
}
}else{
############################### inference #########################
t1=Sys.time()
gamma=gamma
theta=theta
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma1!=0)
indvalid=which(gamma1==0)
res=by-matrixVectorMultiply(bX,theta)-matrixVectorMultiply(LD,gamma)
ThetaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(bX)
cat("Bootstrapping starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time){
indicator <- FALSE
setTxtProgressBar(pb, j)
tryCatch({
bXj=bX0j=bX*0
for(ii in 1:p){
rsamples=susie_effect_resampling(LD=LD,alpha=xQTLfitList[[ii]]$alpha,mu=xQTLfitList[[ii]]$mu,mu2=xQTLfitList[[ii]]$mu2,sampling=1,method="probabilistic")
bXj[,ii]=rsamples$bx
bX0j[,ii]=rsamples$bx0
}
bXj=bXj*byseinv
bX0j=bX0j*byseinv
emptyy=fix_empty_resamples(bXj,bX0j,dBtB)
bXj=emptyy$bXj
bX0j=emptyy$bX0j
Btj <- matrixMultiply(t(bXj),Theta)
BtBj <- matrixMultiply(Btj,bXj)
BtBj=(t(BtBj)+BtBj)/2
dBtBj=diag(BtBj)/m
thetaj=theta*runif(1,0.95,1.05)
gammaj=gamma1j=gamma
byj=by
deltaj=gammaj*0
fit.thetaj=fit.theta
errorj=1
for(jiter in 1:sampling.iter){
theta_prevj=thetaj
indvalidj <- which(gamma1j==0)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
res.thetaj=byj-matrixVectorMultiply(LD,gammaj)
XtXj=BtBj
Xtyj=matrixVectorMultiply(Btj,res.thetaj)
ytyj=sum(res.thetaj*(Theta%*%res.thetaj))
fit.thetaj=susie_suff_stat(XtX=XtXj,Xty=Xtyj,yty=ytyj,n=m,L=Lvec[vstar],estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=sampling.iter,s_init=fit.thetaj,coverage=max(0.5,coverage.causal*0.9))
thetaj=coef.susie(fit.thetaj)[-1]*(fit.thetaj$pip>pip.thres)
theta.csj=group.pip.filter(pip.summary=summary(fit.thetaj)$var,xQTL.cred.thres=0.95,xQTL.pip.thres=pip.thres)
pip.alivej=theta.csj$ind.keep
thetaj[-pip.alivej]=0
Diffj=generate_block_matrix(summary(fit.thetaj)$vars,1/dBtBj,thetaj)
indthetaj=which(thetaj!=0)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,1+p]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
XtXj=XtXj[indthetaj,indthetaj]+ridge*Diffj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
Xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,1+p]
thetaj[indthetaj]=c(solve(XtXj)%*%Xtyj)
}
gammaj=as.vector(Thetarho%*%(byj-matrixVectorMultiply(bXj,thetaj)-deltaj+rho*gamma1j))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
if(jiter>2) errorj=norm(thetaj-theta_prevj,"2")
if(errorj<max.eps) break
}
ThetaList[j, ] <- thetaj
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
theta.se=colSDMAD(ThetaList)*sqrt((m-length(indtheta))/(m-length(indtheta)-length(indgamma)))
theta.cov=spearmancov(ThetaList)*(m-length(indtheta))/(m-length(indtheta)-length(indgamma))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}
A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$theta.z=A$theta/A$theta.se
A$theta.pip=fit.theta$pip
A$theta.list=ThetaList
A$Bic=Bbic
A$reliability.adjust=r
A$susie.theta=fit.theta
A$theta.pratt=getPratt(bX=bX,by=by,bXse=bXse,byse=byse,Theta=Theta,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma,Theta=Theta,LD=LD)
A$var.inf=var_inf
A$var.error=var_error
A$Diff=Diff
return(A)
}
