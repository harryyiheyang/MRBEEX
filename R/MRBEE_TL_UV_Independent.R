MRBEE_TL_UV_Independent=function(by,bX,byse,bXse,Rxy,theta.source,theta.source.cov,tauvec=seq(3,30,3),admm.rho=3,ebic.delta=1,ebic.gamma=2,transfer.coef=1,susie.iter=200,pip.thres=0.3,max.iter=50,max.eps=1e-4,reliability.thres=0.8,ridge.diff=100,sampling.time=100,sampling.iter=10,LDSC=NULL,Omega=NULL){
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
BtB=sum(bX^2)
r=reliability.adj.uv(bX,bXse*sqrt(Rxy[1,1]),thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,LDSC=LDSC,Omega=Omega)
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
by.complement=as.vector(by-bX*theta.complement-gamma)
XtX=BtB-sum(bXse[indvalid]^2*Rxy[1,1])
Xty=sum(bX*by.complement)-Rxy[1,2]*sum(bXse[indvalid]*byse[indvalid])+sum(bXse[indvalid]^2*theta.complement)
yty=sum(by.complement^2)
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
gamma=as.vector(by-bX*theta-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
dftheta=as.numeric(fit.susie$pip>pip.thres)
e=as.vector(by-bX*theta-as.vector(gamma))
vare=sum(e^2)/(length(indvalid)-dftheta)
Bic[v]=log(vare)+log(n)/n*dftheta+(1+ebic.gamma)*log(n)*(n-length(indvalid))/n
Btheta[v]=theta
Bgamma[v,]=gamma1
e=as.vector(by-bX*theta.source-as.vector(gamma))
vare=sum(e^2)/(sum(gamma==0))
Bic_direct[v]=log(vare)+(1+ebic.gamma)*log(n)*sum(gamma!=0)/n
Bgamma_direct[v,]=gamma
}
#########################################################################
s1=min(Bic)
s2=min(Bic_direct)
if(s1<=s2){
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
by.complement=as.vector(by-bX*theta.complement-gamma)
XtX=BtB-sum(bXse[indvalid]^2*Rxy[1,1])
Xty=sum(bX*by.complement)-Rxy[1,2]*sum(bXse[indvalid]*byse[indvalid])+sum(bXse[indvalid]^2*theta.complement)
yty=sum(by.complement^2)
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
gamma=as.vector(by-bX*theta-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[vstar]/admm.rho)
gamma=gamma*(gamma1!=0)
u=u+admm.rho*(gamma-gamma1)
iter=iter+1
if(iter>5){
error=sqrt(sum((theta-theta1)^2))
}
}
############################### inference #########################
res=gamma*byse1
names(res)=rownames(bX)
ThetaList=c(1:sampling.time)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time) {
setTxtProgressBar(pb, j)
indicator <- FALSE
tryCatch({
indj=sample(n,0.5*n,replace=F)
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
BtBj <- sum(bXj^2)
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
if(length(indvalidj)==nj){
Rxysumj=Rxyallj
}else{
Rxysumj=Rxyallj-biasterm(RxyList=RxyListj,setdiff(1:nj,indvalidj))
}
theta.complementj=ifelse(which.min(c(abs(thetaj),abs(thetaj-theta.source)))==1,0,theta.source)
by.complementj=as.vector(byj-bXj*theta.complementj-gammaj)
XtXj=BtBj-sum(bXsej[indvalidj]^2*Rxy[1,1])
Xtyj=sum(bX*by.complementj)-Rxy[1,2]*sum(bXsej[indvalidj]*bysej[indvalidj])+sum(bXsej[indvalidj]^2*theta.complementj)
ytyj=sum(by.complementj^2)
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
gammaj=as.vector(byj-bXj*thetaj-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
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
theta.cov=var(ThetaList)*n/length(indvalid)
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
A$susie.theta=fit.susie
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
