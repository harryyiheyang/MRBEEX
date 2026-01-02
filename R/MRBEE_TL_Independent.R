MRBEE_TL_Independent=function(by,bX,byse,bXse,Rxy,theta.source,theta.source.cov,tauvec=seq(3,30,3),admm.rho=3,Lvec=c(1:6),ebic.delta=1,ebic.gamma=2,transfer.coef=1,susie.iter=200,pip.thres=0.5,pip.min=0.1,cred.pip.thres=0.95,max.iter=50,max.eps=1e-4,reliability.thres=0.8,ridge.diff=100,sampling.time=100,sampling.iter=10,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,coverage.causal=0.95,LDSC=NULL,Omega=NULL,estimate_residual_method="MoM",sampling.strategy="bootstrap",standardize=T){
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
r=reliability.adj(bX,bXse%*%diag(sqrt(diag(Rxy[1:p,1:p]))),thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=LDSC,Omega=Omega)
Rxyall=biasterm(RxyList=RxyList,c(1:n))
br=as.vector(by-bX%*%theta.source)
BtB=t(bX)%*%bX
Diff_matrix=diag(p)*0
if(group.penalize==T){
Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
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
gamma1=u=gamma*0
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
br.complement=c(br-bX%*%delta.complement-gamma)
addbias=matrixVectorMultiply(Rxysum[1:p,1:p],theta.source+delta.complement)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2+Diff_matrix
Xty=matrixVectorMultiply(t(bX),br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum((br.complement)^2)
tryCatch({
fit.susie=susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=0.9,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie=susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[i],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
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
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]+Diff_matrix[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
gamma1=mcp(gamma+u/admm.rho,tauvec[v]/admm.rho)
u=u+admm.rho*(gamma-gamma1)
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>5){
error=sqrt(sum((delta-delta1)^2))
}
}
e=by-matrixVectorMultiply(bX,theta)-gamma
vare=sum(e^2)/(length(indvalid)-length(inddelta))
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
delta=theta.source-theta
fit.susie=NULL
error=2
iter=0
gamma1=u=gamma*0
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
br.complement=c(br-bX%*%delta.complement-gamma)
addbias=matrixVectorMultiply(Rxysum[1:p,1:p],theta.source+delta.complement)
XtX=BtB-Rxysum[1:p,1:p]
XtX=t(XtX)/2+XtX/2+Diff_matrix
Xty=matrixVectorMultiply(t(bX),br.complement)-Rxysum[1+p,1:p]+addbias
yty=sum((br.complement)^2)
tryCatch({
fit.susie=susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,residual_variance_lowerbound=0.9,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie=susie_ss(XtX=XtX,Xty=Xty,yty=yty,L=Lvec[istar],n=length(indvalid),estimate_prior_method="EM",residual_variance=1,model_init=fit.susie,max_iter=susie.iter,estimate_residual_variance=F,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
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
xtx=XtX[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=xty/xtx
}
if(length(inddelta)>1){
xtx=XtX[inddelta,inddelta]+ridge.diff*Diff[inddelta,inddelta]+Diff_matrix[inddelta,inddelta]
xty=Xty[inddelta]
delta.latent[inddelta]=c(solve(xtx)%*%xty)
}
delta=delta.latent+delta.complement
theta=delta+theta.source
gamma=(by-matrixVectorMultiply(bX,theta)-u+admm.rho*gamma1)/(1+admm.rho)
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
res=gamma*byse1
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
if (sampling.strategy == "bootstrap") {
indj <- sample(1:n, size = n, replace = TRUE)
} else {
indj <- sample(1:n, size = 0.5 * n, replace = FALSE)
}
nj=length(indj)
bXj=bX[indj,]
byj=by[indj]
bXsej=bXse[indj,]
bysej=byse[indj]
brj=br[indj]
thetaj=theta*runif(length(theta),0.95,1.05)
RxyListj=RxyList[indj,,]
Rxyallj=biasterm(RxyList=RxyListj,c(1:nj))
gammaj=gamma[indj]
indvalidj=which(gammaj==0)
deltaj=theta.source-thetaj
BtBj=matrixMultiply(t(bXj),bXj)
gamma1j=uj=gammaj*0
errorj=1
if(sampling.strategy=="bootstrap"){
fit.susiej=fit.susie
}else{
fit.susiej=NULL
}

for(jiter in 1:sampling.iter){
theta_prevj=thetaj
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
XtXj=t(XtXj)/2+XtXj/2+Diff_matrix
Xtyj=matrixVectorMultiply(t(bXj),br.complementj)-Rxysumj[1+p,1:p]+addbiasj
ytyj=sum((br.complementj)^2)
tryCatch({
fit.susiej=susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=ifelse(jiter==1,1000,30),residual_variance_lowerbound=0.9,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susiej=susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,L=Lvec[istar],n=length(indvalidj),estimate_prior_method="EM",residual_variance=1,model_init=fit.susiej,max_iter=ifelse(jiter==1,1000,30),estimate_residual_variance=F,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
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
xtxj=XtXj[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=xtyj/xtxj
}
if(length(inddeltaj)>1){
xtxj=XtXj[inddeltaj,inddeltaj]+ridge.diff*Diffj[inddeltaj,inddeltaj]+Diff_matrix[inddeltaj,inddeltaj]
xtyj=Xtyj[inddeltaj]
delta.latentj[inddeltaj]=c(solve(xtxj)%*%xtyj)
}
deltaj=delta.latentj+delta.complementj
thetaj=deltaj+theta.source
gammaj=(byj-matrixVectorMultiply(bXj,thetaj)-uj+admm.rho*gamma1j)/(1+admm.rho)
gamma1j=mcp(gammaj+uj/admm.rho,tauvec[vstar]/admm.rho)
uj=uj+admm.rho*(gammaj-gamma1j)
gammaj=gammaj*(gamma1j!=0)
if(jiter>3) errorj=norm(thetaj-theta_prevj,"2")
if(errorj<max.eps) break
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
