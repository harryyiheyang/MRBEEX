MRBEE_IPOD_SuSiE=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),Lvec=c(1:min(10,nrow(bX))),pip.thres=0.5,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=1,ebic.gamma=2,reliability.thres=0.8,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=10,theta.ini=F,gamma.ini=F,ridge.diff=1e5,verbose=T,pip.min=0.1,cred.pip.thres=0.95,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,coverage.causal=0.95,LDSC=NULL,Omega=NULL,estimate_residual_variance=T,prob_shrinkage_coef=0.5,prob_shrinkage_size=4,estimate_residual_method="MoM",sampling.strategy="bootstrap",standardize=T){
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
BtB=matrixMultiply(t(bX),bX)
BtB=t(BtB)/2+BtB/2
dBtB=diag(BtB)
tilde.y=by
Thetarho=diag(m)*1/(1+rho)
Thetarho=Matrix(Thetarho,sparse=T)
}
r=reliability.adj(bX,bXse/sqrt(diag(Rxy[1:p,1:p])),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=LDSC,Omega=Omega)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
Diff_matrix=diag(p)*0
if(group.penalize==T){
Diff_matrix=group.diff*generate_group_matrix(group_index=group.index,COV=BtB)
}
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=susie_ss(XtX=BtB+Diff_matrix,Xty=matrixVectorMultiply(t(bXinv),by),yty=sum(by*(Theta%*%by)),L=10,n=m,coverage = coverage.causal,standardize=standardize)
theta.ini=theta.ini1=coef(fit0)[-1]
gamma.ini=gamma.ini1=by*0
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
XtX=BtB+Diff_matrix-Rxysum[1:p,1:p]
XtX=XtX/2+t(XtX)/2
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
tryCatch({
fit.theta=susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],estimate_prior_method="EM",max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.theta=susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[v],estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.min)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.theta)$vars,m/dBtB,theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]
xty=Xty[indtheta]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]
Xty=Xty[indtheta]
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
XtX=BtB+Diff_matrix-Rxysum[1:p,1:p]
XtX=XtX/2+t(XtX)/2
Xty=matrixVectorMultiply(Bt,res.theta)-Rxysum[1:p,1+p]
yty=sum(res.theta*(Theta%*%res.theta))
tryCatch({
fit.theta=susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],estimate_prior_method="EM",max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.theta=susie_ss(XtX=XtX,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=susie.iter,model_init=fit.theta,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.min)
theta.cs=group.pip.filter(pip.summary=summary(fit.theta)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive=theta.cs$ind.keep
theta[-pip.alive]=0
indtheta=which(theta!=0)
Diff=generate_block_matrix(summary(fit.theta)$vars,m/dBtB,theta)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]
xty=Xty[indtheta]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]+ridge.diff*Diff[indtheta,indtheta]
Xty=Xty[indtheta]
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
cat("Bootstrapping starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
cluster.index <- as.integer(factor(cluster.index))
cluster_prob <- cluster_prob(cluster.index,LD,alpha=prob_shrinkage_coef,group_size=prob_shrinkage_size)
if(isLD) {
cluster_cache <- precompute_cluster_blocks(
bX = bX,
bXse = bXse,
by = by,
byse = byse,
LD = LD,
Theta = Theta,
Thetarho = Thetarho,
cluster.index = cluster.index
)
}
j=1
while(j<=sampling.time){
indicator <- FALSE
setTxtProgressBar(pb, j)
tryCatch({
if(isLD==T){
if (sampling.strategy == "bootstrap") {
cluster.sampling <- sample(1:max(cluster.index),
                           size = max(cluster.index),
                           replace = TRUE,
                           prob = cluster_prob)
} else {
cluster.sampling <- sample(1:max(cluster.index),
                           size = 0.5 * max(cluster.index),
                           replace = FALSE,
                           prob = cluster_prob)
}
cluster.sampling=sort(cluster.sampling)
sampled_blocks <- cluster_cache[cluster.sampling]
indj <- unlist(lapply(sampled_blocks, function(b) b$idx))
indj <- sort(indj)
mj <- length(indj)
LDj <- bdiag(lapply(sampled_blocks, function(b) b$LD))
Thetaj <- bdiag(lapply(sampled_blocks, function(b) b$Theta))
Thetarhoj <- bdiag(lapply(sampled_blocks, function(b) b$Thetarho))
bXj <- do.call(rbind, lapply(sampled_blocks, function(b) b$bX))
bXsej <- do.call(rbind, lapply(sampled_blocks, function(b) b$bXse))
byj <- unlist(lapply(sampled_blocks, function(b) b$by))
bysej <- unlist(lapply(sampled_blocks, function(b) b$byse))
BtBj <- Reduce('+', lapply(sampled_blocks, function(b) b$BtB))
BtBj <- (t(BtBj) + BtBj) / 2
dBtBj <- diag(BtBj)
Btj <- do.call(cbind, lapply(sampled_blocks, function(b) b$Bt))
}else{
if (sampling.strategy == "bootstrap") {
indj <- sample(1:m, size = m, replace = TRUE)
} else {
indj <- sample(1:m, size = 0.5 * m, replace = FALSE)
}
indj=sort(indj)
LDj=Matrix(diag(length(indj)),sparse=T)
bXj=bX[indj,]
bXsej=bXse[indj,]
byj=by[indj]
bysej=byse[indj]
Btj <- t(bXj)
BtBj=matrixMultiply(Btj,bXj)
BtBj=(t(BtBj)+BtBj)/2
dBtBj=diag(BtBj)
Thetaj=LDj
mj=length(indj)
}
thetaj=theta*runif(p,0.95,1.05)
gammaj=gamma1j=gamma[indj]*runif(1,0.975,1.025)
deltaj=gammaj*0
errorj=1
if(sampling.strategy=="bootstrap"){
fit.thetaj=fit.theta
}else{
fit.thetaj=NULL
}

for(jiter in 1:sampling.iter){
theta_prevj=thetaj
indvalidj <- which(gamma1j==0)
if(length(indvalidj)<(0.55*length(indj))){
indvalidj=sample(1:length(indj),0.6*length(indj))
gamma1j[indvalidj]=gammaj[indvalidj]=0
}
Rxysumj <- biasterm(RxyList = RxyList, indj[indvalidj])
res.thetaj=byj-as.vector(LDj%*%gammaj)
XtXj=BtBj+Diff_matrix/2-Rxysumj[1:p,1:p]
XtXj=XtXj/2+t(XtXj)/2
Xtyj=matrixVectorMultiply(Btj,res.thetaj)-Rxysumj[1:p,p+1]
ytyj=sum(res.thetaj*(Thetaj%*%res.thetaj))
tryCatch({
fit.thetaj=susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,n=mj,L=Lvec[vstar],estimate_prior_method="EM",max_iter=ifelse(jiter==1,1000,30),model_init=fit.thetaj,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.thetaj=susie_ss(XtX=XtXj,Xty=Xtyj,yty=ytyj,n=mj,L=Lvec[vstar],estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),max_iter=ifelse(jiter==1,1000,30),model_init=fit.thetaj,coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
thetaj=coef.susie(fit.thetaj)[-1]*(fit.thetaj$pip>max(pip.min/sqrt(2),0.1))
theta.csj=group.pip.filter(pip.summary=summary(fit.thetaj)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=max(pip.thres/sqrt(2),0.1))
pip.alivej=theta.csj$ind.keep
thetaj[-pip.alivej]=0
indthetaj=which(thetaj!=0)
Diffj=generate_block_matrix(summary(fit.thetaj)$vars,m/dBtBj,thetaj)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
XtXj=XtXj[indthetaj,indthetaj]+ridge.diff*Diffj[indthetaj,indthetaj]
Xtyj=Xtyj[indthetaj]
thetaj[indthetaj]=c(solve(XtXj)%*%Xtyj)
}
if((norm(thetaj, "2") / norm(theta.ini1, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini1, "2")
}
if(isLD){
gammaj=as.vector(Thetarhoj%*%(byj-matrixVectorMultiply(bXj,thetaj)-deltaj+rho*gamma1j))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}else{
gammaj=as.vector((byj-matrixVectorMultiply(bXj,thetaj)-deltaj+rho*gamma1j))/(1+rho)
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
gammaj=gammaj*(gamma1j!=0)
if(jiter>4) errorj=norm(thetaj-theta_prevj,"2")
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
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Bootstrapping ends: ",time_to_print," secs\n"))
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
