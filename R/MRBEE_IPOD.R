MRBEE_IPOD=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),tauvec=seq(3,50,by=5),max.iter=100,max.eps=0.001,ebic.gamma=1,reliability.thres=0.5,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=5,theta.ini=F,gamma.ini=F,ebic.theta=1,verbose=T,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,sampling.strategy="bootstrap"){
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
cluster.index <- as.integer(factor(cluster.index))
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
Thetarho=solve(LD+rho*diag(m))
}else{
isLD=F
LD=Theta=TC=Matrix(diag(m),sparse=T)
RC=diag(m)
bXinv=tilde.X=bX
Bt=t(bX)
BtB=matrixMultiply(bX,bX,transA=TRUE)
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
############################ Initial Estimate #######################
theta.ini.missing=length(theta.ini)==1&&((is.logical(theta.ini)&&!theta.ini)||
                  (is.numeric(theta.ini)&&is.finite(theta.ini)&&theta.ini==0))
if(theta.ini.missing){
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak",pv.thres=0.01,group.penalize=group.penalize,group.index=group.index,group.diff=group.diff)
gamma.ini=gamma.ini1=fit0$gamma/byse
theta.ini=theta.ini1=fit0$theta
}else{
if(length(theta.ini)!=p){
stop("theta.ini must have one value for each exposure.")
}
gamma.ini.missing=length(gamma.ini)==1&&is.logical(gamma.ini)&&!gamma.ini
if(gamma.ini.missing){
gamma.ini=rep(0,m)
}else{
if(length(gamma.ini)!=m){
stop("gamma.ini must have one value for each variant.")
}
gamma.ini=gamma.ini/byse
}
theta.ini1=theta.ini
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
Btheta=array(0,c(p,w))
Bgamma=array(0,c(m,w))
Bbic=c(1:w)
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1
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
g=matrixVectorMultiply(Bt,as.vector(by-LD%*%gamma))-Rxysum[1:p,p+1]
theta=c(CppMatrix::matrixSolve(BtB-Rxysum[1:p,1:p]+Diff_matrix,g))
if(theta.ini.norm>0&&(norm(theta,"2")/theta.ini.norm)>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*theta.ini.norm
}
########################### update gamma ############################
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
Btheta[,j]=theta
Bgamma[,j]=gamma
df1=sum(gamma1!=0)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
rss=sum(res*(Theta%*%res))/(m-df1-p)
Bbic[j]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+p*(log(m)+ebic.theta*log(p))
}
Bbic=Bbic/m
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Estimation ends: ",time_to_print," secs\n"))
}
######################## Inference #################################
t1=Sys.time()
jstar=which.min(Bbic)
theta=Btheta[,jstar]
gamma=Bgamma[,jstar]
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma!=0)
indvalid=which(gamma==0)
Rxysum=biasterm(RxyList=RxyList,indvalid=indvalid)

if(sampling.time>0){
ThetaList=matrix(0,sampling.time,p)
colnames(ThetaList)=colnames(bX)
GammaList=matrix(0,sampling.time,m)
colnames(GammaList)=rownames(bX)
cat("Resampling starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time){
indicator <- FALSE
setTxtProgressBar(pb, j)
tryCatch({
if(isLD==T){
indj <- sort(sample.int(m, size = max(1L, floor(0.5 * m)), replace = FALSE))
mj <- length(indj)
LDj <- LD[indj, indj, drop = FALSE]
Thetaj <- solve(LDj)
Thetarhoj <- solve(LDj + rho * diag(mj))
bXj <- bX[indj,,drop=FALSE]
bXsej <- bXse[indj,,drop=FALSE]
byj <- by[indj]
bysej <- byse[indj]
Btj <- as.matrix(t(bXj)%*%Thetaj)
BtBj <- matrixMultiply(Btj,bXj)
BtBj <- (t(BtBj) + BtBj) / 2
dBtBj <- diag(BtBj)
}else{
if (is_bootstrap_sampling(sampling.strategy)) {
indj <- sample(1:m, size = m, replace = TRUE)
} else {
indj <- sample(1:m, size = 0.5 * m, replace = FALSE)
}
indj=sort(indj)
LDj=Matrix(diag(length(indj)),sparse=T)
Thetaj=LDj
bXj=bX[indj,,drop=FALSE]
bXsej=bXse[indj,,drop=FALSE]
byj=by[indj]
bysej=byse[indj]
Btj <- t(bXj)
BtBj=matrixMultiply(Btj,bXj)
BtBj=(t(BtBj)+BtBj)/2
dBtBj=diag(BtBj)
}
Rxyallj <- biasterm(RxyList = RxyList, indj)
thetaj=theta*runif(p,0.95,1.05)
gammaj=gamma1j=gamma[indj]*runif(1,0.975,1.025)
deltaj=0*gammaj
errorj=1
for(jiter in 1:sampling.iter){
theta_prevj=thetaj
indvalidj <- which(gamma1j==0)
if(length(indvalidj)<(0.55*length(indj))){
indvalidj=sample(1:length(indj),0.6*length(indj))
gamma1j[indvalidj]=gammaj[indvalidj]=0
}
invalidj <- which(gamma1j!=0)
if(length(invalidj)<length(indvalidj)){
Rxysumj <- Rxyallj-biasterm(RxyList = RxyList, indj[invalidj])
}else{
Rxysumj <- biasterm(RxyList = RxyList, indj[indvalidj])
}
g <- matrixVectorMultiply(Btj, byj - as.vector(LDj%*%gammaj)) - Rxysumj[1:p, p + 1]
thetaj <- c(CppMatrix::matrixSolve(BtBj - Rxysumj[1:p, 1:p]+Diff_matrix/2, g))
if(theta.ini.norm>0&&(norm(thetaj,"2")/theta.ini.norm)>maxdiff){
thetaj=thetaj/norm(thetaj,"2")*maxdiff*theta.ini.norm
}
if(isLD){
gammaj=as.vector(Thetarhoj%*%(byj-matrixVectorMultiply(bXj,thetaj)-deltaj+rho*gamma1j))
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}else{
gammaj=gamma1j=mcp(by[indj]-matrixVectorMultiply(bXj,thetaj),tauvec[jstar])
deltaj=deltaj+rho*(gammaj-gamma1j)
}
gammaj=gammaj*(gamma1j!=0)
if(jiter>4) errorj=norm(thetaj-theta_prevj,"2")
if(errorj<max.eps) break
}
ThetaList[j, ] <- thetaj
j=j+1
}, error = function(e) {
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
})
if (indicator) {
next  # Retry the current iteration
}
}
close(pb)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Resampling ends: ",time_to_print," secs\n"))
}
theta.se=colSDMAD(ThetaList)
theta.cov=covmad(ThetaList)
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)
}else{
adjf=m/(length(indvalid)-p)
if(sum(gamma!=0)>0){
bZ=as.matrix(cbind(bX,LD[,which(gamma!=0)]))
}else{
bZ=as.matrix(bX)
}
H=matrixMultiply(bZ,as.matrix(Theta%*%bZ),transA=TRUE)
H[1:p,1:p]=H[1:p,1:p]-Rxysum[1:p,1:p]+Diff_matrix
H=(H+t(H))/2
e=res
E=-as.matrix(Theta%*%bZ)*e
for(i in 1:length(indvalid)){
E[indvalid[i],1:p]=E[indvalid[i],1:p]+RxyList[indvalid[i],1:p,p+1]-as.vector(RxyList[indvalid[i],1:p,1:p]%*%theta)
}
V=crossprod(E)
npar=ncol(bZ)
Dmcp=rep(0,npar)
if(length(indgamma)>0){
ig=p+c(1:length(indgamma))
gabs=abs(gamma[indgamma])
Dmcp[ig]=ifelse(gabs<3*tauvec[jstar]/rho,-rho/3,0)
}
H.mcp=H+diag(Dmcp,npar)
Hinv.mcp=CppMatrix::matrixSolve(H.mcp,diag(npar))
covtheta=(Hinv.mcp%*%V%*%Hinv.mcp)*adjf
theta.se=sqrt(diag(covtheta[1:p,1:p]))
theta.cov=covtheta[1:p,1:p]
ThetaList=GammaList=NULL
}


A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini*byse1
A$reliability.adjust=r
A$thetalist=ThetaList
A$gammalist=GammaList
A$Group_Penalty=Diff_matrix
return(A)
}
