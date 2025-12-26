MRBEE_IPOD_UV=function(by,bX,byse,bXse,LD=LD,Rxy,cluster.index,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,ebic.gamma=1,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=5,theta.ini=F,gamma.ini=F,reliability.thres=0.8,LDSC=NULL,Omega=NULL,prob_shrinkage_coef=0.5,prob_shrinkage_size=4,sampling.strategy="bootstrap"){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=length(bX)
cluster.index <- as.integer(factor(cluster.index))
if(is.null(LDSC)==T){
LDSC=by*0
Omega=diag(2)*0
}
if(LD[1]!="identity"){
isLD=T
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
RC=as.matrix(TC%*%LD)
byinv=as.vector(Theta%*%by)
bXinv=as.vector(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.vector(TC%*%bX)
Bt=t(bXinv)
BtB=sum(Bt*bX)
Thetarho=solve(LD+rho*diag(m))
}else{
stop("please use MRBEE.IMPR.UV in MRBEE")
}
r=reliability.adj.uv(bX,bXse/sqrt(Rxy[1,1]),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=MRBEE.IMRP.UV(by=by,bx=bX,byse=byse,bxse=bXse,Rxy=Rxy)
gamma.ini=fit0$delta/byse
theta.ini=fit0$theta
}
############################## Tuning Parameter ######################
w=length(tauvec)
Btheta=c(1:w)
Bgamma=matrix(0,m,w)
Bbic=tauvec
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)<(0.55*m)) indvalid=sample(m,0.6*m)
Hinv=1/(BtB-sum(bXse[indvalid]^2)*Rxy[1,1]-sum(LDSC[indvalid]*byseinv[indvalid]^2)*Omega[1,1])
g=sum(Bt*(by-as.vector(LD%*%gamma)))-sum(bXse[indvalid])*Rxy[2,1]-sum(LDSC[indvalid]*byseinv[indvalid]^2)*Omega[1,2]
theta=g*Hinv
########################### update gamma ############################
gamma=as.vector(Thetarho%*%(by-bX*theta-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
gamma=gamma*(gamma1!=0)
iter=iter+1
if(iter>3){
error=abs(theta-theta1)
}
}
Btheta[j]=theta
Bgamma[,j]=gamma1
df1=sum(gamma1!=0)
res=by-bX*theta-as.vector(LD%*%gamma)
rss=sum(res*(Theta%*%res))/(m-df1-1)
Bbic[j]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+log(m)
}
Bbic=Bbic/m
######################## Inference #################################
jstar=last_min(Bbic)
theta=Btheta[jstar]
gamma=Bgamma[,jstar]
error=1
iter=1
theta=theta
gamma=gamma
names(gamma)=rownames(bX)
indgamma=which(gamma!=0)
indvalid=which(gamma==0)
res=by-bX*theta-as.vector(LD%*%gamma)
var_error=sum(res*(Theta%*%res))/(length(indvalid)-1)
if(sampling.time>0){
ThetaList=c(1:sampling.time)
cluster_prob <- cluster_prob(cluster.index,LD,alpha=prob_shrinkage_coef,group_size=prob_shrinkage_size)
cluster_cache <- precompute_cluster_blocks_uv(
bX = bX,
bXse = bXse,
by = by,
byse = byse,
LD = LD,
Theta = Theta,
Thetarho = Thetarho,
cluster.index = cluster.index,
rho = rho
)
for(j in 1:sampling.time){
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
LDj <- bdiag(lapply(sampled_blocks, function(b) b$LD))
Thetaj <- bdiag(lapply(sampled_blocks, function(b) b$Theta))
Thetarhoj <- bdiag(lapply(sampled_blocks, function(b) b$Thetarho))
Bt <- unlist(lapply(sampled_blocks, function(b) b$Bt))
BtB <- sum(sapply(sampled_blocks, function(b) b$BtB))
gammaj=gamma1j=gamma*runif(1,0.95,1.05)
deltaj=0*gamma1j
errorj=1
thetaj=theta*runif(1,0.95,1.05)
for(iterj in 1:sampling.iter){
theta_prevj=thetaj
indvalidj=which(gamma1j==0)
indvalidj=intersect(indj,indvalidj)
if(length(indvalidj)<(0.55*length(indj))) indvalidj=sample(indj,0.6*length(indj))
Hinv=1/(BtB-sum(bXse[indvalidj]^2)*Rxy[1,1]-sum(LDSC[indvalidj]*byseinv[indvalidj]^2)*Omega[1,1])
g=sum(Bt*(by[indj]-as.vector(LD[indj,]%*%gammaj)))-sum(bXse[indvalidj])*Rxy[2,1]-sum(LDSC[indvalidj]*byseinv[indvalidj]^2)*Omega[1,2]
thetaj=g*Hinv
resgammaj=as.vector(Thetarhoj%*%(by[indj]-bX[indj]*thetaj-deltaj[indj]+rho*gamma1j[indj]))
gammaj[indj]=resgammaj
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=(deltaj+rho*(gammaj-gamma1j))
gammaj=gammaj*(gamma1j!=0)
errorj=abs(thetaj-theta_prevj)
if(iterj>5&errorj<max.eps) break
}
ThetaList[j]=thetaj
}
theta.se=robust_sd(ThetaList)
}else{
ThetaList=NULL
if(sum(gamma!=0)==0){
e=tilde.y-tilde.X*theta
adjf=m/(m-1)
Theta_valid=solve(LD[indvalid,indvalid])
tilde.X=as.vector(chol(Theta_valid)%*%bX[indvalid])
BtB=sum(tilde.X^2)
h=(BtB-sum(bXse[indvalid]^2)*Rxy[1,1]-sum(LDSC[indvalid]*byseinv[indvalid]^2)*Omega[1,1])
E=-tilde.X*e[indvalid]+bXse[indvalid]*byse[indvalid]-bXse[indvalid]^2*theta+LDSC[indvalid]*byseinv[indvalid]^2*Omega[1,2]-LDSC[indvalid]*byseinv[indvalid]^2*Omega[1,1]*theta
vartheta=sum(E^2)/h^2*adjf
theta.se=sqrt(vartheta)
}else{
adjf=m/(length(indvalid)-1)
bZ=as.matrix(cbind(bX,LD[,which(gamma!=0)]))
H=matrixMultiply(t(bZ),as.matrix(Theta%*%bZ))
H[1,1]=H[1,1]-sum(bXse[indvalid]^2)*Rxy[1,1]-sum(LDSC[indvalid]*byseinv[indvalid]^2)*Omega[1,1]
e=res
Hinv=solve(H)
E=-as.matrix(Theta%*%bZ)*e
for(i in 1:length(indvalid)){
E[indvalid[i],1]=E[indvalid[i],1]+bXse[indvalid[i]]*Rxy[1,2]+(LDSC[indvalid[i]]*byseinv[indvalid[i]]^2)*Omega[1,2]-bXse[indvalid[i]]^2*Rxy[1,1]*theta-LDSC[indvalid[i]]*byseinv[indvalid[i]]^2*Omega[1,1]*theta
}
V=t(E)%*%E
covtheta=(Hinv%*%V%*%Hinv)*adjf
theta.se=sqrt(covtheta[1,1])
}
}

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini
A$theta.bootstrap=ThetaList
A$reliability.adjust=r
A$var_error=var_error
return(A)
}
