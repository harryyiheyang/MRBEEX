MRBEE_IPOD_UV=function(by,bX,byse,bXse,LD=LD,Rxy,cluster.index,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,ebic.gamma=1,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=5,theta.ini=F,gamma.ini=F,reliability.thres=0.8,LDSC=NULL,Omega=NULL){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=length(bX)
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
isLD=F
LD=Theta=TC=Matrix(diag(m),sparse=T)
RC=diag(m)
byinv=by
bXinv=bX
tilde.y=by
tilde.X=bX
Bt=t(bX)
BtB=sum(bX*bX)
Thetarho=Matrix(diag(m)/(1+rho),sparse=T)
}
r=reliability.adj.uv(bX,bXse,Theta=Theta,thres=reliability.thres)
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
if(isLD){
gamma=as.vector(Thetarho%*%(by-bX*theta-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
}else{
gamma=as.vector((by-bX*theta-delta+rho*gamma1))/(rho+1)
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
}
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

ThetaList=c(1:sampling.time)
for(j in 1:sampling.time){
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
LDj=LD[indj,indj]
Thetaj=solve(LDj)
Thetarhoj=solve(LDj+rho*diag(length(indj)))
Bt=as.vector(Thetaj%*%bX[indj])
BtB=sum(bX[indj]*(Thetaj%*%bX[indj]))
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
if(isLD){
resgammaj=as.vector(Thetarhoj%*%(by[indj]-bX[indj]*thetaj-deltaj[indj]+rho*gamma1j[indj]))
gammaj[indj]=resgammaj
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=(deltaj+rho*(gammaj-gamma1j))
}else{
resgammaj=as.vector((by[indj]-bX[indj]*thetaj-deltaj+rho*gamma1j))/(1+rho)
gammaj[indj]=resgammaj
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=(deltaj+rho*(gammaj-gamma1j))
}
gammaj=gammaj*(gamma1j!=0)
errorj=abs(thetaj-theta_prevj)
if(iterj>3&errorj<max.eps) break
}
ThetaList[j]=thetaj
}
theta.se=mad(ThetaList)


A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini
A$theta.bootstrap=ThetaList
A$reliability.adjust=r
A$theta.pratt=getPratt.uv(bX=bX,by=by,bXse=bXse,byse=byse,Theta=Theta,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma,Theta=Theta,LD=LD)
A$var_error=var_error
return(A)
}
