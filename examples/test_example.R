library(susieR)
library(CppMatrix)
library(varbvs)
library(mixtools)
library(Matrix)

devtools::load_all()
ARcov=function(p,rho){
  s=c(1:p)
  for(i in 1:p){
    s[i]=rho^(i-1)
  }
  return(toeplitz(s))
}

CScov=function(p,rho){
  return(matrix(rho,p,p)+(1-rho)*diag(p))
}
m=600 # number of IVs
p=10 # number of exposure
n1=n0=1e5 # outcome sample size
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
LD=kronecker(diag(150),ARcov(4,0.7)) # LD matrix
Theta=solve(LD)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.01,p),.00667)*15 # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=Bse=array(0,c(500,p,9))
Btime=matrix(0,500,9)
cluster.index=kronecker(c(1:150),rep(1,4))
theta0=c(1,-0.5,rep(0,6),-0.5,1)
UHP.var=1
UHP.frac=0.05*1
UHP.frac=0.01
CHP.frac=0.05*2
iter=1

A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,Nxy=Nxy,non.zero.frac=rep(0.8,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
Rxy=A$Rxy
par(mfrow=c(2,2))
plot(A$bX[-c(1:480),]%*%c(0,0,0,0,1,-1,rep(0,4)),A$by[-c(1:480)])
plot(A$bX[c(1:480),]%*%theta0,A$by[c(1:480)])
plot(A$bX0[-c(1:480),]%*%c(0,0,0,0,1,-1,rep(0,4)),A$by0[-c(1:480)])
plot(A$bX0[c(1:480),]%*%theta0,A$by0[c(1:480)])
Lvec=c(1:min(10,nrow(bX)));pip.thres=0.5;tauvec=seq(3,50,by=2);max.iter=100;max.eps=0.001;susie.iter=100;ebic.theta=1;ebic.gamma=2;reliability.thres=0.8;ADMM.rho=2;maxdiff=1.5;sampling.time=100;sampling.iter=10;theta.ini=F;gamma.ini=F
main.cluster.thres=0.49

t1=Sys.time()
fit.ipod=MRBEEX(Method="IPOD",use.susie=F,by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rnn*Ruv,cluster.index=cluster.index,reliability.thres=0.8,tauvec=c(2.5,3,3.5,4:10,seq(12,30,2)),ADMM.rho=2,sampling.time=200,maxdiff=3,sampling.iter=5)
t2=Sys.time()
ipod.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.susie=MRBEEX(Method="IPOD",use.susie=T,by,bX,byse,bXse,LD=LD,Lvec=c(1:6),pip.thres=0.2,Rxy=Ruv*Rnn,ADMM.rho=2,tauvec=c(2.5,3,3.5,4:10,seq(12,30,2)),ebic.gamma=1,cluster.index=cluster.index,sampling.time=200,maxdiff=3,reliability.thres=0.5,theta.ini=fit.ipod$theta,gamma.ini=fit.ipod$gamma,sampling.iter=5)
t2=Sys.time()
susie.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.mixture=MRBEEX(Method="Mixture",use.susie=F,by,bX,byse,bXse,LD=LD,Rxy=Ruv*Rnn,ebic.gamma=1,cluster.index=cluster.index,sampling.time=200,maxdiff=3,reliability.thres=0.8,pip.thres=0.15,sampling.iter=5)
t2=Sys.time()
mixture.time=difftime(t2, t1, units = "secs")
if(is.null(fit.mixture$IsMixture)==1){
  ThetaVec=cbind(fit.mixture$theta1,fit.mixture$theta2)
  ThetaVecSE=cbind(fit.mixture$theta.se1,fit.mixture$theta.se2)
  ThetaNorm=colMeans((ThetaVec-cbind(theta0,theta0))^2)
  mixture.theta=ThetaVec[,which.min(ThetaNorm)]
  mixture.theta.se=ThetaVecSE[,which.min(ThetaNorm)]
}else{
  fit.mixture=fit.ipod
  mixture.theta=fit.mixture$theta
  mixture.theta.se=fit.mixture$theta.se
}

t1=Sys.time()
fit.mixture.susie=MRBEEX(Method="Mixture",use.susie=T,by,bX,byse,bXse,LD=LD,Rxy=Ruv*Rnn,ebic.gamma=1,cluster.index=cluster.index,sampling.time=200,maxdiff=3,reliability.thres=0.8,sampling.iter=5)
t2=Sys.time()
mixture.susie.time=difftime(t2, t1, units = "secs")
if(is.null(fit.mixture.susie$IsMixture)==1){
  ThetaVec=cbind(fit.mixture.susie$theta1,fit.mixture.susie$theta2)
  ThetaVecSE=cbind(fit.mixture.susie$theta.se1,fit.mixture.susie$theta.se2)
  ThetaNorm=colMeans((ThetaVec-cbind(theta0,theta0))^2)
  mixture.susie.theta=ThetaVec[,which.min(ThetaNorm)]
  mixture.susie.theta.se=ThetaVecSE[,which.min(ThetaNorm)]
}else{
  fit.mixture.susie=fit.susie
  mixture.susie.theta=fit.mixture.susie$theta
  mixture.susie.theta.se=fit.mixture.susie$theta.se
}

cbind(theta0,fit.ipod$theta,fit.susie$theta,fit.mixture$theta1,fit.mixture.susie$theta1)
c(ipod.time,susie.time,mixture.time,mixture.susie.time)
par(mfrow=c(2,2))
barplot(fit.ipod$gamma/byse)
barplot(fit.susie$gamma/byse)
barplot(fit.mixture$gamma/byse)
barplot(fit.mixture.susie$gamma/byse)
cor(cbind(A$pleiotropy,fit.ipod$gamma,fit.susie$gamma,fit.mixture$gamma,fit.mixture.susie$gamma))
c(min(fit.ipod$Bic),min(fit.susie$Bic),min(fit.mixture$Bic),min(fit.mixture.susie$Bic))
