library(susieR)
library(CppMatrix)
library(varbvs)
library(mixtools)
library(Matrix)
library(devtools)
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
Bcov=function(B,Bse,theta0){
Cov=B
q=ncol(B)
for(i in 1:q){
b=B[,i]
se=Bse[,i]
diff=abs(b-theta0)/2
Cov[,i]=ifelse(diff<=se,1,0)
}
return(Cov)
}
m=1500 # number of IVs
p=10 # number of exposure
n1=n0=2e5 # outcome sample size
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
LD=kronecker(diag(125*3),ARcov(4,0.7)) # LD matrix
Theta=matrixInverse(LD)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.01,p),.00667)*15 # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=Bse=array(0,c(500,p,9))
Btime=matrix(0,500,9)
cluster.index=kronecker(c(1:(125*3)),rep(1,4))
theta0=c(1,-0.5,rep(0,6),-0.5,1)/2
UHP.var=1
UHP.frac=0.05*1
#UHP.frac=0.00
CHP.frac=0
iter=1

Btheta=Bse=BCov=array(0,c(100,p,3))
while(iter<=100){
A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,Nxy=Nxy,non.zero.frac=rep(0.8,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="ash",cluster.index=cluster.index)
bX=A$bX
#bX[,5]=bX[,5]-LD%*%A$bX0[,5]+LD%*%A$bX0[,1]
by=A$by
bXse=A$bXse
byse=A$byse
Rxy=A$Rxy
#par(mfrow=c(2,2))
#plot(A$bX[-c(1:480),]%*%c(0,0,0,0,1,-1,rep(0,4)),A$by[-c(1:480)])
#plot(A$bX[c(1:480),]%*%theta0,A$by[c(1:480)])
#plot(A$bX0[-c(1:480),]%*%c(0,0,0,0,1,-1,rep(0,4)),A$by0[-c(1:480)])
#plot(A$bX0[c(1:480),]%*%theta0,A$by0[c(1:480)])
Lvec=c(1:min(10,nrow(bX)));pip.thres=0.5;tauvec=seq(3,50,by=2);max.iter=100;max.eps=0.001;susie.iter=100;ebic.theta=1;ebic.gamma=2;reliability.thres=0.8;admm.rho=2;maxdiff=1.5;sampling.time=10;sampling.iter=5;theta.ini=F;gamma.ini=F
main.cluster.thres=0.45
ridge.diff=1e5;verbose=T;pip.min=0.1;cred.pip.thres=0.95;group.penalize=F;group.index=c(1:ncol(bX)[1]);group.diff=10;coverage.causal=0.95;LDSC=NULL;Omega=NULL;estimate_residual_variance=T;prob.shrinkage=0.5;estimate_residual_method="MoM";sampling.strategy="bootstrap"

t1=Sys.time()
fit.ipod=MRBEEX(method="IPOD",use.susie=F,by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rnn*Ruv,cluster.index=cluster.index,reliability.thres=0.8,tauvec=c(2.5,3,3.5,4,5,6),admm.rho=1,sampling.time=100,maxdiff=3,sampling.iter=8)
t2=Sys.time()
ipod.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.susie=MRBEEX(method="IPOD",use.susie=T,by,bX,byse,bXse,LD=LD,Lvec=c(1:6),pip.thres=0.2,Rxy=Ruv*Rnn,admm.rho=2,tauvec=c(2.5,3,3.5,4:10),ebic.gamma=1,cluster.index=cluster.index,sampling.time=100,maxdiff=3,reliability.thres=0.5,theta.ini=fit.ipod$theta,gamma.ini=fit.ipod$gamma,sampling.iter=10,coverage.causal=0.9)
t2=Sys.time()
susie.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.susie1=MRBEEX(method="IPOD",use.susie=T,by,bX,byse,bXse,LD=LD,Lvec=c(1:6),pip.thres=0.2,Rxy=Ruv*Rnn,admm.rho=2,tauvec=c(2.5,3,3.5,4:10),ebic.gamma=1,cluster.index=cluster.index,sampling.time=100,maxdiff=3,reliability.thres=0.5,theta.ini=fit.ipod$theta,gamma.ini=fit.ipod$gamma,sampling.iter=10,coverage.causal=0.9,sampling.strategy="subsampling")
t2=Sys.time()
susie.time=difftime(t2, t1, units = "secs")

# t1=Sys.time()
# fit.mixture=MRBEEX(method="Mixture",use.susie=F,by,bX,byse,bXse,LD=LD,Rxy=Ruv*Rnn,ebic.gamma=1,cluster.index=cluster.index,sampling.time=100,maxdiff=3,reliability.thres=0.8,pip.thres=0.15,sampling.iter=5)
# t2=Sys.time()
# mixture.time=difftime(t2, t1, units = "secs")
# ThetaVec=cbind(fit.mixture$theta1,fit.mixture$theta2)
# ThetaVecSE=cbind(fit.mixture$theta.se1,fit.mixture$theta.se2)
# ThetaNorm=colMeans((ThetaVec-cbind(theta0,theta0))^2)
# mixture.theta=ThetaVec[,which.min(ThetaNorm)]
# mixture.theta.se=ThetaVecSE[,which.min(ThetaNorm)]

# t1=Sys.time()
# fit.mixture.susie=MRBEEX(method="Mixture",use.susie=T,by,bX,byse,bXse,LD=LD,Rxy=Ruv*Rnn,ebic.gamma=1,cluster.index=cluster.index,sampling.time=100,maxdiff=3,reliability.thres=0.8,sampling.iter=10,sampling.strategy="subsampling")
# t2=Sys.time()
# mixture.susie.time=difftime(t2, t1, units = "secs")
# ThetaVec=cbind(fit.mixture.susie$theta1,fit.mixture.susie$theta2)
# ThetaVecSE=cbind(fit.mixture.susie$theta.se1,fit.mixture.susie$theta.se2)
# ThetaNorm=colMeans((ThetaVec-cbind(theta0,theta0))^2)
# mixture.susie.theta=ThetaVec[,which.min(ThetaNorm)]
# mixture.susie.theta.se=ThetaVecSE[,which.min(ThetaNorm)]
#
# t1=Sys.time()
# fit.mixture.susie1=MRBEEX(method="Mixture",use.susie=T,by,bX,byse,bXse,LD=LD,Rxy=Ruv*Rnn,ebic.gamma=1,cluster.index=cluster.index,sampling.time=100,maxdiff=3,reliability.thres=0.8,sampling.iter=10,sampling.strategy="subsampling")
# t2=Sys.time()
# mixture.susie.time=difftime(t2, t1, units = "secs")
# ThetaVec=cbind(fit.mixture.susie1$theta1,fit.mixture.susie1$theta2)
# ThetaVecSE=cbind(fit.mixture.susie1$theta.se1,fit.mixture.susie1$theta.se2)
# ThetaNorm=colMeans((ThetaVec-cbind(theta0,theta0))^2)
# mixture.susie.theta1=ThetaVec[,which.min(ThetaNorm)]
# mixture.susie.theta.se1=ThetaVecSE[,which.min(ThetaNorm)]

# cbind(theta0,fit.ipod$theta,fit.susie$theta,fit.mixture$theta1,fit.mixture.susie$theta1)
# print(c(min(fit.ipod$Bic),min(fit.susie$Bic),min(fit.mixture$Bic),min(fit.mixture.susie$Bic)))
# print(c(colMeans(fit.mixture$Voting$Weight),colMeans(fit.mixture.susie$Voting$Weight)))
# Btheta[iter,,]=cbind(fit.ipod$theta,fit.susie$theta,fit.mixture$theta1,fit.mixture.susie$theta1)
Btheta[iter,,]=cbind(fit.ipod$theta,fit.susie$theta,fit.susie1$theta)
Bse[iter,,]=cbind(fit.ipod$theta.se,fit.susie$theta.se,fit.susie1$theta.se)
BCov[iter,,]=Bcov(Btheta[iter,,],Bse[iter,,],theta0)
if(iter %% 10 ==0){
par(mfrow=c(2,2))
boxplot(Btheta[1:iter,1,])
lines(c(0:4),rep(theta0[1],5))
boxplot(Btheta[1:iter,2,])
lines(c(0:4),rep(theta0[2],5))
boxplot(Btheta[1:iter,9,])
lines(c(0:4),rep(theta0[9],5))
boxplot(Btheta[1:iter,10,])
lines(c(0:4),rep(theta0[10],5))
}
iter=iter+1
}
