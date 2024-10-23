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
m=200 # number of IVs
p=10 # number of exposure
n1=3000
n0=1e6 # outcome sample size
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
LD=kronecker(diag(50),ARcov(4,0.7)) # LD matrix
Theta=solve(LD)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.3,p),.000667*15) # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=Bse=array(0,c(100,p,2))
cluster.index=kronecker(c(1:50),rep(1,4))
theta0=c(1,-0.5,rep(0,6),-0.5,1)
UHP.var=1
UHP.frac=0.05*1
UHP.frac=0.01*0
CHP.frac=0*2

iter=1
while(iter<=100){
A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,Nxy=Nxy,non.zero.frac=rep(0.02,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
cluster.index=4;
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.2;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy

fit1=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),ridge = 0)
fit2=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=F,xQTL.Nvec=rep(n1,p),eQTLfitList=fit1$eQTLfitList,ridge = 0)
#fit3=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),cluster.index=5,block.rho=0.25,eQTLfitList=fit1$eQTLfitList)
#fit4=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),cluster.index=4,block.rho=0.5,eQTLfitList=fit1$eQTLfitList)
#fit5=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),cluster.index=2,block.rho=0.25,eQTLfitList=fit1$eQTLfitList)
#fit6=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),cluster.index=1,block.rho=0.25,eQTLfitList=fit1$eQTLfitList)
#fit7=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),cluster.index=10,block.rho=0.25,eQTLfitList=fit1$eQTLfitList)
#fit8=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,use.susie=T,xQTL.Nvec=rep(n1,p),robust.sandwith=F,eQTLfitList=fit1$eQTLfitList)

Btheta[iter,,]=cbind(fit1$theta,fit2$theta)#,fit3$theta,fit4$theta,fit5$theta,fit6$theta,fit7$theta,fit8$theta)
Bse[iter,,]=cbind(fit1$theta.se,fit2$theta.se)#,fit3$theta.se,fit4$theta.se,fit5$theta.se,fit6$theta.se,fit7$theta.se,fit8$theta.se)
cbind(fit1$theta.se,fit2$theta.se)
iter=iter+1
if(iter %% 25==0){print(iter)}
}

cbind(colSD(Btheta[1:100,1,]),colMeans(Bse[1:100,1,]))
