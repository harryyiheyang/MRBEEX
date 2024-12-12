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
LD1=kronecker(diag(100),ARcov(2,0.9)) # LD matrix
Theta1=solve(LD1)
LD2=kronecker(diag(100),ARcov(2,0.7)) # LD matrix
Theta2=solve(LD2)
LD3=kronecker(diag(100),ARcov(2,0.5)) # LD matrix
Theta3=solve(LD3)
LD4=kronecker(diag(100),ARcov(2,0.3)) # LD matrix
Theta4=solve(LD4)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.3,p),.000667*15) # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=Bse=array(0,c(100,p,4))
cluster.index=kronecker(c(1:50),rep(1,4))
theta0=c(1,-0.5,rep(0,6),-0.5,1)
UHP.var=1
UHP.frac=0.05*1
UHP.frac=0.01*0
CHP.frac=0*2

iter=1
while(iter<=100){
A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD1,Nxy=Nxy,non.zero.frac=rep(0.02,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.05;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
eQTLfitList=NULL
xQTL.cred.thres=0.95

fit1=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD1,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,pip.thres=0.1)

A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD2,Nxy=Nxy,non.zero.frac=rep(0.02,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.05;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
eQTLfitList=NULL
xQTL.cred.thres=0.95

fit2=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD2,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,pip.thres=0.1)

A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD3,Nxy=Nxy,non.zero.frac=rep(0.02,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.05;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
eQTLfitList=NULL
xQTL.cred.thres=0.95

fit3=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD3,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,pip.thres=0.1)

A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD4,Nxy=Nxy,non.zero.frac=rep(0.02,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.05;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
eQTLfitList=NULL
xQTL.cred.thres=0.95

fit4=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD4,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,pip.thres=0.1)

Btheta[iter,,]=cbind(fit1$theta,fit2$theta,fit3$theta,fit4$theta)#,fit5$theta,fit6$theta,fit7$theta,fit8$theta)
Bse[iter,,]=cbind(fit1$theta.se,fit2$theta.se,fit3$theta.se,fit4$theta.se)#,fit5$theta.se,fit6$theta.se,fit7$theta.se,fit8$theta.se)
iter=iter+1
if(iter %% 25==0){print(iter)}
}

cbind(colSD(Btheta[1:100,1,]),colMeans(Bse[1:100,1,]))
