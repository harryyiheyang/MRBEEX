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

coverage=function(theta,theta.se,theta0){
  delta=abs(theta-theta0)/2
  cover=ifelse(delta<=theta.se,1,0)
  return(cover)
}
m=200 # number of IVs
p=10 # number of exposure
n1=300
n0=1e6 # outcome sample size
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
LD=kronecker(CScov(4,0.7),ARcov(50,0.7)) # LD matrix
Theta=solve(LD)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.3,p),.000667*15) # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=matrix(0,100,p)
Bse=Bcov=array(0,c(100,p,3))
cluster.index=kronecker(c(1:50),rep(1,4))
theta0=c(1,rep(0,8),0)
UHP.var=1
UHP.frac=0.05*1
#UHP.frac=0.01*0
CHP.frac=0*2

A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,Nxy=Nxy,non.zero.frac=rep(0.01,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX
by=A$by
bXse=A$bXse
byse=A$byse
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.8;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.5;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
xQTLfitList=NULL
xQTL.cred.thres=0.95
xQTL.pip.min=0.2
outlier.switch=F
Annotation=cbind(1,rowSums(A$bX0!=0))
Annotation=NULL
output.labels=NULL;
carma.iter=5;carma.inner.iter=5;xQTL.max.num=10;
carma.epsilon.threshold=1e-3;

reliability.thres=0.75;Lvec=c(1:5);causal.pip.thres=0.2;
xQTL.selection.rule="top_K";
top_K=1;xQTL.pip.min=0.2;
xQTL.max.L=10;xQTL.cred.thres=0.95;xQTL.pip.thres=0.5;
xQTL.Nvec=Nvec=rep(n1,p);tauvec=seq(3,30,by=3);xQTL.weight=NULL;
admm.rho=2;ridge.diff=1e3;
max.iter=100;max.eps=0.001;susie.iter=500;
ebic.theta=0;ebic.gamma=1;
theta.ini=F;gamma.ini=F;xQTLfitList=NULL;
sampling.iter=10;sampling.time=1000;sampling.size=0.5;
batch.size=1;verbose=T

fit1=CisMRBEE_UV(by=by,bX=bX[,1],byse=byse,bXse=bXse[,1],LD=LD,Rxy=Rxy[c(1,11),c(1,11)],xQTL.N=n1)


