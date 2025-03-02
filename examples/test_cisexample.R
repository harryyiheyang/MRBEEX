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
n1=30000
n0=1e6 # outcome sample size
Rbb=ARcov(p,-0.5) # exposure covariance
Ruv=ARcov(p+1,-0.3) # estimation error covariance
LD=kronecker(CScov(4,0.7),ARcov(50,0.7)) # LD matrix
Theta=solve(LD)
Nxy=c(rep(n1,p),n0) # sample size vector
Hxy=c(rep(.3,p),.000667*15) # H2 vector
Rnn=CScov(p=p+1,1)
Btheta=Bse=array(0,c(100,p,4))
cluster.index=kronecker(c(1:50),rep(1,4))
theta0=c(1,rep(0,9))
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
estimate_residual_variance=T;residual_variance=1;
reliability.thres=0.9;Lvec=c(1:5);pip.thres=0.2;
xQTL.max.L=10;xQTL.sampling=1000;
xQTL.pip.thres=0.5;xQTL.Nvec=rep(n1,p);
block.rho=0;robust.sandwith=T;
tauvec=seq(3,30,by=3);rho=2;ridge=0.05;
max.iter=100;max.eps=0.001;susie.iter=100;
ebic.theta=1;ebic.gamma=2;maxdiff=3;
theta.ini=F;gamma.ini=F
Rxy=A$Rxy
eQTLfitList=NULL
xQTL.cred.thres=0.95
xQTL.pip.min=0.2
outlier.switch=F
Annotation=cbind(1,rowSums(A$bX0!=0))
Annotation=NULL
output.labels=NULL;
carma.iter=5;carma.inner.iter=5;xQTL.max.num=10;
carma.epsilon.threshold=1e-3;

fit1=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,causal.pip.thres=0.1)
fit1$theta
#fit2=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,causal.pip.thres=0.1,eQTL.method="CARMA")
#fit2$theta

bX[,2]=LD%*%A$bX0[,1]+(A$bX[,2]-LD%*%A$bX0[,2])
bX[,2]=-bX[,2]
Rxy[-2,2]=-Rxy[-2,2];Rxy[2,-2]=-Rxy[2,-2];
bX[,3]=LD%*%A$bX0[,1]+(A$bX[,3]-LD%*%A$bX0[,3])

fit3=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,causal.pip.thres=0.1,top_K=3)
fit3$theta
#fit4=CisMRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,xQTL.Nvec=rep(n1,p),ridge=100,causal.pip.thres=0.1,eQTL.method="CARMA")
#fit4$theta

Btheta[iter,,]=cbind(fit1$theta,fit2$theta,fit3$theta,fit4$theta)#,fit5$theta,fit6$theta,fit7$theta,fit8$theta)
Bse[iter,,]=cbind(fit1$theta.se,fit2$theta.se,fit3$theta.se,fit4$theta.se)#,fit5$theta.se,fit6$theta.se,fit7$theta.se,fit8$theta.se)
iter=iter+1
if(iter %% 25==0){print(iter)}
}

cbind(colSD(Btheta[1:100,1,]),colMeans(Bse[1:100,1,]))
