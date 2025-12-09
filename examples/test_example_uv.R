library(susieR)
library(CppMatrix)
library(varbvs)
library(mixtools)
library(Matrix)
library(devtools)
devtools::document()
is_coverage=function(theta,theta0,theta.se){
d=theta-theta0
return(ifelse(abs(d/2)<=theta.se,1,0))
}
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
m=1500 # number of IVs
p=10 # number of exposure
n1=n0=1e5 # outcome sample size
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
theta0=c(1,rep(0,9))
UHP.var=1
UHP.frac=0.05*1
#UHP.frac=0.00
CHP.frac=0.05*0
iter=1

Btheta=Bse=Bcov=matrix(0,100,4)
while(iter<101){
A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,Nxy=Nxy,non.zero.frac=rep(0.8,p),UHP.frac=UHP.frac,CHP.frac=CHP.frac,UHP.var=UHP.var,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,UHP.dis="normal",cluster.index=cluster.index)
bX=A$bX[,1]
by=A$by
bXse=A$bXse[,1]
byse=A$byse
Rxy=A$Rxy[c(1,11),c(1,11)]
Lvec=c(1:min(10,nrow(bX)));pip.thres=0.5;tauvec=seq(3,50,by=2);max.iter=100;max.eps=0.001;susie.iter=100;ebic.theta=1;ebic.gamma=2;reliability.thres=0.8;admm.rho=2;maxdiff=1.5;sampling.time=10;sampling.iter=5;theta.ini=F;gamma.ini=F
main.cluster.thres=0.45

t1=Sys.time()
fit.ipod=MRBEEX_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,sampling.time=100)
t2=Sys.time()
ipod.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.greedy=MRBEEX_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Method="Greedy",sampling.time=100)
t2=Sys.time()
greedy.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.ipod1=MRBEEX_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,sampling.time=0)
t2=Sys.time()
ipod.time1=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit.greedy1=MRBEEX_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,Method="Greedy",sampling.time=0)
t2=Sys.time()
greedy.time1=difftime(t2, t1, units = "secs")

Btheta[iter,]=cbind(fit.ipod$theta,fit.greedy$theta,fit.ipod1$theta1,fit.greedy1$theta1)
Bse[iter,]=cbind(fit.ipod$theta.se,fit.greedy$theta.se,fit.ipod1$theta.se,fit.greedy1$theta.se)
Bcov[iter,]=cbind(is_coverage(fit.ipod$theta,theta0[1],fit.ipod$theta.se),is_coverage(fit.greedy$theta,theta0[1],fit.greedy$theta.se),is_coverage(fit.ipod1$theta,theta0[1],fit.ipod1$theta.se),is_coverage(fit.greedy1$theta,theta0[1],fit.greedy1$theta.se))
if(iter %% 10 ==0){
colMeans(Bcov[1:iter,])
}
iter=iter+1
}
