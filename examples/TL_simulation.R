ARcov=function(p,rho){
  s=c(1:p)
  for(i in 1:p){
    s[i]=rho^(i-1)
  }
  return(toeplitz(s))
}

CScov=function(p,rho){
  A=matrix(rho,p,p)+(1-rho)*diag(p)
  return(A)
}

library(data.table)
library(dplyr)
library(MASS)
library(MRBEEX)
library(devtools)
library(MendelianRandomization)

document()
theta_eur=c(1,0,0.5,0,0.2,0,0.2,rep(0,3))/2
theta_eas=c(1,0,0.5,0,0.2,0,0.2,rep(0,3))/2
theta_eas=c(1,0,0.3,0,0.2,0,0,rep(0,3))/2
n_eas=0.5e5
n_eur=5e5
m_eas=2000
m_eur=2000
h_eur=0.1
h_eas=0.1
SB_eur=kronecker(CScov(2,0.5),ARcov(5,0.5))
SB_eas=kronecker(CScov(2,0.6),ARcov(5,0.4))
SE_eur=SE_eas=CScov(11,0.3)
Rnn=CScov(11,0.95)
UHP.frac=0.1
UHP.var=1.5
iter=1

A_eur=summary_generation(theta=theta_eur,m=m_eur,Rbb=SB_eur,Ruv=SE_eur,Rnn=Rnn,LD="identity",Nxy=rep(n_eur,11),non.zero.frac=rep(0.8,10),UHP.frac=UHP.frac,CHP.frac=0,UHP.var=UHP.var,Hxy=rep(h_eur,11),UHP.dis="ash")
A_eas=summary_generation(theta=theta_eas,m=m_eas,Rbb=SB_eas,Ruv=SE_eas,Rnn=Rnn,LD="identity",Nxy=rep(n_eas,11),non.zero.frac=rep(0.8,10),UHP.frac=UHP.frac,CHP.frac=0,UHP.var=UHP.var,Hxy=rep(h_eas,11),UHP.dis="ash")
bX_eur=A_eur$bX
by_eur=A_eur$by
bXse_eur=A_eur$bXse
byse_eur=A_eur$byse
Rxy_eur=A_eur$Rxy
bX_eas=A_eas$bX
by_eas=A_eas$by
bXse_eas=A_eas$bXse
byse_eas=A_eas$byse
Rxy_eas=A_eas$Rxy

fit_eur_mrbee=MRBEEX::MRBEEX(by=by_eur,bX=bX_eur,byse=byse_eur,bXse=bXse_eur,Rxy=Rxy_eur,LD="identity")

t1=Sys.time()
fit_mrbee=MRBEE_IMRP(by=by_eas,bX=bX_eas,byse=byse_eas,bXse=bXse_eas,Rxy=Rxy_eas)
t2=Sys.time()
imrp.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_mrbee_susie=MRBEEX::MRBEEX(by=by_eas,bX=bX_eas,byse=byse_eas,bXse=bXse_eas,Rxy=Rxy_eas,LD="identity")
t2=Sys.time()
mrbee.susie.time=difftime(t2, t1, units = "secs")

MVINPUT=mr_mvinput(by=by_eas,bx=bX_eas,byse=byse_eas,bxse=bXse_eas)

# t1=Sys.time()
# fit_median=mr_mvmedian(MVINPUT,iterations=1000)
# t2=Sys.time()
# median.time=difftime(t2, t1, units = "secs")
#
# t1=Sys.time()
# fit_lasso=mr_mvlasso(MVINPUT)
# t2=Sys.time()
# lasso.time=difftime(t2, t1, units = "secs")
#
# t1=Sys.time()
# fit_mvcML=mr_mvcML(MVINPUT,n=n_eas,DP=F,rho_mat=Rxy_eas,K_vec=c(0:sum(A_eas$pleiotropy!=0)))
# t2=Sys.time()
# mvcML.time=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_tran_mrbee=MRBEE_TL(by=by_eas,bX=bX_eas,byse=byse_eas,bXse=bXse_eas,Rxy=Rxy_eas,theta.source=fit_eur_mrbee$theta,theta.source.cov=fit_eur_mrbee$theta.cov)
t2=Sys.time()
mrbee.tran.time=difftime(t2, t1, units = "secs")

Estimate=cbind(fit_tran_mrbee$theta,fit_mrbee_susie$theta,fit_mrbee$theta,fit_median@Estimate,fit_lasso@Estimate)
SE=cbind(fit_tran_mrbee$theta.se,fit_mrbee_susie$theta.se,fit_mrbee$theta.se,fit_median@StdError,fit_lasso@StdError)
time=c(mrbee.tran.time,mrbee.susie.time,imrp.time,mvcML.time,median.time,lasso.time)
print(apply(Estimate-theta_eas,2,norm,"2"))
Estimate
SE
time
#by=by_eas;bX=bX_eas;byse=byse_eas;bXse=bXse_eas;Rxy=Rxy_eas;Lvec=c(1:3);theta.source=fit_eur_mrbee$theta;theta.source.cov=fit_eur_mrbee$theta.cov
#susie.iter=500;pip.thres=0.2;max.iter=100;max.eps=1e-4;pv.thres=0.05;var.est="variance";FDR=T;adjust.method="Sidak";reliability.thres=0.8;ridge.diff=100
#Lvec=c(1:8);transfer.coef=1;ebic.theta=ebic.gamma=1;tauvec=seq(5,50,5);admm.rho=5;ebic.delta=2
#sampling.time=100;sampling.iter=10
