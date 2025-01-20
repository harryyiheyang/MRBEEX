library(data.table)
library(dplyr)
library(MRBEEX)
library(devtools)
load_all()
options(bitmapType="cairo")
EUR_Clumped <- readRDS("data/EUR_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bx=BETA[,1]
bxse=SE[,1]
by=BETA[,10]
byse=SE[,10]
NAME=colnames(BETA)
NAM=NAME[c(1,10)]
fit_EUR=MRBEEX_UV(by=by,bX=bx,byse=byse,bXse=bxse,Rxy=Rxy[NAM,NAM])
fit_EUR1=MRBEE.IMRP.UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy[NAM,NAM])

EAS_Clumped <- readRDS("data/EAS_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bx=BETA[,1]
bxse=SE[,1]
by=BETA[,10]
byse=SE[,10]
NAME=colnames(BETA)
NAM=NAME[c(1,10)]
fit_EAS=MRBEEX_UV(by=by,bX=bx,byse=byse,bXse=bxse,Rxy=Rxy[NAM,NAM])
fit_EAS1=MRBEE.IMRP.UV(by=by,bx=bx,byse=byse,bxse=bxse,Rxy=Rxy[NAM,NAM])

by=by;bX=bx;byse=byse;bXse=bxse;Rxy=Rxy[NAM,NAM];LD="identity";
cluster.index=c(1:length(by));theta.source=0.82
theta.source.cov=0.05^2
tauvec=seq(3,30,3);
admm.rho=3;ebic.delta=1;ebic.gamma=2;
transfer.coef=1;susie.iter=200;pip.thres=0.3;max.iter=50;max.eps=1e-4;
reliability.thres=0.8;ridge.diff=100;sampling.time=100;sampling.iter=10
fit_EAS_Tran=MRBEE_TL_UV(by=by,bX=bx,byse=byse,bXse=bxse,Rxy=Rxy[NAM,NAM],theta.source=fit_EUR$theta,theta.source.cov=fit_EUR$theta.se^2,LD="identity")
#c(fit_EAS_Tran$theta,fit_EAS$theta,fit_EUR$theta)
#par(mfrow=c(2,1))
#barplot(fit_EAS$gamma)
#barplot(fit_EAS_Tran$gamma)
