library(data.table)
library(dplyr)
library(devtools)
load_all()
options(bitmapType="cairo")
EUR_Clumped <- readRDS("data/EUR_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-12]
bXse=SE[,-12]
by=BETA[,12]
byse=SE[,12]
NAM=c(colnames(bX),"Ischemic Stroke")
fit_EUR=MRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy[NAM,NAM],ebic.theta=0,LD="identity",use.susie=F)
indEUR=which(abs(fit_EUR$theta/fit_EUR$theta.se)>2.807034)
theta.source=fit_EUR$theta;
theta.source.cov=fit_EUR$theta.cov;
theta.source[-indEUR]=0
theta.source.cov[-indEUR,-indEUR]=0

EAS_Clumped <- readRDS("data/EAS_Clumped.rds") %>% list2env(.,envir=.GlobalEnv)
BETA=Z/sqrt(N)
SE=1/sqrt(N)
bX=BETA[,-12]
bXse=SE[,-12]
by=BETA[,12]
byse=SE[,12]
NAM=c(colnames(bX),"Ischemic Stroke")
fit_EAS=MRBEEX(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy[NAM,NAM],ebic.theta=0,LD="identity",use.susie=F)
fit_EAS$theta/fit_EAS$theta.se

cbind(fit_EUR$theta,fit_EAS$theta)
plot(fit_EUR$theta,fit_EAS$theta)
c(cor(fit_EUR$theta,fit_EAS$theta),cor(fit_EUR$theta,fit_EAS$theta,method="kendall"),cor(fit_EUR$theta,fit_EAS$theta,method="spearman"))

tauvec=seq(3,30,3);Lvec=c(1:6);admm.rho=3;ebic.delta=1;ebic.gamma=2;transfer.coef=1;susie.iter=200;pip.thres=0.3;max.iter=50;max.eps=1e-4;reliability.thres=0.8;ridge.diff=100;sampling.time=100;sampling.iter=10

fit_EAS_Tran=MRBEE_TL(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy[NAM,NAM],LD="identity",theta.source=theta.source,theta.source.cov=theta.source.cov)
Estimate=cbind(theta.source,fit_EAS_Tran$theta,fit_EAS_Tran$delta,fit_EAS$theta)
SE=cbind(sqrt(diag(theta.source.cov)),fit_EAS_Tran$theta.se,fit_EAS_Tran$delta.se,fit_EAS$theta.se)
Z=Estimate/SE
colnames(Z)=colnames(Estimate)=c("EUR","EAS_Tran","EAS_Diff","EAS")
print(Z)
print(Estimate)
