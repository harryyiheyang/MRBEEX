library(CppMatrix)
library(Matrix)
library(devtools)
devtools::load_all()

ARcov=function(p,rho){
  s=c(1:p)
  for(i in 1:p) s[i]=rho^(i-1)
  return(toeplitz(s))
}

CScov=function(p,rho){
  return(matrix(rho,p,p)+(1-rho)*diag(p))
}

args=commandArgs(trailingOnly=TRUE)
n.sim=if(length(args)>=1) as.integer(args[1]) else 500L
result.file=if(length(args)>=2) args[2] else file.path("examples","diagnose_ipod_pleiotropy_results.rds")
summary.file=if(length(args)>=3) args[3] else file.path("examples","diagnose_ipod_pleiotropy_summary.csv")
t0=Sys.time()

block.num=50
m=block.num*4
p=10
n1=n0=2e5
Rbb=ARcov(p,-0.5)
Ruv=ARcov(p+1,-0.3)
LD=kronecker(diag(block.num),ARcov(4,0.7))
Nxy=c(rep(n1,p),n0)
Hxy=c(rep(.01,p),.00667)*15
Rnn=CScov(p=p+1,1)
cluster.index=kronecker(c(1:block.num),rep(1,4))
theta0=c(1,-0.5,rep(0,6),-0.5,1)/2
Rxy=Ruv*Rnn

Btheta0=Bse0=matrix(0,n.sim,p)
Btheta5=Bse5=matrix(0,n.sim,p)
Bselect0=Bselect5=matrix(FALSE,n.sim,m)
Btruth5=matrix(FALSE,n.sim,m)

set.seed(2025)
for(iter in 1:n.sim){
  A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,
    Nxy=Nxy,non.zero.frac=rep(0.8,p),UHP.frac=0,CHP.frac=0,
    UHP.var=1,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,
    UHP.dis="ash",cluster.index=cluster.index)
  fit=MRBEE_LDA(use.susie=F,by=A$by,bX=A$bX,byse=A$byse,bXse=A$bXse,
    LD=LD,Rxy=Rxy,cluster.index=cluster.index,reliability.thres=0.6,
    tauvec=c(2.5,3,3.5,4,5,6),admm.rho=1,maxdiff=3,verbose=F)
  Btheta0[iter,]=fit$theta
  Bse0[iter,]=fit$theta.se
  Bselect0[iter,]=fit$gamma!=0
  if(iter%%50==0) cat("UHP.frac=0:",iter,"/",n.sim,"\n")
}

set.seed(2026)
for(iter in 1:n.sim){
  A=MRBEEX::summary_generation(theta=theta0,m=m,Rbb=Rbb,Ruv=Ruv,Rnn=Rnn,LD=LD,
    Nxy=Nxy,non.zero.frac=rep(0.8,p),UHP.frac=0.05,CHP.frac=0,
    UHP.var=1,CHP.effect=c(0,0,0,0,1,-1,rep(0,4)),Hxy=Hxy,
    UHP.dis="ash",cluster.index=cluster.index)
  fit=MRBEE_LDA(use.susie=F,by=A$by,bX=A$bX,byse=A$byse,bXse=A$bXse,
    LD=LD,Rxy=Rxy,cluster.index=cluster.index,reliability.thres=0.6,
    tauvec=c(2.5,3,3.5,4,5,6),admm.rho=1,maxdiff=3,verbose=F)
  Btheta5[iter,]=fit$theta
  Bse5[iter,]=fit$theta.se
  Bselect5[iter,]=fit$gamma!=0
  Btruth5[iter,]=A$pleiotropy!=0
  if(iter%%50==0) cat("UHP.frac=0.05:",iter,"/",n.sim,"\n")
}

S=data.frame()
for(j in 1:p){
  esd=sd(Btheta0[,j])
  S=rbind(S,data.frame(UHP.frac=0,exposure=paste0("X",j),truth=theta0[j],
    mean_theta=mean(Btheta0[,j]),bias=mean(Btheta0[,j])-theta0[j],empirical_sd=esd,
    mean_se=mean(Bse0[,j]),se_sd_ratio=mean(Bse0[,j])/esd,
    coverage=mean(abs(Btheta0[,j]-theta0[j])<=1.96*Bse0[,j])))
}
for(j in 1:p){
  esd=sd(Btheta5[,j])
  S=rbind(S,data.frame(UHP.frac=0.05,exposure=paste0("X",j),truth=theta0[j],
    mean_theta=mean(Btheta5[,j]),bias=mean(Btheta5[,j])-theta0[j],empirical_sd=esd,
    mean_se=mean(Bse5[,j]),se_sd_ratio=mean(Bse5[,j])/esd,
    coverage=mean(abs(Btheta5[,j]-theta0[j])<=1.96*Bse5[,j])))
}

causal=which(theta0!=0)
cat("\nCausal exposure averages\n")
print(aggregate(cbind(empirical_sd,mean_se,se_sd_ratio,coverage)~UHP.frac,
  data=S[S$exposure%in%paste0("X",causal),],mean))
cat("\nPleiotropy selection\n")
cat("False positive SNPs per fit, UHP.frac=0:",mean(rowSums(Bselect0)),"\n")
cat("Sensitivity, UHP.frac=0.05:",sum(Bselect5&Btruth5)/sum(Btruth5),"\n")
cat("False discovery proportion, UHP.frac=0.05:",sum(Bselect5&!Btruth5)/sum(Bselect5),"\n")
cat("Selected SNPs per fit, UHP.frac=0.05:",mean(rowSums(Bselect5)),"\n")

elapsed=as.numeric(difftime(Sys.time(),t0,units="secs"))
saveRDS(list(theta0=theta0,Btheta0=Btheta0,Bse0=Bse0,Bselect0=Bselect0,
  Btheta5=Btheta5,Bse5=Bse5,Bselect5=Bselect5,Btruth5=Btruth5,S=S,
  n.sim=n.sim,elapsed=elapsed),result.file)
write.csv(S,summary.file,row.names=FALSE)
cat("Elapsed seconds:",round(elapsed,1),"\n")
