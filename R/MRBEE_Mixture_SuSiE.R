MRBEE_Mixture_SuSiE=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),main.cluster.thres=0.45,min.cluster.size=5,Lvec=c(1:min(5,ncol(bX))),pip.thres=0.2,ebic.theta=1,reliability.thres=0.8,sampling.time=100,max.iter=30,max.eps=5e-4,sampling.iter=5,susie.iter=100,ridge.diff=1e5,verbose=T,pip.min=0.1,cred.pip.thres=0.95,estimate_residual_variance=T,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,coverage.causal=0.95,LDSC=NULL,Omega=NULL,prob_shrinkage_coef=0.5,prob_shrinkage_size=4,estimate_residual_method="MoM",sampling.strategy="bootstrap",standardize=T){
########################### Basic information #######################
t1=Sys.time()
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
vary=Rxy[p+1,p+1]
min.cluster.size=5
cluster.index <- as.integer(factor(cluster.index))
if(LD[1]!="identity"){
isLD=T
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
}else{
isLD=F
LD=Theta=TC=diag(m)
tilde.y=by
tilde.X=bX
}
r=reliability.adj(bX,bXse/sqrt(diag(Rxy[1:p,1:p])),Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=LDSC,Omega=Omega)
############################ Initial Estimate #######################
fit.init=fit.mixture=regmixEM(y=tilde.y,x=tilde.X,k=2,epsilon=max.eps,maxit=300)
max.cluster=ifelse(sum(fit.init$posterior[,1]>main.cluster.thres)>(m/2),1,2)
cluster.ini.1=which(fit.init$posterior[,max.cluster]>main.cluster.thres)
cluster.ini.2=setdiff(1:m,cluster.ini.1)
if(length(cluster.ini.2)<min.cluster.size){
cluster.ini.2=c(1:min.cluster.size)
}
m1=length(cluster.ini.1)
m2=length(cluster.ini.2)
theta.ini.1=theta.ini.11=fit.init$beta[,max.cluster]
theta.ini.2=theta.ini.22=fit.init$beta[,setdiff(1:2,max.cluster)]
fit.susie.init1=susie(X=tilde.X[cluster.ini.1,],y=tilde.y[cluster.ini.1],L=5,intercept=F)
fit.susie.init2=susie(X=tilde.X[cluster.ini.2,],y=tilde.y[cluster.ini.2],L=5,intercept=F)
theta.ini.1=coef.susie(fit.susie.init1)[-1]*(fit.susie.init1$pip>0.1)
theta.ini.2=coef.susie(fit.susie.init2)[-1]*(fit.susie.init2$pip>0.1)
sigma.ini.1=sum((tilde.y[cluster.ini.1]-tilde.X[cluster.ini.1,]%*%theta.ini.1)^2)/(length(cluster.ini.1)-sum(theta.ini.1!=0))
sigma.ini.2=sum((tilde.y[cluster.ini.2]-tilde.X[cluster.ini.2,]%*%theta.ini.2)^2)/(length(cluster.ini.2)-sum(theta.ini.2!=0))
sigma.ini.2=max(0.25,sigma.ini.2)
sigma.ini.1=max(0.25,sigma.ini.1)
Voting.ini=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta.ini.1,theta2=theta.ini.2,sigma1=sigma.ini.1,sigma2=sigma.ini.2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
cluster.ini.2=which(Voting.ini$Cluster[,2]==1)
cluster.ini.1=which(Voting.ini$Cluster[,1]==1)
cluster.ratio.ini=c(length(cluster.ini.1),length(cluster.ini.2))/m
m1=length(cluster.ini.1)
m2=length(cluster.ini.2)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Initialization ends: ",time_to_print," secs\n"))
}
############################## Tuning Parameter ######################
t1=Sys.time()
q=length(Lvec)
Btheta1=array(0,c(p,q,q))
Btheta2=array(0,c(p,q,q))
Bbic=SIG1=SIG2=M1=M2=matrix(1e6,q,q)
for(v in 1:length(Lvec)){
fit.susie1=NULL
for(l in 1:v){
theta1=theta.ini.1
theta2=theta.ini.2
cluster1=cluster.ini.1
cluster2=cluster.ini.2
m1=length(cluster1)
m2=length(cluster2)
cluster.ratio=cluster.ratio.ini
iter=0
error=1
fit.susie2=NULL
while(iter<max.iter&error>max.eps){
theta11=theta1
theta22=theta2
Rxysum1=biasterm(RxyList=RxyList,cluster1)
XtX1=matrixMultiply(t(tilde.X[cluster1,]),tilde.X[cluster1,])
XtX1=XtX1-Rxysum1[1:p,1:p]
XtX1=t(XtX1)/2+XtX1/2
Xty1=matrixVectorMultiply(t(tilde.X[cluster1,]),tilde.y[cluster1])-Rxysum1[1:p,1+p]
yty1=sum(tilde.y[cluster1]^2)
adjX1=xtx_positive(XtX1,Xty1)
XtX1=adjX1$XtX
Xty1=adjX1$Xty
Diff_matrix1=diag(p)*0
if(group.penalize==T){
Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1)
}
tryCatch({
fit.susie1=susie_ss(XtX=XtX1+Diff_matrix1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[v],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie1=susie_ss(XtX=XtX1+Diff_matrix1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[v],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta1=coef.susie(fit.susie1)[-1]*(fit.susie1$pip>pip.min)
theta.cs1=group.pip.filter(pip.summary=summary(fit.susie1)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive1=theta.cs1$ind.keep
theta1[-pip.alive1]=0
if(length(cluster2)>min.cluster.size){

Rxysum2=biasterm(RxyList=RxyList,cluster2)
XtX2=matrixMultiply(t(tilde.X[cluster2,]),tilde.X[cluster2,])
XtX2=XtX2-Rxysum2[1:p,1:p]
XtX2=t(XtX2)/2+XtX2/2
Xty2=matrixVectorMultiply(t(tilde.X[cluster2,]),tilde.y[cluster2])-Rxysum2[1:p,1+p]
yty2=sum(tilde.y[cluster2]^2)
adjX2=xtx_positive(XtX2,Xty2)
XtX2=adjX2$XtX
Xty2=adjX2$Xty
Diff_matrix2=diag(p)*0
if(group.penalize==T){
Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2)
}
tryCatch({
fit.susie2=susie_ss(XtX=XtX2+Diff_matrix2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[l],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie2=susie_ss(XtX=XtX2+Diff_matrix2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[l],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta2=coef.susie(fit.susie2)[-1]*(fit.susie2$pip>pip.min)
theta.cs2=group.pip.filter(pip.summary=summary(fit.susie2)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive2=theta.cs2$ind.keep
theta2[-pip.alive2]=0
Diff2=generate_block_matrix(summary(fit.susie2)$vars,length(cluster2)/diag(XtX2),theta2)
}else{
theta2=theta1*0
cluster2=c(1:min.cluster.size)
}
indtheta1=which(theta1!=0)
Diff1=generate_block_matrix(summary(fit.susie1)$vars,length(cluster1)/diag(XtX1),theta1)
if(length(indtheta1)==1){
xtx1=XtX1[indtheta1,indtheta1]
xty1=Xty1[indtheta1]
theta1[indtheta1]=xty1/xtx1
}
if(length(indtheta1)>1){
XtX1=XtX1[indtheta1,indtheta1]+ridge.diff*Diff1[indtheta1,indtheta1]+Diff_matrix1[indtheta1,indtheta1]
Xty1=Xty1[indtheta1]
theta1[indtheta1]=as.vector(solve(XtX1)%*%Xty1)
}
indtheta2=which(theta2!=0)
if(length(indtheta2)==1){
xtx2=XtX2[indtheta2,indtheta2]
xty2=Xty2[indtheta2]
theta2[indtheta2]=xty2/xtx2
}
if(length(indtheta2)>1){
XtX2=XtX2[indtheta2,indtheta2]+ridge.diff*Diff2[indtheta2,indtheta2]+Diff_matrix2[indtheta2,indtheta2]
Xty2=Xty2[indtheta2]
theta2[indtheta2]=as.vector(solve(XtX2)%*%Xty2)
}
sigma1=sum((tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
sigma2=sum((tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
sigma2=max(0.25,sigma2)
sigma1=max(0.25,sigma1)
Voting=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
cluster2=which(Voting$Cluster[,2]==1)
if(length(cluster2)==0){
cluster2=c(1:min.cluster.size)
}
cluster1=which(Voting$Cluster[,1]==1)
m1=length(cluster1)
m2=length(cluster2)
cluster.ratio=c(length(cluster1),length(cluster2))/m
iter=iter+1
if(iter>3){
error=min(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
}
}
Btheta1[,v,l]=theta1
Btheta2[,v,l]=theta2
df1=min(sum(theta1!=0),Lvec[v])
df2=min(sum(theta2!=0),Lvec[l])
Bbic[v,l]=MRFit(by=tilde.y,bX=tilde.X,theta1=theta1,theta2=theta2,cluster1=cluster1,cluster2=cluster2,df1=df1,df2=df2)+(df2+df1)*(log(p*2)*ebic.theta+log(m))+log(m)
SIG1[v,l]=sigma1
SIG2[v,l]=sigma2
M1[v,l]=m1
M2[v,l]=m2
}
}
Bbic=Bbic/m
######################## Final Estimate #################################
vstar=bimin(Bbic)[1]
lstar=bimin(Bbic)[2]
theta1=Btheta1[,vstar,lstar]
theta2=Btheta2[,vstar,lstar]
sigma1=SIG1[vstar,lstar]
sigma2=SIG2[vstar,lstar]
m1=M1[vstar,lstar]
m2=M2[vstar,lstar]
Voting=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
cluster2=which(Voting$Cluster[,2]==1)
cluster1=which(Voting$Cluster[,1]==1)
m1=length(cluster1)
m2=length(cluster2)
iter=0
error=1
theta11=theta1*0
theta22=theta2*0
fit.susie1=fit.susie2=NULL
while(iter<(max.iter*2)&error>max.eps){
theta11=theta1
theta22=theta2
Rxysum1=biasterm(RxyList=RxyList,cluster1)
XtX1=matrixMultiply(t(tilde.X[cluster1,]),tilde.X[cluster1,])
XtX1=XtX1-Rxysum1[1:p,1:p]
XtX1=t(XtX1)/2+XtX1/2
Xty1=matrixVectorMultiply(t(tilde.X[cluster1,]),tilde.y[cluster1])-Rxysum1[1:p,1+p]
yty1=sum(tilde.y[cluster1]^2)
adjX1=xtx_positive(XtX1,Xty1)
XtX1=adjX1$XtX
Xty1=adjX1$Xty
Diff_matrix1=diag(p)*0
if(group.penalize==T){
Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1)
}
tryCatch({
fit.susie1=susie_ss(XtX=XtX1+Diff_matrix1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[vstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie1=susie_ss(XtX=XtX1+Diff_matrix1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[vstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta1=coef.susie(fit.susie1)[-1]*(fit.susie1$pip>pip.min)
theta.cs1=group.pip.filter(pip.summary=summary(fit.susie1)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive1=theta.cs1$ind.keep
theta1[-pip.alive1]=0
if(length(cluster2)>min.cluster.size){

Rxysum2=biasterm(RxyList=RxyList,cluster2)
XtX2=matrixMultiply(t(tilde.X[cluster2,]),tilde.X[cluster2,])-Rxysum2[1:p,1:p]
XtX2=t(XtX2)/2+XtX2/2
Xty2=matrixVectorMultiply(t(tilde.X[cluster2,]),tilde.y[cluster2])-Rxysum2[1:p,1+p]
yty2=sum(tilde.y[cluster2]^2)
adjX2=xtx_positive(XtX2,Xty2)
XtX2=adjX2$XtX
Xty2=adjX2$Xty
Diff_matrix2=diag(p)*0
if(group.penalize==T){
Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2)
}
tryCatch({
fit.susie2=susie_ss(XtX=XtX2+Diff_matrix2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[lstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie2=susie_ss(XtX=XtX2+Diff_matrix2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[lstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta2=coef.susie(fit.susie2)[-1]*(fit.susie2$pip>pip.min)
theta.cs2=group.pip.filter(pip.summary=summary(fit.susie2)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive2=theta.cs2$ind.keep
theta2[-pip.alive2]=0
Diff2=generate_block_matrix(summary(fit.susie2)$vars,length(cluster2)/diag(XtX2),theta2)
}else{
theta2=theta1*0
cluster2=c(1:4)
}
indtheta1=which(theta1!=0)
Diff1=generate_block_matrix(summary(fit.susie1)$vars,length(cluster1)/diag(XtX1),theta1)
if(length(indtheta1)==1){
xtx1=XtX1[indtheta1,indtheta1]
xty1=Xty1[indtheta1]
theta1[indtheta1]=xty1/xtx1
}
if(length(indtheta1)>1){
XtX1=XtX1[indtheta1,indtheta1]+ridge.diff*Diff1[indtheta1,indtheta1]+Diff_matrix1[indtheta1,indtheta1]
Xty1=Xty1[indtheta1]
theta1[indtheta1]=c(solve(XtX1)%*%Xty1)
}
indtheta2=which(theta2!=0)
if(length(indtheta2)==1){
xtx2=XtX2[indtheta2,indtheta2]
xty2=Xty2[indtheta2]
theta2[indtheta2]=xty2/xtx2
}
if(length(indtheta2)>1){
XtX2=XtX2[indtheta2,indtheta2]+ridge.diff*Diff2[indtheta2,indtheta2]+Diff_matrix2[indtheta2,indtheta2]
Xty2=Xty2[indtheta2]
theta2[indtheta2]=c(solve(XtX2)%*%Xty2)
}
sigma1=sum((tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
sigma2=sum((tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
sigma2=max(0.25,sigma2)
sigma1=max(0.25,sigma1)
Voting=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
cluster2=which(Voting$Cluster[,2]==1)
if(length(cluster2)==0){
cluster2=c(1:min.cluster.size)
}
cluster1=which(Voting$Cluster[,1]==1)
m1=length(cluster1)
m2=length(cluster2)
cluster.ratio=c(length(cluster1),length(cluster2))/m
iter=iter+1
if(iter>3){
error=min(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
}
}
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Estimation ends: ",time_to_print," secs\n"))
}
############################### inference #########################
t1=Sys.time()
cat("Bootstrapping starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
names(theta1)=names(theta2)=colnames(bX)
ThetaList1=ThetaList2=matrix(0,sampling.time,p)
colnames(ThetaList1)=colnames(ThetaList2)=colnames(bX)
cluster.index <- as.integer(factor(cluster.index))
cluster_prob <- cluster_prob(cluster.index,LD,alpha=prob_shrinkage_coef,group_size=prob_shrinkage_size)
j=1
if(isLD==T){
cluster_cache <- precompute_cluster_blocks_mixture(
bX = bX,
bXse = bXse,
by = by,
byse = byse,
TC = TC,
cluster.index = cluster.index
)
}
while(j<=sampling.time){
indicator <- FALSE
setTxtProgressBar(pb, j)
tryCatch({
if(isLD==T){
if (sampling.strategy == "bootstrap") {
cluster.sampling <- sample(1:max(cluster.index),
                           size = max(cluster.index),
                           replace = TRUE,
                           prob = cluster_prob)
} else {
cluster.sampling <- sample(1:max(cluster.index),
                           size = 0.5 * max(cluster.index),
                           replace = FALSE,
                           prob = cluster_prob)
}
  cluster.sampling=sort(cluster.sampling)
sampled_blocks <- cluster_cache[cluster.sampling]
indj <- unlist(lapply(sampled_blocks, function(b) b$idx))
indj <- sort(indj)
tilde.Xj <- do.call(rbind, lapply(sampled_blocks, function(b) b$tilde.X))
tilde.yj <- unlist(lapply(sampled_blocks, function(b) b$tilde.y))
bXsej <- do.call(rbind, lapply(sampled_blocks, function(b) b$bXse))
bysej <- unlist(lapply(sampled_blocks, function(b) b$byse))
bXj <- bX[indj, ]
byj <- by[indj]
}else{
if (sampling.strategy == "bootstrap") {
indj <- sample(1:m, size = m, replace = TRUE)
} else {
indj <- sample(1:m, size = 0.5 * m, replace = FALSE)
}
indj=sort(indj)
bXj=tilde.Xj=bX[indj,]
bXsej=bXse[indj,]
byj=tilde.yj=by[indj]
bysej=byse[indj]
}
theta1j=theta1*runif(p,0.95,1.05)*0.95
theta2j=theta2*runif(p,0.95,1.05)*0.95
cluster1j=which(indj%in%cluster1)
cluster2j=which(indj%in%cluster2)
m1j=length(cluster1j)
m2j=length(cluster2j)
sigma1j=sigma1
sigma2j=sigma2
errorj=1
if(sampling.strategy=="bootstrap"){
fit.susie1j=fit.susie1
fit.susie2j=fit.susie2
}else{
fit.susie1j=NULL
fit.susie2j=NULL
}

for(jiter in 1:sampling.iter){
theta_prev1j=theta1j
theta_prev2j=theta2j
Rxysum1j=biasterm(RxyList=RxyList,indj[cluster1j])
XtX1j=matrixMultiply(t(tilde.Xj[cluster1j,]),tilde.Xj[cluster1j,])-Rxysum1j[1:p,1:p]
XtX1j=t(XtX1j)/2+XtX1j/2
Xty1j=matrixVectorMultiply(t(tilde.Xj[cluster1j,]),tilde.yj[cluster1j])-Rxysum1j[1:p,1+p]
yty1j=sum(tilde.yj[cluster1j]^2)
adjX1j=xtx_positive(XtX1j,Xty1j)
XtX1j=adjX1j$XtX
Xty1j=adjX1j$Xty
Diff_matrix1=diag(p)*0
if(group.penalize==T){
Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1j)
}
tryCatch({
fit.susie1j=susie_ss(XtX=XtX1j+Diff_matrix1,Xty=Xty1j,yty=yty1j,n=length(cluster1j),L=Lvec[vstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie1j,estimate_prior_method="EM",coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie1j=susie_ss(XtX=XtX1j+Diff_matrix1,Xty=Xty1j,yty=yty1j,n=length(cluster1j),L=Lvec[vstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie1j,estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta1j=coef.susie(fit.susie1j)[-1]*(fit.susie1j$pip>pip.min)
theta.cs1j=group.pip.filter(pip.summary=summary(fit.susie1j)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive1j=theta.cs1j$ind.keep
theta1j[-pip.alive1j]=0
if(length(cluster2j)>(min.cluster.size/2)){

Rxysum2j=biasterm(RxyList=RxyList,indj[cluster2j])
XtX2j=matrixMultiply(t(tilde.Xj[cluster2j,]),tilde.Xj[cluster2j,])-Rxysum2j[1:p,1:p]
XtX2j=XtX2j/2+t(XtX2j)/2
Xty2j=matrixVectorMultiply(t(tilde.Xj[cluster2j,]),tilde.yj[cluster2j])-Rxysum2j[1:p,1+p]
yty2j=sum(tilde.yj[cluster2j]^2)
adjX2j=xtx_positive(XtX2j,Xty2j)
XtX2j=adjX2j$XtX
Xty2j=adjX2j$Xty
Diff_matrix2=diag(p)*0
if(group.penalize==T){
Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2j)
}
tryCatch({
fit.susie2j=susie_ss(XtX=XtX2j+Diff_matrix2,Xty=Xty2j,yty=yty2j,n=length(cluster2j),L=Lvec[lstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie2j,estimate_prior_method="EM",coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
},error = function(e) {
fit.susie2j=susie_ss(XtX=XtX2j+Diff_matrix2,Xty=Xty2j,yty=yty2j,n=length(cluster2j),L=Lvec[lstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie2j,estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
})
theta2j=coef.susie(fit.susie2j)[-1]*(fit.susie2j$pip>pip.min)
theta.cs2j=group.pip.filter(pip.summary=summary(fit.susie2j)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
pip.alive2j=theta.cs2j$ind.keep
theta2j[-pip.alive2j]=0
Diff2j=generate_block_matrix(summary(fit.susie2j)$vars,length(cluster2j)/diag(XtX2j),theta2j)
}else{
theta2j=theta1j*0
cluster2j=c(1:2)
}
indtheta1j=which(theta1j!=0)
Diff1j=generate_block_matrix(summary(fit.susie1j)$vars,length(cluster1j)/diag(XtX1j),theta1j)
if(length(indtheta1j)==1){
xtx1j=XtX1j[indtheta1j,indtheta1j]
xty1j=Xty1j[indtheta1j]
theta1j[indtheta1j]=xty1j/xtx1j
}
if(length(indtheta1j)>1){
XtX1j=XtX1j[indtheta1j,indtheta1j]+ridge.diff*Diff1j[indtheta1j,indtheta1j]+Diff_matrix1[indtheta1j,indtheta1j]/2
Xty1j=Xty1j[indtheta1j]
theta1j[indtheta1j]=c(solve(XtX1j)%*%Xty1j)
}
indtheta2j=which(theta2j!=0)
if(length(indtheta2j)==1){
xtx2j=XtX2j[indtheta2j,indtheta2j]
xty2j=Xty2j[indtheta2j]
theta2j[indtheta2j]=xty2j/xtx2j
}
if(length(indtheta2j)>1){
XtX2j=XtX2j[indtheta2j,indtheta2j]+ridge.diff*Diff2j[indtheta2j,indtheta2j]+Diff_matrix2[indtheta2j,indtheta2j]/2
Xty2j=Xty2j[indtheta2j]
theta2j[indtheta2j]=c(solve(XtX2j)%*%Xty2j)
}
sigma1j=sum((tilde.yj[cluster1j]-tilde.Xj[cluster1j,]%*%theta1j)^2)/(length(cluster1j)-sum(theta1j!=0))
sigma2j=sum((tilde.yj[cluster2j]-tilde.Xj[cluster2j,]%*%theta2j)^2)/(length(cluster2j)-sum(theta2j!=0))
sigma2j=max(0.25,sigma2j)
sigma1j=max(0.25,sigma1j)
Votingj=cluster_voting(by=tilde.yj,bX=tilde.Xj,cluster.index=cluster.index[indj],theta1=theta1j,theta2=theta2j,sigma1=sigma1j,sigma2=sigma2j,main.cluster.thres=main.cluster.thres,m1=m1j,m2=m2j)
cluster2j=which(Votingj$Cluster[,2]==1)
if(length(cluster2j)==0){
cluster2j=c(1:2)
}
cluster1j=which(Votingj$Cluster[,1]==1)
m1j=length(cluster1j)
m2j=length(cluster2j)
if(jiter>4) errorj=max(norm(theta1j-theta_prev1j,"2"),norm(theta2j-theta_prev2j,"2"))
if(errorj<max.eps) break
}
ThetaVecj=cbind(theta1j,theta2j)
ThetaNormj=colMeans((ThetaVecj-cbind(theta1,theta1))^2)
theta1j=ThetaVecj[,which.min(ThetaNormj)]
theta2j=ThetaVecj[,which.max(ThetaNormj)]

ThetaList1[j, ] <- theta1j
ThetaList2[j, ] <- theta2j
j=j+1
}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, "\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
j <<- j - 1  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
}
close(pb)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Bootstrapping ends: ",time_to_print," secs\n"))
}
indtheta1=which(theta1!=0)
indtheta2=which(theta2!=0)

theta.se1=colSDMAD(ThetaList1)
theta.cov1=covmad(ThetaList1)
theta.se2=colSDMAD(ThetaList2)
theta.cov2=covmad(ThetaList2)

colnames(theta.cov2)=rownames(theta.cov2)=names(theta.se2)=colnames(bX)
colnames(theta.cov1)=rownames(theta.cov1)=names(theta.se1)=colnames(bX)

A=list()
A$theta1=theta1
A$theta2=theta2
A$theta.se1=theta.se1
A$theta.cov1=theta.cov1
A$theta.se2=theta.se2
A$theta.cov2=theta.cov2
A$Bic=Bbic
A$reliability.adjust=r
A$thetalist1=ThetaList1
A$thetalist2=ThetaList2
A$Voting=Voting
A$theta.pip1=colMeans(ThetaList1!=0)
A$theta.pip2=colMeans(ThetaList2!=0)
A$susie.theta1=fit.susie1
A$susie.theta2=fit.susie2
A$Diff1=Diff1
A$Diff2=Diff2
A$Group_Penalty1=Diff_matrix1
A$Group_Penalty2=Diff_matrix2
return(A)
}
