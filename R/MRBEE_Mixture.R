MRBEE_Mixture=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),main.cluster.thres=0.45,min.cluster.size=5,reliability.thres=0.8,sampling.time=100,ebic.theta=1,max.iter=30,max.eps=5e-4,sampling.iter=5,verbose=T,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,LDSC=NULL,Omega=NULL,prob_shrinkage_coef=0.5,prob_shrinkage_size=4,sampling.strategy="bootstrap"){
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
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=LDSC,Omega=Omega)
############################ Initial Estimate #######################
fit.init=fit.mixture=regmixEM(y=tilde.y,x=tilde.X,k=2,addintercept=F,maxit=300,epsilon=max.eps)
max.cluster=ifelse(sum(fit.init$posterior[,1]>main.cluster.thres)>(m/2),1,2)
cluster1=which(fit.init$posterior[,max.cluster]>main.cluster.thres)
cluster2=setdiff(1:m,cluster1)
if(length(cluster2)<min.cluster.size){
cluster2=c(1:min.cluster.size)
}
m1=length(cluster1)
m2=length(cluster2)
theta1=fit.init$beta[,max.cluster]
theta2=fit.init$beta[,setdiff(1:2,max.cluster)]
sigma1=sum((tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
sigma2=sum((tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
sigma2=max(0.25,sigma2)
sigma1=max(0.25,sigma1)
Voting=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
cluster2=which(Voting$Cluster[,2]==1)
cluster1=which(Voting$Cluster[,1]==1)
m1=length(cluster1)
m2=length(cluster2)
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Initialization ends: ",time_to_print," secs\n"))
}
############################## Tuning Parameter ######################
t1=Sys.time()
iter=0
error=1
while(iter<max.iter&error>max.eps){
theta11=theta1
theta22=theta2
Rxysum1=biasterm(RxyList=RxyList,cluster1)
XtX1=matrixMultiply(t(tilde.X[cluster1,]),tilde.X[cluster1,])-Rxysum1[1:p,1:p]
Xty1=matrixVectorMultiply(t(tilde.X[cluster1,]),tilde.y[cluster1])-Rxysum1[1:p,p+1]
adjX1=xtx_positive(XtX1,Xty1)
XtX1=adjX1$XtX
Xty1=adjX1$Xty
Diff_matrix1=diag(p)*0
if(group.penalize==T){
Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1)
}
theta1=c(solve(XtX1+Diff_matrix1)%*%Xty1)
if(length(cluster2)>min.cluster.size){
Rxysum2=biasterm(RxyList=RxyList,cluster2)
XtX2=matrixMultiply(t(tilde.X[cluster2,]),tilde.X[cluster2,])-Rxysum2[1:p,1:p]
Xty2=matrixVectorMultiply(t(tilde.X[cluster2,]),tilde.y[cluster2])-Rxysum2[1:p,p+1]
adjX2=xtx_positive(XtX2,Xty2)
XtX2=adjX2$XtX
Xty2=adjX2$Xty
Diff_matrix2=diag(p)*0
if(group.penalize==T){
Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2)
}
theta2=c(solve(XtX2+Diff_matrix2)%*%Xty2)
}else{
theta2=0*theta1
cluster2=c(1:min.cluster.size)
}
sigma1=sum((tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-p)
sigma2=sum((tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-p)
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
iter=iter+1
if(iter>3){
error=min(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
}
Bic=MRFit(by=tilde.y,bX=tilde.X,theta1=theta1,theta2=theta2,cluster1=cluster1,cluster2=cluster2,df1=p,df2=p)+2*p*(log(p*2)*ebic.theta+log(m))+log(m)
Bic=Bic/m
}
t2=Sys.time()
time_to_print=round(difftime(t2, t1, units = "secs"),3)
if(verbose==T){
cat(paste0("Estimation ends: ",time_to_print," secs\n"))
}
############################### inference #########################
t1=Sys.time()
names(theta1)=names(theta2)=colnames(bX)
ThetaList1=ThetaList2=matrix(0,sampling.time,p)
colnames(ThetaList1)=colnames(ThetaList2)=colnames(bX)
cat("Bootstrapping starts:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
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
errorj=1
for(jiter in 1:sampling.iter){
theta_prev1j=theta1j
theta_prev2j=theta2j
Rxysum1j=biasterm(RxyList=RxyList,indj[cluster1j])
XtX1j=matrixMultiply(t(tilde.Xj[cluster1j,]),tilde.Xj[cluster1j,])-Rxysum1j[1:p,1:p]
Xty1j=matrixVectorMultiply(t(tilde.Xj[cluster1j,]),tilde.yj[cluster1j])-Rxysum1j[1:p,p+1]
adjX1j=xtx_positive(XtX1j,Xty1j)
XtX1j=adjX1j$XtX
Xty1j=adjX1j$Xty
Diff_matrix1=diag(p)*0
if(group.penalize==T){
Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1j)
}
theta1j=c(solve(XtX1j+Diff_matrix1/2)%*%Xty1j)
if(length(cluster2j)>(min.cluster.size/2)){
Rxysum2j=biasterm(RxyList=RxyList,indj[cluster2j])
XtX2j=matrixMultiply(t(tilde.Xj[cluster2j,]),tilde.Xj[cluster2j,])-Rxysum2j[1:p,1:p]
Xty2j=matrixVectorMultiply(t(tilde.Xj[cluster2j,]),tilde.yj[cluster2j])-Rxysum2j[1:p,p+1]
adjX2j=xtx_positive(XtX2j,Xty2j)
XtX2j=adjX2j$XtX
Xty2j=adjX2j$Xty
if(group.penalize==T){
Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2j)
}
theta2j=c(solve(XtX2j+Diff_matrix2/2)%*%Xty2j)
}else{
theta2j=0*theta1j
cluster2j=c(1:4)
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
A$reliability.adjust=r
A$thetalist1=ThetaList1
A$thetalist2=ThetaList2
A$cluster1=cluster1
A$cluster2=cluster2
A$fit.mixture=fit.mixture
A$IsIPOD=F
A$Voting=Voting
A$Bic=Bic
A$Group_Penalty1=Diff_matrix1
A$Group_Penalty2=Diff_matrix2
return(A)
}
