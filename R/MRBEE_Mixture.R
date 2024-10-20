MRBEE_Mixture=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),main.cluster.thres=0.48,min.cluster.size=5,reliability.thres=0.8,sampling.time=100,robust.se=T,ebic.theta=1,ebic.gamma=1,max.iter=30,max.eps=5e-4,sampling.iter=5,tau=tau,step.size=0.8){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
if(LD[1]!="identity"){
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
}else{
LD=Theta=TC=diag(m)
}
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
tilde.R=as.matrix(TC%*%LD)
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
fit.init=fit.mixture=regmixEM(y=tilde.y,x=tilde.X,k=2,addintercept=F,maxit=300)
max.cluster=ifelse(sum(fit.init$posterior[,1]>main.cluster.thres)>(m/2),1,2)
cluster1=which(fit.init$posterior[,max.cluster]>main.cluster.thres)
cluster2=setdiff(1:m,cluster1)
########################### Check if applied mixture model #######################
if(length(cluster2)>=min.cluster.size){
theta1=fit.init$beta[,max.cluster]
theta2=fit.init$beta[,setdiff(1:2,max.cluster)]
sigma1=sum((tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
sigma2=sum((tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
sigma2=max(0.1,sigma2)
sigma1=max(0.1,sigma1)
Voting=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres)
cluster2=which(Voting$Cluster[,2]==1)
cluster1=which(Voting$Cluster[,1]==1)
cluster.ratio.ini=c(length(cluster1),length(cluster2))/m
iter=0
error=1
gamma=by*0
while(iter<max.iter&error>max.eps){
theta11=theta1
theta22=theta2
indvalid=which(gamma==0)
valid1=intersect(indvalid,cluster1)
valid2=intersect(indvalid,cluster2)
Rxysum1=biasterm(RxyList=RxyList,valid1)
gamma.y=tilde.y-matrixVectorMultiply(tilde.R,gamma)
XtX1=matrixMultiply(t(tilde.X[cluster1,]),tilde.X[cluster1,])
Xty1=matrixVectorMultiply(t(tilde.X[cluster1,]),gamma.y[cluster1])
theta1=c(solve(XtX1-Rxysum1[1:p,1:p])%*%(Xty1-Rxysum1[1:p,p+1]))
if(length(cluster2)>min.cluster.size){
Rxysum2=biasterm(RxyList=RxyList,valid2)
XtX2=matrixMultiply(t(tilde.X[cluster2,]),tilde.X[cluster2,])
Xty2=matrixVectorMultiply(t(tilde.X[cluster2,]),gamma.y[cluster2])
theta2=c(solve(XtX2-Rxysum2[1:p,1:p])%*%(Xty2-Rxysum2[1:p,p+1]))
}else{
theta2=0*theta1
cluster2=sample(m,min.cluster.size)
}
res.theta=as.vector(tilde.y-tilde.R%*%gamma)
sigma1=sum((res.theta[cluster1]-tilde.X[cluster1,]%*%theta1)^2)/(length(cluster1)-p)
sigma2=sum((res.theta[cluster2]-tilde.X[cluster2,]%*%theta2)^2)/(length(cluster2)-p)
sigma2=max(0.1,sigma2)
sigma1=max(0.1,sigma1)
Voting=cluster_voting(by=res.theta,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres)
cluster2=which(Voting$Cluster[,2]==1)
cluster1=which(Voting$Cluster[,1]==1)
grad=tilde.y
grad[cluster1]=tilde.y[cluster1]-tilde.X[cluster1,]%*%theta1
grad[cluster2]=tilde.y[cluster2]-tilde.X[cluster2,]%*%theta2
grad=matrixVectorMultiply(t(tilde.R),-grad+matrixVectorMultiply(tilde.R,gamma))
gamma=gamma-grad*step.size
gamma=mcp(gamma,tau)
iter=iter+1
if(iter>3){
error=min(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
}
Bic=MRFit(by=as.vector(tilde.y-tilde.R%*%gamma),bX=tilde.X,theta1=theta1,theta2=theta2,cluster1=cluster1,cluster2=cluster2,df1=p,df2=p)+2*p*(log(p*2)*ebic.theta+log(m))+log(m)*(1+ebic.gamma)*sum(gamma!=0)
Bic=Bic/m
}

############################### inference #########################
names(theta1)=names(theta2)=colnames(bX)
ThetaList1=ThetaList2=matrix(0,sampling.time,p)
colnames(ThetaList1)=colnames(ThetaList2)=colnames(bX)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
j=1
while(j<=sampling.time){
indicator <- FALSE
setTxtProgressBar(pb, j)
tryCatch({
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
tilde.yj=tilde.y[indj]
tilde.Xj=tilde.X[indj,]
tilde.Rj=tilde.R[indj,indj]
theta1j=theta1*runif(1,0.95,1.05)*0.95
theta2j=theta2*runif(1,0.95,1.05)*0.95
cluster1j=which(indj%in%cluster1)
cluster2j=which(indj%in%cluster2)
gammaj=tilde.yj*0
for(iterj in 1:sampling.iter){
indvalidj=which(gamma==0)
valid1j=intersect(indvalidj,cluster1j)
valid2j=intersect(indvalidj,cluster2j)
gamma.yj=tilde.y-matrixVectorMultiply(tilde.Rj,gammaj)
Rxysum1j=biasterm(RxyList=RxyList,indj[valid1j])
XtX1j=matrixMultiply(t(tilde.Xj[cluster1j,]),tilde.Xj[cluster1j,])
Xty1j=matrixVectorMultiply(t(tilde.Xj[cluster1j,]),gamma.yj[cluster1j])
theta1j=c(solve(XtX1j-Rxysum1j[1:p,1:p])%*%(Xty1j-Rxysum1j[1:p,p+1]))
if(length(cluster2j)>(min.cluster.size/2)){
Rxysum2j=biasterm(RxyList=RxyList,indj[valid2j])
XtX2j=matrixMultiply(t(tilde.Xj[cluster2j,]),tilde.Xj[cluster2j,])
Xty2j=matrixVectorMultiply(t(tilde.Xj[cluster2j,]),gamma.yj[cluster2j])
theta2j=c(solve(XtX2j-Rxysum2j[1:p,1:p])%*%(Xty2j-Rxysum2j[1:p,p+1]))
}else{
theta2j=0*theta1j
cluster2j=sample(length(indj),max(round(min.cluster.size/2),1))
}
res.thetaj=as.vector(tilde.yj-tilde.Rj%*%gammaj)
sigma1j=sum((res.thetaj[cluster1j]-tilde.Xj[cluster1j,]%*%theta1j)^2)/(length(cluster1j)-sum(theta1j!=0))
sigma2j=sum((res.thetaj[cluster2j]-tilde.Xj[cluster2j,]%*%theta2j)^2)/(length(cluster2j)-sum(theta2j!=0))
sigma2j=max(0.1,sigma2j)
sigma1j=max(0.1,sigma1j)
Votingj=cluster_voting(by=res.thetaj,bX=tilde.Xj,cluster.index=cluster.index[indj],theta1=theta1j,theta2=theta2j,sigma1=sigma1j,sigma2=sigma2j,main.cluster.thres=main.cluster.thres)
cluster2j=which(Votingj$Cluster[,2]==1)
cluster1j=which(Votingj$Cluster[,1]==1)
gradj=tilde.yj
gradj[cluster1j]=tilde.yj[cluster1j]-tilde.Xj[cluster1j,]%*%theta1j
gradj[cluster2j]=tilde.yj[cluster2j]-tilde.Xj[cluster2j,]%*%theta2j
gradj=matrixVectorMultiply(t(tilde.Rj),-gradj+matrixVectorMultiply(tilde.Rj,gammaj))
gammaj=gammaj-gradj*step.size
gammaj=mcp(gammaj,tau)
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
indtheta1=which(theta1!=0)
indtheta2=which(theta2!=0)
theta.se1=colSD(ThetaList1)
theta.cov1=cov(ThetaList1)
theta.se2=colSD(ThetaList2)
theta.cov2=cov(ThetaList2)
if(robust.se==T){
theta.se1=colSDMAD(ThetaList1)
theta.cov1=covmad(ThetaList1)
theta.se2=colSDMAD(ThetaList2)
theta.cov2=covmad(ThetaList2)
}
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
A$theta.pratt1=getPratt(bX=bX[cluster1,],by=by[cluster1],bXse=bXse[cluster1,],byse=byse[cluster1],Theta=solve(LD[cluster1,cluster1]),theta=theta1,Rxy=Rxy)
A$theta.pratt2=getPratt(bX=bX[cluster2,],by=by[cluster2],bXse=bXse[cluster2,],byse=byse[cluster2],Theta=solve(LD[cluster2,cluster2]),theta=theta2,Rxy=Rxy)
A$IsIPOD=F
A$gamma=gamma
A$Bic=Bic
return(A)
}else{
A=list()
A$IsIPOD=T
return(A)
}
}
