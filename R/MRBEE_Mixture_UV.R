MRBEE_Mixture_UV=function(by,bX,byse,bXse,Rxy,LD=LD,cluster.index=cluster.index,reliability.thres=0.8,sampling.time=100,ebic.theta=1,max.iter=30,max.eps=5e-4,sampling.iter=5){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=length(bX)
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
by=as.vector(TC%*%by)
bX=as.vector(TC%*%bX)
tilde.R=as.matrix(TC%*%LD)
r=reliability.adj.uv(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
p=1
############################ Initial Estimate #######################
fit.init=fit.mixture=regmixEM(y=by,x=bX,k=2,addintercept=T)
max.cluster=ifelse(sum(fit.init$posterior[,1]>0.49)>(m/2),1,2)
cluster.ini.1=which(fit.init$posterior[,max.cluster]>0.49)
cluster.ini.2=setdiff(1:m,cluster.ini.1)
theta1=fit.init$beta[2,max.cluster]
theta2=fit.init$beta[2,setdiff(1:2,max.cluster)]
error=1
iter=1
cluster1=cluster.ini.1
cluster2=cluster.ini.2
cluster.ratio=c(length(cluster1),length(cluster2))/m
iter=0
error=1
while(iter<max.iter&error>max.eps){
theta11=theta1
theta22=theta2
Rxysum1=biasterm(RxyList=RxyList,cluster1)
Rxysum2=biasterm(RxyList=RxyList,cluster2)
XtX1=sum(bX[cluster1]^2)
Xty1=sum(bX[cluster1]*by[cluster1])
theta1=(Xty1-Rxysum1[1:p,p+1])/(XtX1-Rxysum1[1:p,1:p])
if(length(cluster2)>0.5){
XtX2=sum(bX[cluster2]^2)
Xty2=sum(bX[cluster2]*by[cluster2])
theta2=(Xty2-Rxysum2[1:p,p+1])/(XtX2-Rxysum2[1:p,1:p])
}else{
theta2=0
cluster2=sample(m,0.5)
}
sigma1=sum((by[cluster1]-bX[cluster1]*theta1)^2)/(length(cluster1)-1)
sigma2=sum((by[cluster2]-bX[cluster2]*theta2)^2)/(length(cluster2)-1)
sigma2=max(0.1,sigma2)
sigma1=max(0.1,sigma1)
Voting=cluster_voting_uv(by=by,bX=bX,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2)
cluster2=which(Voting$Cluster[,2]==1)
cluster1=which(Voting$Cluster[,1]==1)
cluster.ratio=c(length(cluster1),length(cluster2))/m
res1=by[cluster1]-bX[cluster1]*theta1
res2=by[cluster2]-bX[cluster2]*theta2
rss1=sum(res1^2)/(length(cluster1)-2)
rss2=sum(res2^2)/(length(cluster2)-2)
Bic=log(rss1)*length(cluster1)+log(rss2)*length(cluster2)+2*p*(log(2*p)*ebic.theta+log(m))
Bic=Bic/m
iter=iter+1
if(iter>3){
error=min(abs(theta1-theta11),abs(theta2-theta22))
}
}
############################### inference #########################
names(theta1)=names(theta2)=colnames(bX)
ThetaList1=ThetaList2=c(1:sampling.time)*0
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
byj=by[indj]
bXj=bX[indj]
theta1j=theta1
theta2j=theta2
cluster1j=which(indj%in%cluster1)
cluster2j=which(indj%in%cluster2)
for(iterj in 1:sampling.iter){
Rxysum1j=biasterm(RxyList=RxyList,indj[cluster1j])
XtX1j=sum(bXj[cluster1j])-Rxysum1j[1:p,1:p]
Xty1j=sum(bXj[cluster1j]*byj[cluster1j])
theta1j=(Xty1j-Rxysum1j[1:p,p+1])/XtX1j
if(length(cluster2j)>0.5){
Rxysum2j=biasterm(RxyList=RxyList,indj[cluster2j])
XtX2j=sum(bXj[cluster2j])-Rxysum2j[1:p,1:p]
Xty2j=sum(bXj[cluster2j]*byj[cluster2j])
theta2j=(Xty2j-Rxysum2j[1:p,p+1])/XtX2j
}else{
theta2j=0*theta1j
cluster2j=sample(length(indj),0.5)
}
sigma1j=sum((byj[cluster1j]-bXj[cluster1j,]*theta1j)^2)/(length(cluster1j)-1)
sigma2j=sum((byj[cluster2j]-bXj[cluster2j,]*theta2j)^2)/(length(cluster2j)-1)
sigma2j=max(0.1,sigma2j)
sigma1j=max(0.1,sigma1j)
Votingj=cluster_voting_uv(by=byj,bX=bXj,cluster.index=cluster.index[indj],theta1=theta1j,theta2=theta2j,sigma1=sigma1j,sigma2=sigma2j)
cluster2j=which(Votingj$Cluster[,2]==1)
cluster1j=which(Votingj$Cluster[,1]==1)
}
ThetaVecj=c(theta1j,theta2j)
ThetaNormj=abs(ThetaVecj-c(theta1,theta1))
theta1j=ThetaVecj[which.min(ThetaNormj)]
theta2j=ThetaVecj[which.max(ThetaNormj)]
ThetaList1[j] <- theta1j
ThetaList2[j] <- theta2j
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

theta.se1=mad(ThetaList1)
theta.cov1=mad(ThetaList1)^2
theta.se2=mad(ThetaList2)
theta.cov2=mad(ThetaList2)^2

names(theta.cov2)=names(theta.se2)=names(bX)
names(theta.cov1)=names(theta.se1)=names(bX)
A=list()
A$theta1=theta1
A$theta2=theta2
A$theta.se1=theta.se1
A$theta.cov1=theta.cov1
A$theta.se2=theta.se2
A$theta.cov2=theta.cov2
A$reliability.adjust=r
A$thetalist2=ThetaList2
A$thetalist2=ThetaList1
A$cluster1=cluster1
A$cluster2=cluster2
A$fit.mixture=fit.mixture
A$theta.pratt1=getPratt.uv(bX=bX[cluster1],by=by[cluster1],bXse=bXse[cluster1],byse=byse[cluster1],theta=theta1,Rxy=Rxy)
A$theta.pratt2=getPratt.uv(bX=bX[cluster2],by=by[cluster2],bXse=bXse[cluster2],byse=byse[cluster2],theta=theta2,Rxy=Rxy)
A$IsIPOD=F
A$Bic=Bic
return(A)
}
