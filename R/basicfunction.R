MRFit=function(by,bX,theta1,theta2,cluster1,cluster2,df1,df2,ebic.en=1){
m=length(by)
res1=by[cluster1]-bX[cluster1,]%*%theta1
res2=by[cluster2]-bX[cluster2,]%*%theta2
m1=length(cluster1)
m2=length(cluster2)
rss=sum(res1^2)+sum(res2^2)
rss=rss/(m-df1-df2)
alpha1=m1/m
alpha2=m2/m
en=(-m1*log(alpha1)-m2*log(alpha2))*2
if(ebic.en>0){
return(log(rss)*m+(ebic.en+log(m))*log(en))
}else{
return(log(rss)*m)
}
}

susie_xQTL_resampling=function(LD,alpha,mu,mu2,pip.thres=0.5,sampling=100){
L=dim(mu)[1]
n=dim(mu)[2]
a=matrix(0,sampling,n)
for(i in 1:sampling){
N=mu*0
for(j in 1:L){
phi=c(rmultinom(1,1,alpha[j,]))
varr=c(mu2[j,])-c(mu[j,])^2
varr[varr<0]=0
varr=sqrt(varr)
N[j,]=(rnorm(n=n,mean=0,sd=1)*varr+c(mu[j,]))*phi
}
a[i,]=colSums(N)
}
pip=colMeans(a!=0)
a[,which(pip<pip.thres)]=0
b=a%*%LD
return(list(pip=pip,mean_joint=colMeans(a),sd_joint=colSD(a),cov_joint=cov(a),mean_margin=colMeans(b),sd_margin=colSD(b),cov_margin=cov(b)))
}

vec=function(A){
return(as.vector(A))
}

soft=function(a,b){
c=abs(a)-b
c[c<0]=0
c=c*sign(a)
return(c)
}

#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen
positiveadj=function(A,min.eps=0.001){
a=matrixEigen(A)
d=c(a$values)
d[d<min.eps]=min.eps
B=matrixMultiply(a$vectors,t(a$vectors)*d)
return(B)
}

#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen
positiveinv=function(A,min.eps=0){
a=matrixEigen(A)
d=c(a$values)
d1=1/d
d1[d<min.eps]=min.eps
B=matrixMultiply(a$vectors,t(a$vectors)*d1)
B=t(B)/2+B/2
return(B)
}

dmcp=function(x,lam,a=3){
b=lam-abs(x)/a
b[abs(x/a)>lam]=0
return(b)
}

mcp=function(x,lam,a=3){
b=abs(x)
z=soft(x,lam)/(1-1/a)
z[which(b>(a*lam))]=x[which(b>(a*lam))]
return(z)
}

groupmcp=function(x, cluster.index, lambda, a = 3) {
J <- max(cluster.index)
z <- x
a_const <- a / (a - 1)
x_split <- split(x, cluster.index)
z_split <- lapply(x_split, function(xj) {
delta <- sqrt(mean(xj^2))
zj <- max(1 - lambda / delta, 0) * xj
zj <- a_const * zj
if (delta > (a * lambda)) {
zj <- xj
}
return(zj)
})
z <- unsplit(z_split, cluster.index)
return(z)
}


trace=function(A){
a=sum(diag(A))
return(a)
}

groupdmcp=function(x,lambda,cluster.index,a=3){
lambda1=x*0
J=max(cluster.index)
for(j in 1:J){
indj=which(cluster.index==j)
xj=x[indj]
xj=mean(abs(xj))
lambda1[indj]=dmcp(xj,lambda,a=3)
}
return(lambda1)
}

#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen
matrixsqrt=function(A){
fit=matrixEigen(A)
d=c(fit$value)
d1=d*0
d1[d>0]=1/d[d>0]
d=sqrt(d)
d1=sqrt(d1)
A=matrixMultiply(fit$vector,t(fit$vector)*d)
B=matrixMultiply(fit$vector,t(fit$vector)*d1)
C=list(w=A,wi=B,eigenfit=fit)
return(C)
}

IVweight=function(byse,bXse,Rxy){
bZse=cbind(bXse,byse)
p=dim(bZse)[2]
n=dim(bZse)[1]
RxyList=array(0,c(n,p,p))
for(i in 1:n){
s=bZse[i,]
RxyList[i,,]=t(t(Rxy)*s)*s
}
return(RxyList)
}

biasterm=function(RxyList,indvalid){
X=RxyList[1,,]*0
n=dim(RxyList)[1]
for(i in indvalid){
X=X+RxyList[i,,]
}
return(X)
}

standardized.residual=function(res,RxyList,theta,adjust=1){
n=length(res)
tilde.theta=c(theta,-1)
vars=res
p=length(theta)
for(i in 1:n){
Rxyi=RxyList[i,,]
adjusti=sqrt(c(rep(adjust,p),1))
Rxyi=t(t(Rxyi)*adjusti)*adjusti
vars[i]=max(0.05,sum(tilde.theta*c(Rxyi%*%tilde.theta)))
}
return(res/sqrt(vars))
}

checkzero=function(A,nonzero.frac){
p=dim(A)[2]
for(i in 1:p){
a=mean(A[,i]!=0)-nonzero.frac[i]
if(a<0){
b=abs(a)
c=sample(which(A[,i]==0),b*length(which(A[,i]==0)),replace=F)
A[c,i]=1
}
}
return(A)
}

verrorvar=function(Rbb,Suu,Suv,hx,theta,hy,pleiotropy.var=0){
p=dim(Rbb)[2]
Sigmabb=sqrt(hx)%*%Rbb%*%sqrt(hx)
Sigmauu=sqrt(diag(p)-hx)%*%Suu%*%sqrt(diag(p)-hx)
Sigmaxx=Sigmabb+Sigmauu
plevar=as.numeric(pleiotropy.var*t(theta)%*%Sigmabb%*%theta)
c=t(theta)%*%Sigmaxx%*%theta+plevar-t(theta)%*%Sigmabb%*%theta/hy*(1+pleiotropy.var)
c=c[1,1]
b=2*t(theta)%*%sqrt(diag(p)-hx)%*%Suv;b=b[1,1]
svv=(-b+sqrt(b^2-4*c))/2
syy=svv^2+t(theta)%*%Sigmaxx%*%theta+plevar+2*t(theta)%*%sqrt(diag(p)-hx)%*%Suv*svv
syy=syy[1]
return(syy)
}

block_subsampling <- function(cluster_result, rho) {
unique_clusters <- unique(cluster_result)
num_selected_clusters <- ceiling(length(unique_clusters) * rho)
selected_clusters <- sample(unique_clusters, num_selected_clusters)
selected_indices <- which(cluster_result %in% selected_clusters)
return(selected_indices)
}

colSD=function(A){
a=A[1,]
for(i in 1:ncol(A)){
a[i]=sd(A[,i])
}
return(a)
}

colSDMAD=function(A){
a=A[1,]
for(i in 1:ncol(A)){
a[i]=mad(A[,i])
}
return(a)
}

covmad=function(A){
R=cor(A,method="spearman")
R=2*sin(R/6*pi)
R[is.na(R)]=0
D=colSDMAD(A)
S=t(t(R)*D)*D
return(S)
}

subsample_indices <- function(cluster.index, sampling.ratio) {
unique_clusters <- unique(cluster.index)
n_subsample <- round(length(unique_clusters) * sampling.ratio)
subsampled_clusters <- sample(unique_clusters, n_subsample)
subsample_indices <- which(cluster.index %in% subsampled_clusters)
return(subsample_indices)
}

generate_default_cluster <- function(n) {
max_category <- ceiling(n / 2)
cluster_index <- rep(1:max_category, each = 2, length.out = n)
return(cluster_index)
}

bimin=function(mat){
min_element <- min(mat)
min_indices <- which(mat == min_element, arr.ind = TRUE)
if (nrow(min_indices) > 1) {
min_indices <- min_indices[nrow(min_indices), ]
}
return(min_indices)
}

parametric.bootstrap=function(bX,bXse,Rxy,theta,LD,RC,var.inf=0){
n=nrow(bX)
p=ncol(bX)
E=matrix(0,n,p+1)
for(i in 1:n){
sei=c(bXse[i,],1)
Vi=t(t(Rxy)*sei)*sei
e=MASS::mvrnorm(n=1,mu=rep(0,p+1),Sigma=Vi)
E[i,]=e
}
hatX=bX+matrixMultiply(RC,E[,1:p])
haty=matrixVectorMultiply(bX,theta)+matrixVectorMultiply(RC,E[,p+1])
if(var.inf>0){
haty=haty+matrixVectorMultiply(LD,rnorm(n=n,mean=0,sd=sqrt(var.inf)))
}
return(list(hatX=hatX,haty=haty))
}

reliability.adj.uv=function(bx,bxse,Theta="identity",thres=0.7){
if(Theta[1]=="identity"){
total.var=mean(bx^2)
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
r=1
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}else{
r=1
Theta=as.matrix(Theta)
total.var=mean(bx*c(Theta%*%bx))
error.var=mean(bxse^2)
reliability=(total.var-error.var)/total.var
if(reliability<thres){
r=total.var/error.var*(1-thres)
}
r=sqrt(r)
}
return(r)
}

reliability.adj=function(bX,bXse,Theta="identity",thres=0.7){
if(Theta[1]=="identity"){
p=ncol(bX)
r=rep(1,p)
total.var=colMeans(bX^2)
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}else{
p=ncol(bX)
r=rep(1,p)
m=length(bX[,1])
Theta=as.matrix(Theta)
total.var=as.vector(diag(t(bX)%*%Theta%*%bX)/m)
error.var=colMeans(bXse^2)
reliability=(total.var-error.var)/total.var
ind=which(reliability<thres)
if(length(ind)>0){
r[ind]=total.var[ind]/error.var[ind]*(1-thres)
}
r=sqrt(r)
}
return(r)
}

#' @importFrom FDRestimation p.fdr
imrpdetect=function(x,theta,RxyList,indvalid,var.est="variance",FDR=T,adjust.method="Sidak"){
p=length(theta)
if(var.est=="robust"){
varx=stats::mad(x[indvalid])^2
}
if(var.est=="variance"){varx=stats::var(x[indvalid])}
if(var.est=="ordinal"){
varx=x*0
for(i in 1:length(x)){
varx[i]=c(RxyList[i,p+1,p+1]+t(theta)%*%RxyList[i,1:p,1:p]%*%theta-2*sum(theta*RxyList[i,p+1,1:p]))
}
}
pv=stats::pchisq(x^2/varx,1,lower.tail=F)
if(FDR==T){
pv=p.fdr(pvalues=pv,adjust.method=adjust.method)$fdrs
}
return(as.vector(pv))
}

validadj <- function(vector1, vector2, tau){
diff <- length(vector2) / length(vector1)
if (diff < tau) {
missing_indices <- setdiff(1:length(vector1), vector2)
sorted_missing_indices <- missing_indices[order(vector1[missing_indices])]
num_to_add <- ceiling(tau * length(vector1)) - length(vector2)
vector2 <- c(vector2, sorted_missing_indices[1:num_to_add])
}
return(vector2)
}

subjecth2=function(A,h2){
p=ncol(A)
for(i in 1:p){
A[,i]=A[,i]/sqrt(sum(A[,i]^2))*sqrt(h2[i])
}
return(A)
}

colMed=function(A){
p=ncol(A)
a=c(1:p)
for(i in 1:p){
a[i]=median(A[,i])
}
return(a)
}

getPratt=function(bX,by,bXse,byse,Theta="identity",theta,Rxy){
p=length(theta)
n=length(by)
beta=xdy=theta
if(Theta[1]!="identity"){
Thetay=as.vector(Theta%*%by)
yty=sum(by*Thetay)-n*Rxy[p+1,p+1]
for(i in 1:p){
bx=bX[,i]
xtxi=sum(bx*(Theta%*%bx))-Rxy[i,i]*sum(bXse[,i]^2)
xdy[i]=sqrt(xtxi/yty)
beta[i]=(sum(bx*Thetay)-Rxy[i,p+1]*sum(bXse[,i]))/sqrt(xtxi)/sqrt(yty)
}
Pratt=xdy*beta*theta
}else{
yty=sum(by^2)-n*Rxy[p+1,p+1]
for(i in 1:p){
bx=bX[,i]
xtxi=sum(bx^2)-Rxy[i,i]*sum(bXse[,i]^2)
xdy[i]=sqrt(xtxi/yty)
beta[i]=(sum(bx*by)-Rxy[i,p+1]*sum(bXse[,i]))/sqrt(xtxi)/sqrt(yty)
}
Pratt=xdy*beta*theta
}
return(Pratt)
}

pleiotropyPratt=function(by,pleiotropy,Theta="identity",LD="identity"){
if(Theta[1]!="identity"){
n=length(by)
ptp=sum(pleiotropy*(LD%*%pleiotropy))
pty=sum(pleiotropy*by)
yty=sum(by*(Theta%*%by))-n
thetaadj=1*sqrt(ptp)/sqrt(yty)
beta=pty/sqrt(yty)/sqrt(ptp)
}else{
n=length(by)
ptp=sum(pleiotropy^2)
pty=sum(pleiotropy*by)
yty=sum(by^2)-n
thetaadj=1*sqrt(ptp)/sqrt(yty)
beta=pty/sqrt(yty)/sqrt(ptp)
}
return(thetaadj*beta)
}

getPratt.uv=function(bX,by,bXse,byse,Theta="identity",theta,Rxy){
p=length(theta)
n=length(by)
beta=xdy=theta
bx=bX
if(Theta[1]!="identity"){
Thetay=as.vector(Theta%*%by)
yty=sum(by*Thetay)-n
xtxi=sum(bX*(Theta%*%bX))-Rxy[1,1]*sum(bXse^2)
xdy=sqrt(xtxi/yty)
beta=(sum(bx*Thetay)-Rxy[1,2]*sum(bXse))/sqrt(xtxi)/sqrt(yty)
Pratt=xdy*beta*theta
}else{
yty=sum(by^2)-n
xtxi=sum(bX^2)-Rxy[1,1]*sum(bXse^2)
xdy=sqrt(xtxi/yty)
beta=(sum(bX*by)-Rxy[1,2]*sum(bXse))/sqrt(xtxi)/sqrt(yty)
Pratt=xdy*beta*theta
}
return(Pratt)
}

mrvariance=function(bXse,byse,theta,Rxy){
Rxy=as.matrix(Rxy)
bZse=cbind(bXse,byse)
vartheta=c(theta,-1)
n=length(byse)
var.res=byse
for(i in 1:n){
bzse=bZse[i,]
Zxy=t(t(Rxy)*bzse)*bzse
var.res[i]=sum(vartheta*(Zxy%*%vartheta))
}
return(var.res)
}

AR1_block_threshold = function(res, cluster.index, block.rho = block.rho) {
S = diag(res) * 0  # Initialize a zero matrix with dimensions matching 'res'

# Step 1: Compute within-group covariance
for (ii in unique(cluster.index)) {
indii = which(cluster.index == ii)
S[indii, indii] = outer(res[indii], res[indii])
}

# Step 2: Compute between-adjacent-group covariance
# Assuming that 'cluster.index' is sorted in the order of the groups
if(block.rho!=0){
for (ii in unique(cluster.index)) {
indii = which(cluster.index == ii)
indjj = which(cluster.index == ii + 1)
if (length(indjj) > 0) {
S[indii, indjj] = outer(res[indii], res[indjj]) * block.rho
S[indjj, indii] = outer(res[indjj], res[indii]) * block.rho
}
}
}

S = Matrix(S, sparse = TRUE)
}

cluster_voting=function(by,bX,cluster.index,theta1,theta2,sigma1,sigma2,main.cluster.thres=0.49){
res1=as.vector(by-bX%*%theta1)
res2=as.vector(by-bX%*%theta2)
J=max(cluster.index)
Weight=matrix(0,length(by),2)
Vote=matrix(0,length(by),2)
for(j in 1:J){
indj=which(cluster.index==j)
logL1 <- sum(dnorm(res1[indj], mean = 0, sd = sqrt(sigma1), log = TRUE))
logL2 <- sum(dnorm(res2[indj], mean = 0, sd = sqrt(sigma2), log = TRUE))
C = max(logL1, logL2)
logL1_adj <- logL1 - C
logL2_adj <- logL2 - C
L1_adj <- exp(logL1_adj)
L2_adj <- exp(logL2_adj)
w1 <- L1_adj / (L1_adj + L2_adj)
w2 <- L2_adj / (L1_adj + L2_adj)
Weight[indj,1]=w1
Weight[indj,2]=w2
Vote[indj,1]=ifelse(w1>main.cluster.thres,1,0)
Vote[indj,2]=ifelse(w2>main.cluster.thres,1,0)
}
G=colMeans(Weight)
if(G[2]>G[1]){
Weight=Weight[,c(2,1)]
Vote=Vote[,c(2,1)]
}
return(list(Weight=Weight,Cluster=Vote))
}

cluster_voting_uv=function(by,bX,cluster.index,theta1,theta2,sigma1,sigma2,main.cluster.thres=0.49){
res1=as.vector(by-bX*theta1)
res2=as.vector(by-bX*theta2)
J=max(cluster.index)
Weight=matrix(0,length(by),2)
Vote=matrix(0,length(by),2)
for(j in 1:J){
indj=which(cluster.index==j)
logL1 <- sum(dnorm(res1[indj], mean = 0, sd = sqrt(sigma1), log = TRUE))
logL2 <- sum(dnorm(res2[indj], mean = 0, sd = sqrt(sigma2), log = TRUE))
C = max(logL1, logL2)
logL1_adj <- logL1 - C
logL2_adj <- logL2 - C
L1_adj <- exp(logL1_adj)
L2_adj <- exp(logL2_adj)
w1 <- L1_adj / (L1_adj + L2_adj)
w2 <- L2_adj / (L1_adj + L2_adj)
Weight[indj,1]=w1
Weight[indj,2]=w2
Vote[indj,1]=ifelse(w1>main.cluster.thres,1,0)
Vote[indj,2]=ifelse(w2>main.cluster.thres,1,0)
}
G=colMeans(Weight)
if(G[2]>G[1]){
Weight=Weight[,c(2,1)]
Vote=Vote[,c(2,1)]
}
return(list(Weight=Weight,Cluster=Vote))
}

auto_cluster = function(corr_mat, size) {

p <- ncol(corr_mat)
dist_mat <- as.dist(1 - abs(corr_mat))
hc <- hclust(dist_mat, method = "ward.D2")
num_clusters <- round(p / size)
cluster_indices <- cutree(hc, k = num_clusters)

return(cluster_indices)
}
