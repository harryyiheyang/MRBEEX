MRFit=function(by,bX,theta1,theta2,cluster1,cluster2,df1,df2){
m=length(by)
res1=as.vector(by-bX%*%theta1)
res2=as.vector(by-bX%*%theta2)
m1=length(cluster1)
m2=length(cluster2)
rss1=sum(res1[cluster1]^2)/(m1-df1)
rss2=sum(res2[cluster2]^2)/(m2-df2)
rss2=max(rss2,0.25)
rss1=max(rss1,0.25)
alpha1=m1/m
alpha2=m2/m
en=2*log(alpha2)*m2
return(log(rss1)*m1+log(rss2)*m2-en)
}

SampleOver=function(Rnn,Nxy){
p=length(Nxy)
S=diag(p)
for(i in 1:p){
for(j in 1:i){
S[i,j]=S[j,i]=Rnn[i,j]*min(Nxy[i],Nxy[j])/sqrt(Nxy[i])/sqrt(Nxy[j])
}
}
return(S)
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
positiveinv=function(A,min.eps=1e-8){
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
J=max(cluster.index)
z=x
a_const=a / (a - 1)
x_split=split(x, cluster.index)
z_split=lapply(x_split, function(xj) {
delta=sqrt(mean(xj^2))
zj=max(1 - lambda / delta, 0) * xj
zj=a_const * zj
if (delta > (a * lambda)) {
zj=xj
}
return(zj)
})
z=unsplit(z_split, cluster.index)
return(z)
}

trace=function(A){
a=sum(diag(A))
return(a)
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

IVweight=function(byse,bXse,Rxy,byseinv=NULL,LDSC=NULL,Omega=NULL){
bZse=cbind(bXse,byse)
p=dim(bZse)[2]
n=dim(bZse)[1]
RxyList=array(0,c(n,p,p))
if(is.null(LDSC)==T){
for(i in 1:n){
s=bZse[i,]
RxyList[i,,]=t(t(Rxy)*s)*s
}
}else{
for(i in 1:n){
s=bZse[i,]
A=t(t(Rxy)*s)*s+LDSC[i]*Omega*byseinv[i]^2
RxyList[i,,]=A
}
}
return(RxyList)
}

biasterm=function(RxyList,indvalid,Weight=NULL){
if(is.null(Weight)==1){
Weight=rep(1,dim(RxyList)[1])
}
X=RxyList[1,,]*0
n=dim(RxyList)[1]
for(i in indvalid){
X=X+RxyList[i,,]*Weight[i]
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

block_subsampling=function(cluster_result, rho) {
unique_clusters=unique(cluster_result)
num_selected_clusters=ceiling(length(unique_clusters) * rho)
selected_clusters=sample(unique_clusters, num_selected_clusters)
selected_indices=which(cluster_result %in% selected_clusters)
return(selected_indices)
}

colSD=function(A){
a=A[1,]
for(i in 1:ncol(A)){
a[i]=sd(A[,i])
}
return(a)
}

colSDMAD <- function(A, probs = c(0.015, 0.985), consistancy_correction = TRUE){
winsor_sd_one <- function(x) {
lims <- quantile(x, probs = probs, na.rm = TRUE)
x_win <- pmin(pmax(x, lims[1]), lims[2])
sd_raw <- sd(x_win, na.rm = TRUE)
return(sd_raw)
}
sds <- apply(A, 2, winsor_sd_one)
if (consistancy_correction) {
sim_norm <- rnorm(100000)
lims_sim <- quantile(sim_norm, probs = probs)
sim_win <- pmin(pmax(sim_norm, lims_sim[1]), lims_sim[2])
factor <- 1 / sd(sim_win)
sds <- sds * factor
}
return(sds)
}

robust_sd <- function(x, probs = c(0.015, 0.985), consistancy_correction = TRUE) {
lims   <- stats::quantile(x, probs = probs, na.rm = TRUE)
x_win <- pmin(pmax(x, lims[1]), lims[2])
sd_raw <- stats::sd(x_win, na.rm = TRUE)
if (consistancy_correction) {
sim_norm <- stats::rnorm(100000)
lims_sim <- stats::quantile(sim_norm, probs = probs)
sim_win  <- pmin(pmax(sim_norm, lims_sim[1]), lims_sim[2])
factor   <- 1 / stats::sd(sim_win)
sd_raw   <- sd_raw * factor
}
return(sd_raw)
}

covmad=function(A){
R=cor(A,method="spearman")
R=2*sin(R/6*pi)
R[is.na(R)]=0
D=colSDMAD(A)
S=t(t(R)*D)*D
return(S)
}

subsample_indices=function(cluster.index, sampling.ratio) {
unique_clusters=unique(cluster.index)
n_subsample=round(length(unique_clusters) * sampling.ratio)
subsampled_clusters=sample(unique_clusters, n_subsample)
subsample_indices=which(cluster.index %in% subsampled_clusters)
return(subsample_indices)
}

generate_default_cluster=function(n) {
max_category=ceiling(n / 2)
cluster_index=rep(1:max_category, each = 2, length.out = n)
return(cluster_index)
}

bimin <- function(mat) {
min_element <- min(mat, na.rm = TRUE)
min_indices <- which(mat == min_element, arr.ind = TRUE)

if (nrow(min_indices) > 1) {
min_indices <- min_indices[order(min_indices[,1], min_indices[,2]), ]
}

return(min_indices[1, , drop = FALSE])
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
total.var=mean(bx*(Theta%*%bx))
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

validadj=function(vector1, vector2, tau){
diff=length(vector2) / length(vector1)
if (diff < tau) {
missing_indices=setdiff(1:length(vector1), vector2)
sorted_missing_indices=missing_indices[order(vector1[missing_indices])]
num_to_add=ceiling(tau * length(vector1)) - length(vector2)
vector2=c(vector2, sorted_missing_indices[1:num_to_add])
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

cluster_voting=function(by,bX,cluster.index,theta1,theta2,sigma1,sigma2,m1,m2,main.cluster.thres){
res1=as.vector(by-bX%*%theta1)
res2=as.vector(by-bX%*%theta2)
J=max(cluster.index)
Weight=matrix(0,length(by),2)
Vote=matrix(0,length(by),2)
ratio1=m1/length(by)
ratio2=m2/length(by)
for(j in 1:J){
indj=which(cluster.index==j)
logL1=sum(dnorm(res1[indj], mean = 0, sd = sqrt(sigma1), log = TRUE))
logL2=sum(dnorm(res2[indj], mean = 0, sd = sqrt(sigma2), log = TRUE))
C = max(logL1, logL2)
logL1_adj=logL1 - C
logL2_adj=logL2 - C
L1_adj=exp(logL1_adj)*ratio1
L2_adj=exp(logL2_adj)*ratio2
w1=L1_adj / (L1_adj + L2_adj)
w2=L2_adj / (L1_adj + L2_adj)
w1=w1
Weight[indj,1]=w1
Weight[indj,2]=w2
Vote[indj,1]=ifelse(w1>main.cluster.thres,1,0)
Vote[indj,2]=ifelse(w1<=main.cluster.thres,1,0)
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
logL1=sum(dnorm(res1[indj], mean = 0, sd = sqrt(sigma1), log = TRUE))
logL2=sum(dnorm(res2[indj], mean = 0, sd = sqrt(sigma2), log = TRUE))
C = max(logL1, logL2)
logL1_adj=logL1 - C
logL2_adj=logL2 - C
L1_adj=exp(logL1_adj)
L2_adj=exp(logL2_adj)
w1=L1_adj / (L1_adj + L2_adj)
w2=L2_adj / (L1_adj + L2_adj)
Weight[indj,1]=w1
Weight[indj,2]=w2
Vote[indj,1]=ifelse(w1>main.cluster.thres,1,0)
Vote[indj,2]=ifelse(w2<=main.cluster.thres,1,0)
}
G=colMeans(Weight)
if(G[2]>G[1]){
Weight=Weight[,c(2,1)]
Vote=Vote[,c(2,1)]
}
return(list(Weight=Weight,Cluster=Vote))
}

auto_cluster = function(corr_mat, size) {
p=ncol(corr_mat)
dist_mat=as.dist(1 - abs(corr_mat))
hc=hclust(dist_mat, method = "ward.D2")
num_clusters=round(p / size)
cluster_indices=cutree(hc, k = num_clusters)
return(cluster_indices)
}

generate_D_matrix <- function(s, sign_vec) {
p <- length(s)
if (length(sign_vec) != p) {
stop("Length of sign_vec must match length of s.")
}

if (p == 1) {
D <- 0
} else {
num_pairs <- p*(p-1)/2
D_all <- matrix(0, nrow = num_pairs, ncol = p)

row_idx <- 1
for (i in 1:(p-1)) {
for (j in (i+1):p) {
D_all[row_idx, i] <-  sign_vec[i] / s[i]
D_all[row_idx, j] <- -sign_vec[j] / s[j]
row_idx <- row_idx + 1
}
}

D <- t(D_all)%*%D_all
}

return(D)
}

generate_block_matrix <- function(vars_df, s, theta) {
ind=which(theta==0)
vars_df$cs[which(vars_df$variable%in%ind)]=-1
concerned_vars <- vars_df[vars_df$cs != -1, ]
cs_values <- unique(concerned_vars$cs)
max_var_index <- max(vars_df$variable)
final_matrix <- matrix(0, nrow = max_var_index, ncol = max_var_index)
for (cs_val in cs_values) {
group_vars <- concerned_vars$variable[concerned_vars$cs == cs_val]
if (length(group_vars) > 1) {
group_s <- s[group_vars]
D <- generate_D_matrix(group_s,sign(theta[group_vars]))
final_matrix[group_vars, group_vars]=D
}
}
return(final_matrix)
}

generate_group_matrix <- function(group_index,COV) {
J=max(group_index)
p=length(COV[1,])
G=diag(p)*0
for(j in 1:J){
indj=which(group_index==j)
s=diag(COV)[indj]
s=median(s)/s
theta=COV[indj[1],indj]
G[indj,indj]=generate_D_matrix(s,sign(theta))
}
return(G)
}

generate_block_matrix_CARMA <- function(sumstat.result, s, theta) {
ind=which(theta==0)
concerned_vars <- sumstat.result[sumstat.result$cs >0, ]
cs_values <- unique(concerned_vars$cs)
max_var_index <- nrow(sumstat.result)
final_matrix <- matrix(0, nrow = max_var_index, ncol = max_var_index)
for (cs_val in cs_values) {
group_vars <- concerned_vars$variable[concerned_vars$cs == cs_val]
if (length(group_vars) > 1) {
group_s <- s[group_vars]
D <- generate_D_matrix(group_s,sign(theta[group_vars]))
final_matrix[group_vars, group_vars]=D
}
}
return(final_matrix)
}

group.pip.filter=function(pip.summary,xQTL.cred.thres=0.95,xQTL.pip.thres=0.1){
ind=which(pip.summary$cs>0)
if(length(ind)>0){
J=max(pip.summary$cs[ind])
pip.summary$cs.pip=pip.summary$variable_prob
for(i in 1:J){
indi=which(pip.summary$cs==i)
summaryi=pip.summary[indi,]
pip.cred=sum(summaryi$variable_prob)
pip.summary$cs.pip[indi]=pip.cred
}
ind.keep=which(pip.summary$cs.pip>=xQTL.cred.thres&pip.summary$variable_prob>=xQTL.pip.thres)
cs=pip.summary$cs
cs.pip=pip.summary$cs.pip
cs->cs[pip.summary$variable]
cs.pip->cs.pip[pip.summary$variable]
cs[which(cs==-1)]=0
}else{
ind.keep=NULL
cs=pip.summary$cs.pip*0
cs.pip=pip.summary$cs.pip*0
}
return(list(ind.keep=pip.summary$variable[ind.keep],cs=cs,cs.pip=cs.pip,result=pip.summary))
}

inf.test=function(res.inf,LD,LD2,Theta,A,var.res=1){
if(is.matrix(A)==T){
Varinf=LD
Var2=LD2
V=Theta/var.res
AT=matrixMultiply(t(A),V)
P=V-matrixMultiply(t(AT),matrixMultiply(CppMatrix::matrixGeneralizedInverse(matrixMultiply(AT,A)),AT))
u=sum(res.inf^2)/2
PVar2=matrixMultiply(P,Var2)
e=sum(diag(PVar2))/2
h=sum(diag(matrixMultiply(PVar2,PVar2)))/2
kappa=h/2/e
v=2*e^2/h
}
if(is.vector(A)==T){
A=c(A)
Var2=LD2
V=Theta/var.res
delta=1/sum(A*matrixVectorMultiply(V,A))
AAT=matrixMultiply(t(t(A)),t(A))/delta
P=V-matrixListProduct(list(V,AAT,V))
u=sum(res.inf^2)/2/var.res
PVar2=matrixMultiply(P,Var2)
e=sum(diag(PVar2))/2
h=sum(diag(matrixMultiply(PVar2,PVar2)))/2
kappa=h/2/e
v=2*e^2/h
}
if(is.null(A)==1){
Varinf=LD
Var2=LD2
V=Theta/var.res
u=sum(res.inf^2)/2/var.res
VVar2=LD
e=sum(diag(VVar2))/2
h=sum(diag(LD2))/2
kappa=h/2/e
v=2*e^2/h
}
pv=pchisq(u/kappa,v,lower.tail=F)
return(pv)
}

findUniqueNonZeroRows <- function(M) {
nonZeroCounts <- colSums(M != 0)
uniqueCols <- which(nonZeroCounts == 1)

if(length(uniqueCols) == 0) {
uniqueRows=NULL
}
if(length(uniqueCols)==1){
uniqueRows=which(M[,uniqueCols]!=0)
}
if(length(uniqueCols)>1){
nonZeroCounts=rowSums(M[,uniqueCols]!=0)
uniqueRows <- unique(which(nonZeroCounts>0))
}
return(uniqueRows)
}

center.classifying=function(theta,theta.source){
p=length(theta)
complement=cluster=theta*0
for(i in 1:p){
s=which.min(c(abs(theta[i]),abs(theta[i]-theta.source[i])))
complement[i]=ifelse(s==1,0,theta.source[i])
cluster[i]=ifelse(s==1,1,2)
}
return(list(complement=complement,cluster=cluster))
}

top_K_pip=function(susie_summary,top_K=1,pip.min.thres=0.01,xQTL.pip.thres=0.5){
ind=which(susie_summary$cs>0&susie_summary$variable_prob>=pip.min.thres)
if(length(ind)>0){
susie_summary=susie_summary[ind,]
J=max(susie_summary$cs)
index=c()
for(j in 1:J){
indj=which(susie_summary$cs==j)
g=susie_summary[indj,]
if(length(indj)<=top_K){
index=c(index,g$variable)
}
if(length(indj)>top_K){
index=c(index,g$variable[top_K_indices(g$variable_prob,k=top_K)])
}
}
}
if(length(ind)==0){
index=which(susie_summary$variable_prob>=xQTL.pip.thres)
}
return(index)
}

top_K_indices <- function(vec, k=1) {
return(order(vec, decreasing = TRUE)[1:k])
}

fix_empty_resamples = function(bXj, bX0j, dBtB) {
bXj=as.matrix(bXj)
bX0j=as.matrix(bX0j)
stopifnot(all(dim(bXj) == dim(bX0j)))
sampling = nrow(bXj)
p = ncol(bXj)
zero_rows = which(rowSums(abs(bXj)) == 0 & rowSums(abs(bX0j)) == 0)
if (length(zero_rows) > 0) {
for (i in zero_rows) {
sd = sqrt(dBtB[i])
bXj[i, ] = rnorm(p, mean = 0, sd = sd)
bX0j[i, ] = rnorm(p, mean = 0, sd = sd)
}
}
return(list(bXj = bXj, bX0j = bX0j))
}

susie_effect_resampling = function(LD, alpha, mu, mu2) {
L = nrow(mu)
p = ncol(mu)
# Step 0: Compute variance and standard deviation matrix
varr = mu2 - mu^2
varr[varr < 0] = 0
sd_mat = sqrt(varr)
# Step 1: Compute weighted means and standard deviations
mu_sel = numeric(L)
sd_sel = numeric(L)
for (j in 1:L) {
mu_sel[j] = sum(alpha[j, ] * mu[j, ])
sd_sel[j] = sqrt(sum(alpha[j, ] * sd_mat[j, ]^2))
}
# Step 2: Sample one effect size per component
noise = rnorm(L)
beta_sample = mu_sel + noise * sd_sel  # length L
beta_matrix=matrixVectorMultiply(t(alpha),beta_sample)
# Step 4: Compute LD-propagated signal
bx_mean = matrixVectorMultiply(LD, beta_matrix)
a_mean = beta_matrix
return(list(bx = bx_mean, bx0 = a_mean))
}


spearmancov=function(A){
p=dim(A)[2]
s=c(1:p)
for(i in 1:p){
s[i]=median(abs(median(A[,i])-A[,i]))*1.483
}
R=cor(A,method="spearman")
R=2*sin(R*pi/6)
S=t(R*s)*s
return(S)
}

last_min <- function(x, na.rm = FALSE) {
if (length(x) == 0) return(NA_integer_)
min_val <- min(x, na.rm = na.rm)
idx <- which(x == min_val)
if (length(idx) == 0) return(NA_integer_) # all NA case
return(max(idx))
}

#' @export
cluster_prob <- function(cluster.index, R, alpha = 0, group_size = 4) {
ids <- sort(unique(cluster.index))
w <- sapply(ids, function(i) {
idx <- which(cluster.index == i)
Ri  <- R[idx, idx, drop = FALSE]
ev  <- eigen(Ri, only.values = TRUE)$values
sr  <- sum(ev)^2 / sum(ev^2)
if (!is.finite(sr) || sr <= 0) sr <- 1e-8
sr
})
p0 <- w / sum(w)
if (alpha <= 0) {
names(p0) <- ids
return(p0)
}
ord       <- order(p0)
p_sorted  <- p0[ord]
n         <- length(p_sorted)
group_id  <- ceiling(seq_len(n) / group_size)
med <- tapply(p_sorted, group_id, median)
alpha <- max(min(alpha, 1), 0)
p_smooth_sorted <- (1 - alpha) * p_sorted + alpha * med[group_id]
p_smooth <- p_smooth_sorted[order(ord)]
p_smooth <- p_smooth / sum(p_smooth)
names(p_smooth) <- ids
return(p_smooth)
}

xtx_positive <- function(XtX, Xty, eps = 1e-10) {
d <- diag(XtX)
ind <- which(d <= eps | is.na(d) | !is.finite(d))
if (length(ind) > 0) {
for (i in ind) {
XtX[i, ] <- 0
XtX[, i] <- 0
XtX[i, i] <- 1
Xty[i] <- 0
}
}
return(list(XtX = XtX, Xty = Xty))
}

precompute_cluster_blocks <- function(bX, bXse, by, byse, LD, Theta, Thetarho, cluster.index) {
unique_clusters <- sort(unique(cluster.index))
n_clusters <- length(unique_clusters)
cluster_cache <- vector("list", n_clusters)
names(cluster_cache) <- as.character(unique_clusters)
for (i in 1:n_clusters) {
g <- unique_clusters[i]
idx <- which(cluster.index == g)
block_size <- length(idx)
LDg <- LD[idx, idx, drop = FALSE]
Thetag <- Theta[idx, idx, drop = FALSE]
Thetarhog <- Thetarho[idx, idx, drop = FALSE]
bXg <- bX[idx, , drop = FALSE]
bXseg <- bXse[idx, , drop = FALSE]
byg <- by[idx]
byseg <- byse[idx]
Btg <- as.matrix(t(bXg) %*% Thetag)
BtBg <- matrixMultiply(Btg, bXg)
BtBg <- (t(BtBg) + BtBg) / 2
dBtBg <- diag(BtBg)
cluster_cache[[i]] <- list(
idx = idx,
block_size = block_size,
LD = LDg,
Theta = Thetag,
Thetarho = Thetarhog,
Bt = Btg,
BtB = BtBg,
dBtB = dBtBg,
bX = bXg,
bXse = bXseg,
by = byg,
byse = byseg
)
}
return(cluster_cache)
}

precompute_cluster_blocks_mixture <- function(bX, bXse, by, byse, TC, cluster.index) {
unique_clusters <- sort(unique(cluster.index))
n_clusters <- length(unique_clusters)
cluster_cache <- vector("list", n_clusters)
names(cluster_cache) <- as.character(unique_clusters)
for (i in 1:n_clusters) {
g <- unique_clusters[i]
idx <- which(cluster.index == g)
TCg <- TC[idx, idx, drop = FALSE]
bXg <- bX[idx, , drop = FALSE]
byg <- by[idx]
tilde.Xg <- as.matrix(TCg %*% bXg)
tilde.yg <- as.vector(TCg %*% byg)
cluster_cache[[i]] <- list(
idx = idx,
tilde.X = tilde.Xg,
tilde.y = tilde.yg,
bXse = bXse[idx, , drop = FALSE],
byse = byse[idx]
)
}
return(cluster_cache)
}

precompute_cluster_blocks_uv <- function(bX, bXse, by, byse, LD, Theta, Thetarho, cluster.index, rho) {
unique_clusters <- sort(unique(cluster.index))
n_clusters <- length(unique_clusters)
cluster_cache <- vector("list", n_clusters)
names(cluster_cache) <- as.character(unique_clusters)
for (i in 1:n_clusters) {
g <- unique_clusters[i]
idx <- which(cluster.index == g)
LDg <- LD[idx, idx, drop = FALSE]
Thetag <- Theta[idx, idx, drop = FALSE]
Thetarhog <- Thetarho[idx, idx, drop = FALSE]
bXg <- bX[idx]
bXseg <- bXse[idx]
byg <- by[idx]
byseg <- byse[idx]
Btg <- as.vector(Thetag %*% bXg)
BtBg <- sum(bXg * Btg)
cluster_cache[[i]] <- list(
idx = idx,
LD = LDg,
Theta = Thetag,
Thetarho = Thetarhog,
Bt = Btg,
BtB = BtBg,
bX = bXg,
bXse = bXseg,
by = byg,
byse = byseg
)
}

return(cluster_cache)
}
