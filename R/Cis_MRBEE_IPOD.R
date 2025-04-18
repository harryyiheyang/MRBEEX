#' @export
Cis_MRBEE_IPOD=function(by,bX,byse,bXse,LD="identity",Rxy,tauvec=seq(3,50,by=5),max.iter=100,max.eps=0.001,ebic.gamma=1,reliability.thres=0.8,rho=2,maxdiff=1.5,theta.ini=F,gamma.ini=F,ebic.theta=1,pleiotropy.rm=NULL){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
LD=as.matrix(LD)
Theta=solve(LD)
TC=chol(Theta)
bXinv=as.matrix(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
Bt=t(bXinv)
BtB=matrixMultiply(Bt,bX)
Thetarho=solve(LD+rho*diag(m))
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
if(theta.ini[1]==F){
if(length(tilde.y)<2000){
RC=as.matrix(TC%*%LD)
fit0=varbvs(X=RC,Z=tilde.X,y=tilde.y,verbose=F,maxiter=100)
gamma.ini=gamma.ini1=fit0$beta*(fit0$pip>0.5)
theta.ini=theta.ini1=fit0$beta.cov[-1]
}else{
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak",pv.thres=0.01)
gamma.ini=gamma.ini1=fit0$gamma
theta.ini=theta.ini1=fit0$theta
}
}else{
gamma.ini=gamma.ini/byse
theta.ini=theta.ini
}
############################## Tuning Parameter ######################
w=length(tauvec)
Btheta=array(0,c(p,w))
Bgamma=array(0,c(m,w))
Bbic=c(1:w)
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
Hinv=matrixInverse(BtB-Rxysum[1:p,1:p])
g=matrixVectorMultiply(Bt,as.vector(by-LD%*%gamma))-Rxysum[1:p,p+1]
theta=c(matrixVectorMultiply(Hinv,g))
if((norm(theta,"2")/norm(theta.ini,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini,"2")
}
########################### update gamma ############################
gamma=as.vector(Thetarho%*%(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1))
gamma[pleiotropy.rm]=0
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
gamma1[pleiotropy.rm]=0
delta=delta+rho*(gamma-gamma1)
delta[pleiotropy.rm]=0
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta[,j]=theta
Bgamma[,j]=gamma1
df1=sum(gamma1!=0)
res=c(by-matrixVectorMultiply(bX,theta)-as.vector(LD%*%gamma))
rss=sum(res*(Theta%*%res))/(m-df1-p)
Bbic[j]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+p*(log(m)+ebic.theta*log(p))
}
Bbic=Bbic/m
######################## Inference #################################
jstar=which.min(Bbic)
theta=Btheta[,jstar]
gamma=Bgamma[,jstar]
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma!=0)
adjf=m/(length(indvalid)-p)
indvalid=which(gamma==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res=by-matrixVectorMultiply(bX,theta)-matrixVectorMultiply(LD,gamma)
var_error=sum(res*(Theta%*%res))/(m-p-length(indgamma))
var_error=max(1,var_error)

bZ=cbind(bX,as.matrix(LD[,indgamma]))
TCbZ=as.matrix(TC%*%bZ)
Hinv=t(bZ)%*%(Theta%*%bZ)
Hinv[1:p,1:p]=Hinv[1:p,1:p]-Rxysum[1:p,1:p]
Hinv=as.matrix(Hinv)
Hinv=positiveinv(as.matrix(Hinv))
Hinv1=var_error*matrixListProduct(list(t(bZ),Theta,bZ))
COV=Hinv%*%Hinv1%*%Hinv
theta.cov=COV[1:p,1:p]
theta.se=sqrt(diag(theta.cov))

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$theta.z=A$theta/A$theta.se
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini*byse1
A$reliability.adjust=r
A$theta.pratt=getPratt(bX=bX,by=by,bXse=bXse,byse=byse,Theta=Theta,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma,Theta=Theta,LD=LD)
A$var.error=var_error
return(A)
}
