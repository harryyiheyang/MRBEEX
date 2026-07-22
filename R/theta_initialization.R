theta_mcp_initialization <- function(XtX,Xty,yty,n,dfmax){
XtX=as.matrix(XtX)
Xty=as.numeric(Xty)
p=ncol(XtX)
if(nrow(XtX)!=p||length(Xty)!=p){
stop("XtX and Xty have incompatible dimensions for theta initialization.")
}
if(any(!is.finite(XtX))||any(!is.finite(Xty))||length(yty)!=1||!is.finite(yty)){
stop("Theta initialization inputs must be finite.")
}
n=as.integer(n)[1]
if(!is.finite(n)||n<p+2L){
stop("Theta initialization requires at least p + 2 IVs.")
}
dfmax=min(as.integer(dfmax)[1],p)
if(!is.finite(dfmax)||dfmax<1L){
stop("dfmax must allow at least one nonzero theta coefficient.")
}

XtX=XtX/2+t(XtX)/2
fit.eigen=CppMatrix::matrixEigen(XtX)
evals=as.numeric(fit.eigen$values)
escale=max(abs(evals),1)
evals=pmax(evals,sqrt(.Machine$double.eps)*escale)
XtX=CppMatrix::matrixListProduct(list(
fit.eigen$vectors,diag(evals,nrow=p),t(fit.eigen$vectors)
))
XtX=XtX/2+t(XtX)/2

theta.ls=c(CppMatrix::matrixSolve(XtX,Xty))
explained=sum(Xty*theta.ls)
residual=max(yty-explained,sqrt(.Machine$double.eps)*max(abs(yty),abs(explained),1))

A=matrix(0,n,p+2L)
A[,1]=1
A[cbind(seq_len(p+1L),2:(p+2L))]=1
Q=qr.Q(qr(A))[,2:(p+2L),drop=FALSE]
Xhalf=t(fit.eigen$vectors)*sqrt(evals)
X0=CppMatrix::matrixMultiply(Q[,seq_len(p),drop=FALSE],Xhalf)
y0=CppMatrix::matrixVectorMultiply(X0,theta.ls)+Q[,p+1L]*sqrt(residual)

fit=ncvreg::ncvreg(X0,y0,penalty="MCP",gamma=3,alpha=1,
                   nlambda=100,dfmax=dfmax,returnX=FALSE)
k.mid=(2+log(n))/2
ic=stats::AIC(fit,k=k.mid)
B=fit$beta[-1,,drop=FALSE]
nz=colSums(abs(B)>sqrt(.Machine$double.eps))
j=if(any(is.finite(ic))) which.min(replace(ic,!is.finite(ic),Inf)) else integer(0)

if(length(j)==0L||nz[j]==0L){
fit=ncvreg::ncvreg(X0,y0,penalty="lasso",alpha=0.5,
                   nlambda=100,dfmax=dfmax,returnX=FALSE)
ic=stats::AIC(fit,k=k.mid)
B=fit$beta[-1,,drop=FALSE]
nz=colSums(abs(B)>sqrt(.Machine$double.eps))
j=if(any(is.finite(ic))) which.min(replace(ic,!is.finite(ic),Inf)) else integer(0)
if(length(j)==0L||nz[j]==0L){
candidates=which(nz>0L)
if(length(candidates)==0L){
stop(paste0("Initialization failed: the exposure signal is too weak ",
            "to identify a nonzero theta support."))
}
finite.candidates=candidates[is.finite(ic[candidates])]
if(length(finite.candidates)>0L){
j=finite.candidates[which.min(ic[finite.candidates])]
}else{
j=candidates[1]
}
}
}

support=which(abs(B[,j])>sqrt(.Machine$double.eps))
theta=rep(0,p)
theta[support]=c(CppMatrix::matrixSolve(
XtX[support,support,drop=FALSE],Xty[support]
))
if(any(!is.finite(theta))){
stop("Initialization failed: the selected theta refit is not finite.")
}
names(theta)=colnames(XtX)
return(theta)
}
