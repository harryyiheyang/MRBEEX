MRBEE_Mixture_SuSiE=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),main.cluster.thres=0.45,min.cluster.size=5,Lvec=c(1:min(5,ncol(bX))),pip.thres=0.2,ebic.theta=1,reliability.thres=0.8,sampling.time=100,max.iter=30,max.eps=5e-4,sampling.iter=5,susie.iter=100,ridge.diff=1e5,projection.eigen.floor=1,verbose=T,pip.min=0.1,cred.pip.thres=0.95,estimate_residual_variance=T,group.penalize=F,group.index=c(1:ncol(bX)[1]),group.diff=10,coverage.causal=0.95,LDSC=NULL,Omega=NULL,estimate_residual_method="MoM",standardize=T,tau=5,step.size=0.5,theta.ini.1=NULL,theta.ini.2=NULL){
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
  r=reliability.adj(bX,bXse%*%diag(sqrt(diag(Rxy[1:p,1:p]))),Theta=Theta,thres=reliability.thres)
  r=c(r,1)
  Rxy=t(t(Rxy)*r)*r
  RxyList=IVweight(byse,bXse,Rxy,byseinv=byseinv,LDSC=LDSC,Omega=Omega)

  if(is.null(theta.ini.1) | is.null(theta.ini.2)){
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
    fit.susie.init1=susie(X=tilde.X[cluster.ini.1,,drop=FALSE],y=tilde.y[cluster.ini.1],L=5,intercept=F)
    fit.susie.init2=susie(X=tilde.X[cluster.ini.2,,drop=FALSE],y=tilde.y[cluster.ini.2],L=5,intercept=F)
    theta.ini.1=coef.susie(fit.susie.init1)[-1]*(fit.susie.init1$pip>0.1)
    theta.ini.2=coef.susie(fit.susie.init2)[-1]*(fit.susie.init2$pip>0.1)
    sigma.ini.1=sum((tilde.y[cluster.ini.1]-tilde.X[cluster.ini.1,,drop=FALSE]%*%theta.ini.1)^2)/(length(cluster.ini.1)-sum(theta.ini.1!=0))
    sigma.ini.2=sum((tilde.y[cluster.ini.2]-tilde.X[cluster.ini.2,,drop=FALSE]%*%theta.ini.2)^2)/(length(cluster.ini.2)-sum(theta.ini.2!=0))
    sigma.ini.2=max(0.25,sigma.ini.2)
    sigma.ini.1=max(0.25,sigma.ini.1)
  }else{
    res1=as.vector(tilde.y-tilde.X%*%theta.ini.1)
    res2=as.vector(tilde.y-tilde.X%*%theta.ini.2)
    mad1=mad(res1)^2
    mad2=mad(res2)^2
    sigma.ini.1=sigma.ini.2=max(0.25,min(mad1,mad2))
    m1=floor(m*0.8);m2=m-m1
  }
  Voting.ini=cluster_voting(by=tilde.y,bX=tilde.X,cluster.index=cluster.index,theta1=theta.ini.1,theta2=theta.ini.2,sigma1=sigma.ini.1,sigma2=sigma.ini.2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
  cluster.ini.2=which(Voting.ini$Cluster[,2]==1)
  cluster.ini.1=which(Voting.ini$Cluster[,1]==1)
  cluster.ratio.ini=c(length(cluster.ini.1),length(cluster.ini.2))/m
  m1=length(cluster.ini.1)
  m2=length(cluster.ini.2)

  eta = matrixVectorMultiply(bX,theta.ini.1)
  eta[cluster.ini.2] = matrixVectorMultiply(bX[cluster.ini.2, , drop = FALSE],theta.ini.2)
  if(isLD){
    initial_res = as.vector(Theta %*% (by - eta))
  } else {
    initial_res = by - eta
  }
  gamma.ini = soft(ifelse(abs(initial_res) > quantile(abs(initial_res), 0.975), initial_res, 0),1) # Shrinkage trick

  t2=Sys.time()
  time_to_print=round(difftime(t2, t1, units = "secs"),3)
  if(verbose==T){
    cat(paste0("Initialization ends: ",time_to_print," secs\n"))
  }

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
      gamma=gamma.ini
      iter=0
      error=1
      fit.susie2=NULL
      project_XtX1 <- new_adj_projector()
      project_XtX2 <- new_adj_projector()
      while(iter<max.iter&error>max.eps){
        theta11=theta1
        theta22=theta2

        if(sum(gamma!=0)){
          if(isLD){
            tilde.res=as.vector(TC%*%(by-LD%*%gamma))
          }else{
            tilde.res=by-gamma
          }
        }else{
          if(isLD){
            tilde.res=tilde.y
          }else{
            tilde.res=by
          }
        }

        Rxysum1=biasterm(RxyList=RxyList,cluster1)
        Cmat1=Rxysum1[1:p,1:p]
        XtX1=matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.X[cluster1,,drop=FALSE],transA=TRUE)
        XtX1=XtX1-Cmat1
        XtX1=t(XtX1)/2+XtX1/2
        Xty1=c(matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.res[cluster1],transA=TRUE))-Rxysum1[1:p,1+p]
        yty1=sum(tilde.res[cluster1]^2)
        adjX1=xtx_positive(XtX1,Xty1)
        XtX1=adjX1$XtX
        Xty1=adjX1$Xty
        XtX1.raw=XtX1
        Diff_matrix1=diag(p)*0
        if(group.penalize==T){
          Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1)
        }
        projection.eigen.floor1=projection.eigen.floor*length(cluster1)/m
        XtX1=project_XtX1(XtX1+Cmat1+Diff_matrix1, Cmat1, cluster1, projection.eigen.floor1)
        fit.susie1=tryCatch({
          susie_ss(XtX=XtX1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[v],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
        },error = function(e) {
          susie_ss(XtX=XtX1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[v],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
        })
        theta1=coef.susie(fit.susie1)[-1]*(fit.susie1$pip>pip.min)
        theta.cs1=group.pip.filter(pip.summary=summary(fit.susie1)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
        pip.alive1=theta.cs1$ind.keep
        if(length(pip.alive1)>0){
          theta1[-pip.alive1]=0
        }else{
          theta1=theta1*0
        }
        if(length(cluster2)>min.cluster.size){
          Rxysum2=biasterm(RxyList=RxyList,cluster2)
          Cmat2=Rxysum2[1:p,1:p]
          XtX2=matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.X[cluster2,,drop=FALSE],transA=TRUE)
          XtX2=XtX2-Cmat2
          XtX2=t(XtX2)/2+XtX2/2
          Xty2=c(matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.res[cluster2],transA=TRUE))-Rxysum2[1:p,1+p]
          yty2=sum(tilde.res[cluster2]^2)
          adjX2=xtx_positive(XtX2,Xty2)
          XtX2=adjX2$XtX
          Xty2=adjX2$Xty
          XtX2.raw=XtX2
          Diff_matrix2=diag(p)*0
          if(group.penalize==T){
            Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2)
          }
          projection.eigen.floor2=projection.eigen.floor*length(cluster2)/m
          XtX2=project_XtX2(XtX2+Cmat2+Diff_matrix2, Cmat2, cluster2, projection.eigen.floor2)
          fit.susie2=tryCatch({
            susie_ss(XtX=XtX2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[l],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
          },error = function(e) {
            susie_ss(XtX=XtX2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[l],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
          })
          theta2=coef.susie(fit.susie2)[-1]*(fit.susie2$pip>pip.min)
          theta.cs2=group.pip.filter(pip.summary=summary(fit.susie2)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
          pip.alive2=theta.cs2$ind.keep
          if(length(pip.alive2)>0){
            theta2[-pip.alive2]=0
          }else{
            theta2=theta2*0
          }
          Diff2=generate_block_matrix(summary(fit.susie2)$vars,length(cluster2)/diag(XtX2),theta2)
        }else{
          theta2=theta1*0
          cluster2=c(1:min.cluster.size)
        }
        indtheta1=which(theta1!=0)
        Diff1=generate_block_matrix(summary(fit.susie1)$vars,length(cluster1)/diag(XtX1),theta1)
        XtX1.reols=XtX1.raw+Diff_matrix1
        if(length(indtheta1)==1){
          xtx1=project_select_xtx(XtX1.reols[indtheta1,indtheta1,drop=FALSE],eigen.floor=projection.eigen.floor1)
          xtx1=xtx1[1,1]
          xty1=Xty1[indtheta1]
          theta1[indtheta1]=xty1/xtx1
        }
        if(length(indtheta1)>1){
          XtX1=project_select_xtx(XtX1.reols[indtheta1,indtheta1,drop=FALSE]+ridge.diff*Diff1[indtheta1,indtheta1,drop=FALSE],eigen.floor=projection.eigen.floor1)
          Xty1=Xty1[indtheta1]
          theta1[indtheta1]=as.vector(CppMatrix::matrixSolve(XtX1,Xty1))
        }
        indtheta2=which(theta2!=0)
        if(length(indtheta2)==1){
          XtX2.reols=XtX2.raw+Diff_matrix2
          xtx2=project_select_xtx(XtX2.reols[indtheta2,indtheta2,drop=FALSE],eigen.floor=projection.eigen.floor2)
          xtx2=xtx2[1,1]
          xty2=Xty2[indtheta2]
          theta2[indtheta2]=xty2/xtx2
        }
        if(length(indtheta2)>1){
          XtX2.reols=XtX2.raw+Diff_matrix2
          XtX2=project_select_xtx(XtX2.reols[indtheta2,indtheta2,drop=FALSE]+ridge.diff*Diff2[indtheta2,indtheta2,drop=FALSE],eigen.floor=projection.eigen.floor2)
          Xty2=Xty2[indtheta2]
          theta2[indtheta2]=as.vector(CppMatrix::matrixSolve(XtX2,Xty2))
        }
        sigma1=sum((tilde.res[cluster1]-tilde.X[cluster1,,drop=FALSE]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
        sigma2=sum((tilde.res[cluster2]-tilde.X[cluster2,,drop=FALSE]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
        sigma2=max(0.25,sigma2)
        sigma1=max(0.25,sigma1)
        Voting=cluster_voting(by=tilde.res,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
        cluster2=which(Voting$Cluster[,2]==1)
        if(length(cluster2)==0){
          cluster2=c(1:min.cluster.size)
        }
        cluster1=which(Voting$Cluster[,1]==1)
        m1=length(cluster1)
        m2=length(cluster2)
        cluster.ratio=c(length(cluster1),length(cluster2))/m
        gamma=update_gamma_mcp(tau, by, bX, theta1, theta2, cluster1, cluster2, gamma, LD, isLD, step = step.size)
        iter=iter+1
        if(iter>3){
          error=max(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
        }
      }
      Btheta1[,v,l]=theta1
      Btheta2[,v,l]=theta2
      df1=min(sum(theta1!=0),Lvec[v])
      df2=min(sum(theta2!=0),Lvec[l])
      Bbic[v,l]=MRFit(by=tilde.res,bX=tilde.X,theta1=theta1,theta2=theta2,cluster1=cluster1,cluster2=cluster2,df1=df1,df2=df2)+(df2+df1)*(log(p*2)*ebic.theta+log(m))+log(m)
      SIG1[v,l]=sigma1
      SIG2[v,l]=sigma2
      M1[v,l]=m1
      M2[v,l]=m2
    }
  }
  Bbic=Bbic/m

  vstar=bimin(Bbic)[1]
  lstar=bimin(Bbic)[2]
  theta1=Btheta1[,vstar,lstar]
  theta2=Btheta2[,vstar,lstar]
  sigma1=SIG1[vstar,lstar]
  sigma2=SIG2[vstar,lstar]
  m1=M1[vstar,lstar]
  m2=M2[vstar,lstar]

  gamma=gamma.ini
  if(sum(gamma!=0)){
    if(isLD){
      tilde.res.final=as.vector(TC%*%(by-LD%*%gamma))
    }else{
      tilde.res.final=by-gamma
    }
  }else{
    tilde.res.final=tilde.y
  }
  Voting=cluster_voting(by=tilde.res.final,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
  cluster2=which(Voting$Cluster[,2]==1)
  cluster1=which(Voting$Cluster[,1]==1)
  m1=length(cluster1)
  m2=length(cluster2)
  iter=0
  error=1
  theta11=theta1*0
  theta22=theta2*0
  fit.susie1=fit.susie2=NULL
  project_XtX1 <- new_adj_projector()
  project_XtX2 <- new_adj_projector()
  while(iter<(max.iter*2)&error>max.eps){
    theta11=theta1
    theta22=theta2

    if(sum(gamma!=0)){
      if(isLD){
        tilde.res=as.vector(TC%*%(by-LD%*%gamma))
      }else{
        tilde.res=by-gamma
      }
    }else{
      if(isLD){
        tilde.res=tilde.y
      }else{
        tilde.res=by
      }
    }

    Rxysum1=biasterm(RxyList=RxyList,cluster1)
    Cmat1=Rxysum1[1:p,1:p]
    XtX1=matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.X[cluster1,,drop=FALSE],transA=TRUE)
    XtX1=XtX1-Cmat1
    XtX1=t(XtX1)/2+XtX1/2
    Xty1=c(matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.res[cluster1],transA=TRUE))-Rxysum1[1:p,1+p]
    yty1=sum(tilde.res[cluster1]^2)
    adjX1=xtx_positive(XtX1,Xty1)
    XtX1=adjX1$XtX
    Xty1=adjX1$Xty
    XtX1.raw=XtX1
    Diff_matrix1=diag(p)*0
    if(group.penalize==T){
      Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1)
    }
    projection.eigen.floor1=projection.eigen.floor*length(cluster1)/m
    XtX1=project_XtX1(XtX1+Cmat1+Diff_matrix1, Cmat1, cluster1, projection.eigen.floor1)
    fit.susie1=tryCatch({
      susie_ss(XtX=XtX1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[vstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
    },error = function(e) {
      susie_ss(XtX=XtX1,Xty=Xty1,yty=yty1,n=length(cluster1),L=Lvec[vstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie1,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
    })
    theta1=coef.susie(fit.susie1)[-1]*(fit.susie1$pip>pip.min)
    theta.cs1=group.pip.filter(pip.summary=summary(fit.susie1)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
    pip.alive1=theta.cs1$ind.keep
    if(length(pip.alive1)>0){
      theta1[-pip.alive1]=0
    }else{
      theta1=theta1*0
    }
    if(length(cluster2)>min.cluster.size){

      Rxysum2=biasterm(RxyList=RxyList,cluster2)
      Cmat2=Rxysum2[1:p,1:p]
      XtX2=matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.X[cluster2,,drop=FALSE],transA=TRUE)-Cmat2
      XtX2=t(XtX2)/2+XtX2/2
      Xty2=c(matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.res[cluster2],transA=TRUE))-Rxysum2[1:p,1+p]
      yty2=sum(tilde.res[cluster2]^2)
      adjX2=xtx_positive(XtX2,Xty2)
      XtX2=adjX2$XtX
      Xty2=adjX2$Xty
      XtX2.raw=XtX2
      Diff_matrix2=diag(p)*0
      if(group.penalize==T){
        Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2)
      }
      projection.eigen.floor2=projection.eigen.floor*length(cluster2)/m
      XtX2=project_XtX2(XtX2+Cmat2+Diff_matrix2, Cmat2, cluster2, projection.eigen.floor2)
      fit.susie2=tryCatch({
        susie_ss(XtX=XtX2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[lstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
      },error = function(e) {
        susie_ss(XtX=XtX2,Xty=Xty2,yty=yty2,n=length(cluster2),L=Lvec[lstar],max_iter=susie.iter,estimate_prior_method="EM",model_init=fit.susie2,estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
      })
      theta2=coef.susie(fit.susie2)[-1]*(fit.susie2$pip>pip.min)
      theta.cs2=group.pip.filter(pip.summary=summary(fit.susie2)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
      pip.alive2=theta.cs2$ind.keep
      if(length(pip.alive2)>0){
        theta2[-pip.alive2]=0
      }else{
        theta2=theta2*0
      }
      Diff2=generate_block_matrix(summary(fit.susie2)$vars,length(cluster2)/diag(XtX2),theta2)
    }else{
      theta2=theta1*0
      cluster2=c(1:4)
    }
    indtheta1=which(theta1!=0)
    Diff1=generate_block_matrix(summary(fit.susie1)$vars,length(cluster1)/diag(XtX1),theta1)
    XtX1.reols=XtX1.raw+Diff_matrix1
    if(length(indtheta1)==1){
      xtx1=project_select_xtx(XtX1.reols[indtheta1,indtheta1,drop=FALSE],eigen.floor=projection.eigen.floor1)
      xtx1=xtx1[1,1]
      xty1=Xty1[indtheta1]
      theta1[indtheta1]=xty1/xtx1
    }
    if(length(indtheta1)>1){
      XtX1=project_select_xtx(XtX1.reols[indtheta1,indtheta1,drop=FALSE]+ridge.diff*Diff1[indtheta1,indtheta1,drop=FALSE],eigen.floor=projection.eigen.floor1)
      Xty1=Xty1[indtheta1]
      theta1[indtheta1]=c(CppMatrix::matrixSolve(XtX1,Xty1))
    }
    indtheta2=which(theta2!=0)
    if(length(indtheta2)==1){
      XtX2.reols=XtX2.raw+Diff_matrix2
      xtx2=project_select_xtx(XtX2.reols[indtheta2,indtheta2,drop=FALSE],eigen.floor=projection.eigen.floor2)
      xtx2=xtx2[1,1]
      xty2=Xty2[indtheta2]
      theta2[indtheta2]=xty2/xtx2
    }
    if(length(indtheta2)>1){
      XtX2.reols=XtX2.raw+Diff_matrix2
      XtX2=project_select_xtx(XtX2.reols[indtheta2,indtheta2,drop=FALSE]+ridge.diff*Diff2[indtheta2,indtheta2,drop=FALSE],eigen.floor=projection.eigen.floor2)
      Xty2=Xty2[indtheta2]
      theta2[indtheta2]=c(CppMatrix::matrixSolve(XtX2,Xty2))
    }
    sigma1=sum((tilde.res[cluster1]-tilde.X[cluster1,,drop=FALSE]%*%theta1)^2)/(length(cluster1)-sum(theta1!=0))
    sigma2=sum((tilde.res[cluster2]-tilde.X[cluster2,,drop=FALSE]%*%theta2)^2)/(length(cluster2)-sum(theta2!=0))
    sigma2=max(0.25,sigma2)
    sigma1=max(0.25,sigma1)
    Voting=cluster_voting(by=tilde.res,bX=tilde.X,cluster.index=cluster.index,theta1=theta1,theta2=theta2,sigma1=sigma1,sigma2=sigma2,main.cluster.thres=main.cluster.thres,m1=m1,m2=m2)
    cluster2=which(Voting$Cluster[,2]==1)
    if(length(cluster2)==0){
      cluster2=c(1:min.cluster.size)
    }
    cluster1=which(Voting$Cluster[,1]==1)
    m1=length(cluster1)
    m2=length(cluster2)
    cluster.ratio=c(length(cluster1),length(cluster2))/m
    gamma=update_gamma_mcp(tau, by, bX, theta1, theta2, cluster1, cluster2, gamma, LD, isLD, step = step.size)
    iter=iter+1
    if(iter>3){
      error=max(norm(theta1-theta11,"2"),norm(theta2-theta22,"2"))
    }
  }
  t2=Sys.time()
  time_to_print=round(difftime(t2, t1, units = "secs"),3)
  if(verbose==T){
    cat(paste0("Estimation ends: ",time_to_print," secs\n"))
  }
  if(sum(gamma!=0)){
    if(isLD){
      tilde.res.ref=as.vector(TC%*%(by-LD%*%gamma))
    }else{
      tilde.res.ref=by-gamma
    }
  }else{
    if(isLD){
      tilde.res.ref=tilde.y
    }else{
      tilde.res.ref=by
    }
  }
  Rxysum1.ref=biasterm(RxyList=RxyList,cluster1)
  Cmat1.ref=Rxysum1.ref[1:p,1:p]
  XtX1.ref=matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.X[cluster1,,drop=FALSE],transA=TRUE)-Cmat1.ref
  XtX1.ref=XtX1.ref/2+t(XtX1.ref)/2
  Xty1.ref=c(matrixMultiply(tilde.X[cluster1,,drop=FALSE],tilde.res.ref[cluster1],transA=TRUE))-Rxysum1.ref[1:p,1+p]
  adjX1.ref=xtx_positive(XtX1.ref,Xty1.ref)
  XtX1.ref=adjX1.ref$XtX
  Diff_matrix1=diag(p)*0
  if(group.penalize==T){
    Diff_matrix1=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1.ref)
  }
  Veigen1=FProject_basis(XtX1.ref+Diff_matrix1)
  Rxysum2.ref=biasterm(RxyList=RxyList,cluster2)
  Cmat2.ref=Rxysum2.ref[1:p,1:p]
  XtX2.ref=matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.X[cluster2,,drop=FALSE],transA=TRUE)-Cmat2.ref
  XtX2.ref=XtX2.ref/2+t(XtX2.ref)/2
  Xty2.ref=c(matrixMultiply(tilde.X[cluster2,,drop=FALSE],tilde.res.ref[cluster2],transA=TRUE))-Rxysum2.ref[1:p,1+p]
  adjX2.ref=xtx_positive(XtX2.ref,Xty2.ref)
  XtX2.ref=adjX2.ref$XtX
  Diff_matrix2=diag(p)*0
  if(group.penalize==T){
    Diff_matrix2=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2.ref)
  }
  Veigen2=FProject_basis(XtX2.ref+Diff_matrix2)
  t1=Sys.time()
  cat("Resampling starts:\n")
  pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
  names(theta1)=names(theta2)=colnames(bX)
  ThetaList1=ThetaList2=matrix(0,sampling.time,p)
  colnames(ThetaList1)=colnames(ThetaList2)=colnames(bX)
  cluster.index <- as.integer(factor(cluster.index))
  j=1
  consec_error=0
  while(j<=sampling.time){
    indicator <- FALSE
    setTxtProgressBar(pb, j)
      tryCatch({
      indj <- sort(sample.int(m, size = max(1L, floor(0.5 * m)), replace = FALSE))
      bXj=bX[indj,,drop=FALSE]
      bXsej=bXse[indj,,drop=FALSE]
      byj=by[indj]
      bysej=byse[indj]
      LDj=LD[indj,indj]
      if(isLD==T){
        TCj=chol(solve(LDj))
        tilde.Xj=as.matrix(TCj%*%bXj)
        tilde.yj=as.vector(TCj%*%byj)
      }else{
        tilde.Xj=bXj
        tilde.yj=byj
      }
      Rxyallj <- biasterm(RxyList = RxyList, indj)
      theta1j=theta1*runif(p,0.95,1.05)*0.95
      theta2j=theta2*runif(p,0.95,1.05)*0.95
      cluster1j=which(indj%in%cluster1)
      cluster2j=which(indj%in%cluster2)
      m1j=length(cluster1j)
      m2j=length(cluster2j)
      sigma1j=sigma1
      sigma2j=sigma2
      errorj=1
      gammaj=gamma[indj]
      tauj=ifelse(consec_error>10, 1e4, tau)
      fit.susie1j=NULL
      fit.susie2j=NULL
      project_XtX1j <- new_FProjector(Veigen1)
      project_XtX2j <- new_FProjector(Veigen2)

      for(jiter in 1:sampling.iter){

        if(sum(gammaj!=0)){
          if(isLD){
            tilde.resj=as.vector(TCj%*%(byj-LDj%*%gammaj))
          }else{
            tilde.resj=byj-gammaj
          }
        }else{
          if(isLD){
            tilde.resj=tilde.yj
          }else{
            tilde.resj=byj
          }
        }

        theta_prev1j=theta1j
        theta_prev2j=theta2j
        cluster_partitionj <- length(intersect(cluster1j, cluster2j)) == 0 &&
          length(unique(c(cluster1j, cluster2j))) == length(indj)
        if(cluster_partitionj){
          if(length(cluster1j) <= length(cluster2j)){
            Rxysum1j=biasterm(RxyList=RxyList,indj[cluster1j])
            Rxysum2j=Rxyallj-Rxysum1j
          }else{
            Rxysum2j=biasterm(RxyList=RxyList,indj[cluster2j])
            Rxysum1j=Rxyallj-Rxysum2j
          }
        }else{
          Rxysum1j=biasterm(RxyList=RxyList,indj[cluster1j])
        }
        Cmat1j=Rxysum1j[1:p,1:p]
        XtX1j.raw=matrixMultiply(tilde.Xj[cluster1j,,drop=FALSE],tilde.Xj[cluster1j,,drop=FALSE],transA=TRUE)-Cmat1j
        XtX1j.raw=t(XtX1j.raw)/2+XtX1j.raw/2
        Xty1j=c(matrixMultiply(tilde.Xj[cluster1j,,drop=FALSE],tilde.resj[cluster1j],transA=TRUE))-Rxysum1j[1:p,1+p]
        yty1j=sum(tilde.resj[cluster1j]^2)
        adjX1j=xtx_positive(XtX1j.raw,Xty1j)
        XtX1j.raw=adjX1j$XtX
        Xty1j=adjX1j$Xty
        Diff_matrix1j=diag(p)*0
        if(group.penalize==T){
          Diff_matrix1j=group.diff*generate_group_matrix(group_index=group.index,COV=XtX1j.raw)
        }
        projection.eigen.floor1j=projection.eigen.floor*length(cluster1j)/m
        XtX1j.susie=project_XtX1j(XtX1j.raw, Diff_matrix1j, cluster1j, projection.eigen.floor1j)
        fit.susie1j=tryCatch({
          susie_ss(XtX=XtX1j.susie,Xty=Xty1j,yty=yty1j,n=length(cluster1j),L=Lvec[vstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie1j,estimate_prior_method="EM",coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
        },error = function(e) {
          susie_ss(XtX=XtX1j.susie,Xty=Xty1j,yty=yty1j,n=length(cluster1j),L=Lvec[vstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie1j,estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
        })
        theta1j=coef.susie(fit.susie1j)[-1]*(fit.susie1j$pip>pip.min)
        theta.cs1j=group.pip.filter(pip.summary=summary(fit.susie1j)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
        pip.alive1j=theta.cs1j$ind.keep
        if(length(pip.alive1j)>0){
          theta1j[-pip.alive1j]=0
        }else{
          theta1j=theta1j*0
        }
        if(length(cluster2j)>(min.cluster.size/2)){

          if(!cluster_partitionj){
          Rxysum2j=biasterm(RxyList=RxyList,indj[cluster2j])
          }
          Cmat2j=Rxysum2j[1:p,1:p]
          XtX2j.raw=matrixMultiply(tilde.Xj[cluster2j,,drop=FALSE],tilde.Xj[cluster2j,,drop=FALSE],transA=TRUE)-Cmat2j
          XtX2j.raw=XtX2j.raw/2+t(XtX2j.raw)/2
          Xty2j=c(matrixMultiply(tilde.Xj[cluster2j,,drop=FALSE],tilde.resj[cluster2j],transA=TRUE))-Rxysum2j[1:p,1+p]
          yty2j=sum(tilde.resj[cluster2j]^2)
          adjX2j=xtx_positive(XtX2j.raw,Xty2j)
          XtX2j.raw=adjX2j$XtX
          Xty2j=adjX2j$Xty
          Diff_matrix2j=diag(p)*0
          if(group.penalize==T){
            Diff_matrix2j=group.diff*generate_group_matrix(group_index=group.index,COV=XtX2j.raw)
          }
          projection.eigen.floor2j=projection.eigen.floor*length(cluster2j)/m
          XtX2j.susie=project_XtX2j(XtX2j.raw, Diff_matrix2j, cluster2j, projection.eigen.floor2j)
          fit.susie2j=tryCatch({
            susie_ss(XtX=XtX2j.susie,Xty=Xty2j,yty=yty2j,n=length(cluster2j),L=Lvec[lstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie2j,estimate_prior_method="EM",coverage = coverage.causal,estimate_residual_variance=estimate_residual_variance,residual_variance=max(0.9,vary),estimate_residual_method=estimate_residual_method,standardize=standardize)
          },error = function(e) {
            susie_ss(XtX=XtX2j.susie,Xty=Xty2j,yty=yty2j,n=length(cluster2j),L=Lvec[lstar],max_iter=ifelse(jiter==1,1000,30),model_init=fit.susie2j,estimate_prior_method="EM",estimate_residual_variance=F,residual_variance=max(0.9,vary),coverage = coverage.causal,estimate_residual_method=estimate_residual_method,standardize=standardize)
          })
          theta2j=coef.susie(fit.susie2j)[-1]*(fit.susie2j$pip>pip.min)
          theta.cs2j=group.pip.filter(pip.summary=summary(fit.susie2j)$var,xQTL.cred.thres=cred.pip.thres,xQTL.pip.thres=pip.thres)
          pip.alive2j=theta.cs2j$ind.keep
          if(length(pip.alive2j)>0){
            theta2j[-pip.alive2j]=0
          }else{
            theta2j=theta2j*0
          }
          Diff2j=generate_block_matrix(summary(fit.susie2j)$vars,length(cluster2j)/diag(XtX2j.susie),theta2j)
        }else{
          theta2j=theta1j*0
          cluster2j=c(1:2)
        }
        indtheta1j=which(theta1j!=0)
        Diff1j=generate_block_matrix(summary(fit.susie1j)$vars,length(cluster1j)/diag(XtX1j.susie),theta1j)
        XtX1j.reols=XtX1j.raw+Diff_matrix1j
        if(length(indtheta1j)==1){
          xtx1j=project_select_xtx(XtX1j.reols[indtheta1j,indtheta1j,drop=FALSE],eigen.floor=projection.eigen.floor1j)
          xtx1j=xtx1j[1,1]
          xty1j=Xty1j[indtheta1j]
          theta1j[indtheta1j]=xty1j/xtx1j
        }
        if(length(indtheta1j)>1){
          XtX1j=project_select_xtx(XtX1j.reols[indtheta1j,indtheta1j,drop=FALSE]+ridge.diff*Diff1j[indtheta1j,indtheta1j,drop=FALSE],eigen.floor=projection.eigen.floor1j)
          Xty1j=Xty1j[indtheta1j]
          theta1j[indtheta1j]=c(CppMatrix::matrixSolve(XtX1j,Xty1j))
        }
        indtheta2j=which(theta2j!=0)
        if(length(indtheta2j)==1){
          XtX2j.reols=XtX2j.raw+Diff_matrix2j
          xtx2j=project_select_xtx(XtX2j.reols[indtheta2j,indtheta2j,drop=FALSE],eigen.floor=projection.eigen.floor2j)
          xtx2j=xtx2j[1,1]
          xty2j=Xty2j[indtheta2j]
          theta2j[indtheta2j]=xty2j/xtx2j
        }
        if(length(indtheta2j)>1){
          XtX2j.reols=XtX2j.raw+Diff_matrix2j
          XtX2j=project_select_xtx(XtX2j.reols[indtheta2j,indtheta2j,drop=FALSE]+ridge.diff*Diff2j[indtheta2j,indtheta2j,drop=FALSE],eigen.floor=projection.eigen.floor2j)
          Xty2j=Xty2j[indtheta2j]
          theta2j[indtheta2j]=c(CppMatrix::matrixSolve(XtX2j,Xty2j))
        }
        sigma1j=sum((tilde.resj[cluster1j]-tilde.Xj[cluster1j,,drop=FALSE]%*%theta1j)^2)/(length(cluster1j)-sum(theta1j!=0))
        sigma2j=sum((tilde.resj[cluster2j]-tilde.Xj[cluster2j,,drop=FALSE]%*%theta2j)^2)/(length(cluster2j)-sum(theta2j!=0))
        sigma2j=max(0.25,sigma2j)
        sigma1j=max(0.25,sigma1j)
        Votingj=cluster_voting(by=tilde.resj,bX=tilde.Xj,cluster.index=cluster.index[indj],theta1=theta1j,theta2=theta2j,sigma1=sigma1j,sigma2=sigma2j,main.cluster.thres=main.cluster.thres,m1=m1j,m2=m2j)
        cluster2j=which(Votingj$Cluster[,2]==1)
        if(length(cluster2j)==0){
          cluster2j=c(1:2)
        }
        cluster1j=which(Votingj$Cluster[,1]==1)
        m1j=length(cluster1j)
        m2j=length(cluster2j)
        gammaj=update_gamma_mcp(tauj, byj, bXj, theta1j, theta2j, cluster1j, cluster2j, gammaj, LDj, isLD, step = step.size)
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
      consec_error=0
    }, error = function(e) {
      indicator <<- TRUE
      consec_error <<- consec_error + 1
    })
    if (indicator) {
      next
    }
  }
  close(pb)
  t2=Sys.time()
  time_to_print=round(difftime(t2, t1, units = "secs"),3)
  if(verbose==T){
    cat(paste0("Resampling ends: ",time_to_print," secs\n"))
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
  A$gamma=gamma
  return(A)
}
