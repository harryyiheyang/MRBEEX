#' @title Sparse Prediction of xQTL Effects using SuSiE and CARMA
#' @description
#' This function performs infomative xQTL selection and sparse prediction using SuSiE and CARMA.
#'
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix.
#' @param xQTL.method The method used in purifying the xQTLs. SuSiE or CARMA can be used here, where the latter can be more accurate but much most computationally costly. Defaults is SuSiE.
#' @param xQTL.max.L When choosing \code{"SuSiE"}, the maximum number of L in estimating the xQTL effects. Defaults to 10.
#' @param xQTL.cred.thres When choosing \code{"SuSiE"}, the minimum empirical posterior inclusion probability (PIP) used in getting credible sets of xQTL selection. Defaults to \code{0.95}.
#' @param unmappable_effects The method for modeling unmappable effects: "none", "inf".
#' @param estimate_residual_method The method used for estimating residual variance. For the original SuSiE model, "MLE" and "MoM" estimation is equivalent, but for the infinitesimal model, "MoM" is more stable. We recommend using "Servin_Stephens" when n < 80 for improved coverage.
#' @param xQTL.Nvec When choosing \code{"SuSiE"}, the vector of sample sizes of exposures.
#' @param xQTL.weight When choosing \code{"SuSiE"}, the vector of weights used in specifying the prior weights of SuSiE. Defaults to \code{NULL}.
#' @param outlier.switch When choosing \code{"CARMA"}, an indicator of whether turning on outlier detection. Defaults to \code{F}.
#' @param Annotation When choosing \code{"CARMA"}, the annotation matrix of SNP. Default is NULL.
#' @param output.labels When choosing \code{"CARMA"}, output directory where output will be written while CARMA is running. Defaults to \code{NULL}, meaning that a temporary folder will be created and automatically deleted upon completion of the computation.
#' @param carma.iter When choosing \code{"CARMA"}, the maximum iterations for EM algorithm to run. Defaults to 5.
#' @param carma.inner.iter When choosing \code{"CARMA"}, the maximum iterations for Shotgun algorithm to run per iteration within EM algorithm. Defaults to 5.
#' @param xQTL.max.num When choosing \code{"CARMA"}, the maximum number of causal variants assumed per locus, which is similar to the number of single effects in SuSiE. Defaults to 10.
#' @param carma.epsilon.threshold When choosing \code{"CARMA"}, the convergence threshold measured by average of Bayes factors. Defaults to \code{1e-3}.
#'
#' @importFrom MASS rlm ginv
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom susieR coef.susie susie susie_rss susie_get_cs
#'
#' @return A list containing the results of the MRBEEX analysis using different methods:
#' \describe{
#'   \item{\code{bXest}}{A matrix where each column represents \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXestse}}{A matrix where each column contains the standard errors of \code{R * beta_PLS} for each exposure.}
#'   \item{\code{bXest0}}{A matrix where each column represents \code{beta_PLS} for each exposure.}
#'   \item{\code{bXestse0}}{A matrix where each column contains the standard errors of \code{beta_PLS} for each exposure.}
#'   \item{\code{xQTLfitList}}{A list containing the xQTL selection results, which can be used in \code{CisMRBEEX} by setting \code{xQTLfitList = xQTLfitList}.}
#' }
#' @export

Sparse_Prediction=function(bX,bXse,LD,xQTL.Nvec,xQTL.method="SuSiE",
                           xQTL.max.L=10,xQTL.weight=NULL,xQTL.cred.thres=0.95,
                           estimate_residual_method="MoM",unmappable_effects="inf",
                           outlier.switch=T,Annotation=NULL,output.labels=NULL,
                           carma.iter=5,carma.inner.iter=5,xQTL.max.num=10,
                           carma.epsilon.threshold=1e-3){
p=ncol(bX)
m=nrow(bX)

if(xQTL.method=="SuSiE"){
xQTLfitList=list()
if(is.null(xQTL.weight)==T){
xQTL.weight=rep(1,m)
}
for(i in 1:p){
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=xQTL.max.L,max_iter=1000,prior_weights=xQTL.weight,estimate_residual_method=estimate_residual_method,unmappable_effects=unmappable_effects,convergence_method=ifelse(unmappable_effects=="none","elbo","pip"))
fit=susie_rss(z=bX[,i]/bXse[,i],R=LD,n=xQTL.Nvec[i],L=length(susie_get_cs(fit,coverage=xQTL.cred.thres)$cs)+1,max_iter=1000,prior_weights=xQTL.weight,estimate_residual_method=estimate_residual_method,unmappable_effects=unmappable_effects,convergence_method=ifelse(unmappable_effects=="none","elbo","pip"))
xQTLfitList[[i]]=fit
}
}

if(xQTL.method=="CARMA"){
if (!requireNamespace("CARMA", quietly = TRUE)) {
stop("Package 'CARMA' is required for xQTL.method='CARMA', but is not installed. ",
"Please install it with install.packages('CARMA').")
}
z.list=ld.list=w.list=lambda.list=list()
for(i in 1:p){
z.list[[i]]=bX[,i]/bXse[,i]
ld.list[[i]]=LD
w.list[[i]]=Annotation
lambda.list[[i]]=1
}
is.detect=F
if(is.null(output.labels)==T){
output.labels=tempfile("my_tmpdir_")
dir.create(output.labels)
is.delect=T
}
if(is.null(Annotation)==F){
fitxQTL=CARMA::CARMA(z.list,ld.list,w.list,lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}else{
fitxQTL=CARMA::CARMA(z.list,ld.list,lambda.list = lambda.list,outlier.switch=outlier.switch,num.causal=xQTL.max.num,printing.log=F,all.iter=carma.iter,all.inner.iter=carma.inner.iter,epsilon.threshold=carma.epsilon.threshold,output.labels=output.labels)
}
if(is.delect==T){
unlink(output.labels, recursive = TRUE, force = TRUE)
}
xQTLfitList=fitxQTL
}

return(xQTLfitList)
}
