#' Construct Sparse Blockwise LD Matrix
#'
#' This function constructs a block-diagonal LD matrix and its related inverses
#' using block-wise bootstrap based on cluster assignments.
#'
#' @param LD A square sparse correlation/covariance matrix of class \code{dgCMatrix}.
#' @param cluster.index An integer vector indicating the cluster assignment of each variable (1-based).
#' @param cluster.sampling An integer vector of sampled cluster IDs (with replacement).
#' @param admm.rho A positive numeric value for ADMM regularization (used in ridge inverse).
#'
#' @return A list with components:
#' \describe{
#'   \item{indj}{An integer vector of the re-ordered variable indices (1-based).}
#'   \item{LDj}{Block-diagonal sparse matrix (class \code{dgCMatrix}).}
#'   \item{Thetaj}{Inverse of LDj for each block (block-diagonal).}
#'   \item{Thetarhoj}{Inverse of LDj + rho * I for each block (block-diagonal).}
#'   \item{TCj}{Cholesky factor of LDj for each block (block-diagonal).}
#' }
#' @export
construct_sparse_blockwise_LD <- function(LD, cluster.index, cluster.sampling, admm.rho=0) {
  .Call(`_MRBEEX_construct_sparse_blockwise_LD_cpp`, LD, cluster.index, cluster.sampling, admm.rho)
}
