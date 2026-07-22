#' Cluster pruning by Local Pratt Index (retain up to 5 per eligible cluster)
#'
#' Given \code{cluster_info} produced by \code{clump_cluster} and two row-aligned
#' effect-size vectors \code{marginal_effect} and \code{direct_effect}, this function
#' prunes rows within clusters whose size exceeds \code{cutoff} by selecting, per such
#' cluster, at most the top 5 variants ranked by the local Pratt index accumulating their normalized contributions
#' until reaching \code{contribution}; if the threshold is not reached, the top 5 are kept.
#' Clusters with size \emph{less than or equal to} \code{cutoff} are fully retained.
#'
#' @param cluster_info A data frame from \code{clump_cluster}, containing at least a
#'   numeric column \code{cluster} giving the cluster index for each row; must be
#'   row-aligned with the effect vectors.
#' @param marginal_effect Numeric vector of GWAS marginal effect sizes, row-aligned
#'   with \code{cluster_info}. Must be standardized to a comparable scale
#'   (e.g., z-scores) before calling this function.
#' @param direct_effect Numeric vector of direct effect sizes (from fine-mapping or PRS),
#'   row-aligned with \code{cluster_info}. Must be standardized to a comparable
#'   scale before calling. If estimated by SBayesRC, a common standardization is
#'   \eqn{\beta*\sqrt{2*f*(1-f)}}, where \eqn{f} is the allele frequency.
#' @param cutoff Integer. Maximum cluster size threshold for pruning: clusters with
#'   size \code{<= cutoff} are entirely kept; clusters with size \code{> cutoff} are
#'   pruned using the Local Pratt rule.
#' @param contribution Numeric in \eqn{(0,1]}. Target cumulative contribution (based on
#'   normalized Local Pratt weights) used to decide how many top variants to retain
#'   per pruned cluster; no more than 5 will be kept.
#' @param max_size The maximum size of each cluster when using the Pratt index to cut off.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{cluster_info_updated}: the pruned subset of \code{cluster_info}
#'   \item \code{selected_global_idx}: global row indices retained
#'   \item \code{dropped_global_idx}: global row indices removed
#' }
#'
#' @export

block_cutoff <- function(cluster_info, marginal_effect, direct_effect, cutoff = 5, max_size=5, contribution = 0.9){
  stopifnot(nrow(cluster_info) == length(marginal_effect), length(marginal_effect) == length(direct_effect))
  cnt <- table(cluster_info$cluster)
  eligible <- as.integer(names(cnt)[cnt > cutoff])
  smalls <- as.integer(names(cnt)[cnt <= cutoff])
  n <- nrow(cluster_info)
  keep_idx <- cluster_info$cluster %in% smalls
  for (cl in eligible) {
    ids <- which(cluster_info$cluster == cl)
    pr <- marginal_effect[ids] * direct_effect[ids]
    w <- abs(pr)
    if (length(ids) == 0) next
    if (sum(w) > 0) {
      w <- w / sum(w)
      ord <- order(w, decreasing = TRUE)
      cumw <- cumsum(w[ord])
      k_hit <- which(cumw >= contribution)[1]
      k <- min(max_size, ifelse(is.na(k_hit), length(ids), k_hit))
      sel <- ord[seq_len(k)]
    } else {
      k <- min(max_size, length(ids))
      sel <- seq_len(k)
    }
    keep_idx[ids[sel]] <- TRUE
  }
  list(
    cluster_info_updated = cluster_info[keep_idx, , drop = FALSE],
    selected_global_idx = which(keep_idx),
    dropped_global_idx = which(!keep_idx)
  )
}
