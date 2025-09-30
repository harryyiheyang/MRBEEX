#' Build block-diagonal LD matrix by SNP clustering with VIF-based pruning (no shrinkage)
#'
#' @param SNP_DF data.frame with columns: SNP, CHR, BP, P.
#' @param window_size numeric, window (bp) for defining independent SNP centers (default 1.5e6).
#' @param R numeric LD/correlation matrix; row/colnames must be SNP IDs; square.
#' @param kappa_thr numeric, condition-number threshold to trigger VIF pruning (default 30).
#' @param vif_thr numeric, VIF threshold on diag(solve(R_block)) to mark variants high-VIF (default 5).
#'
#' @return list with:
#'   - R1: block-diagonal matrix (no shrinkage) after pruning
#'   - cluster.index: data.frame with SNP, CHR, BP, cluster (after pruning)
#'   - removed_snp: character vector of removed SNP IDs
#' @export
build_blockdiag_ld <- function(SNP_DF,
                               window_size = 1.5e6,
                               R,
                               kappa_thr = 30,
                               vif_thr = 5) {
  if (!is.matrix(R)) stop("R must be a numeric matrix.")
  if (nrow(R) != ncol(R)) stop("R must be square.")
  if (is.null(rownames(R)) || is.null(colnames(R))) stop("R must have row/col names (SNP IDs).")
  if (!all(rownames(R) == colnames(R))) stop("R row/col names must be identical and in the same order.")
  if (isTRUE(anyNA(R))) stop("R contains NA; please handle missing entries first.")
  req_cols <- c("SNP","CHR","BP","P")
  miss <- setdiff(req_cols, names(SNP_DF))
  if (length(miss)) stop("SNP_DF is missing columns: ", paste(miss, collapse=", "))

  SNP_DF <- SNP_DF[order(SNP_DF$CHR, SNP_DF$BP), ]
  keep <- intersect(SNP_DF$SNP, colnames(R))
  if (!length(keep)) stop("No SNP overlap between SNP_DF and R.")
  SNP_DF <- SNP_DF[match(keep, SNP_DF$SNP), , drop = FALSE]

  IndependentSNP <- cluster_snps(df = SNP_DF[, c("SNP","CHR","BP","P")], window_size = window_size)
  IndependentSNP <- IndependentSNP[, c("SNP","CHR","BP","P"), drop = FALSE]

  cluster_df <- clump_cluster(
    df1 = IndependentSNP[, c("SNP","CHR","BP","P")],
    df2 = SNP_DF[,  c("SNP","CHR","BP","P")]
  )
  cluster.index <- cluster_df[, c("SNP","CHR","BP","cluster")]
  rownames(cluster.index) <- cluster.index$SNP

  snps_in_R <- intersect(cluster.index$SNP, colnames(R))
  if (!length(snps_in_R)) stop("After clustering, no SNPs left in R.")
  ord <- match(SNP_DF$SNP, snps_in_R)
  ord <- which(!is.na(ord))
  snp_order <- SNP_DF$SNP[ord]
  snp_order <- snp_order[snp_order %in% snps_in_R]
  R_sub <- R[snp_order, snp_order, drop = FALSE]

  R1 <- R_sub * 0
  to_drop_pos <- integer(0)
  cl_map <- split(snp_order, cluster.index[snp_order, "cluster"])

  for (cid in names(cl_map)) {
    snpj <- cl_map[[cid]]
    indexj <- match(snpj, rownames(R_sub))
    indexj <- indexj[!is.na(indexj)]
    if (length(indexj) < 2) next
    R1[indexj, indexj] <- R_sub[indexj, indexj]
    kap <- tryCatch(kappa(R1[indexj, indexj, drop = FALSE]), error = function(e) Inf)
    if (!is.finite(kap) || kap > kappa_thr) {
      Theta <- tryCatch(solve(R1[indexj, indexj, drop = FALSE]), error = function(e) NULL)
      if (is.null(Theta)) next
      vif <- diag(Theta)
      high_vif_idx <- which(vif > vif_thr)
      if (length(high_vif_idx) > 0) {
        sorted_idx <- high_vif_idx[order(vif[high_vif_idx], decreasing = TRUE)]
        n_drop <- max(length(high_vif_idx) - 1, 1)
        indx <- sorted_idx[seq_len(n_drop)]
        to_drop_pos <- c(to_drop_pos, indexj[indx])
      }
    }
  }
  to_drop_pos <- sort(unique(to_drop_pos))
  kept_snps <- if (length(to_drop_pos)) snp_order[-to_drop_pos] else snp_order
  if (!length(kept_snps)) stop("All SNPs were dropped by VIF pruning; adjust thresholds.")

  R_kept <- R[kept_snps, kept_snps, drop = FALSE]
  cl_kept <- cluster.index[kept_snps, , drop = FALSE]

  R1_final <- R_kept * 0
  cl_map2 <- split(kept_snps, cl_kept[kept_snps, "cluster"])
  for (cid in names(cl_map2)) {
    snpj <- cl_map2[[cid]]
    idx <- match(snpj, rownames(R_kept))
    idx <- idx[!is.na(idx)]
    if (length(idx)) R1_final[idx, idx] <- R_kept[idx, idx]
  }

  list(
    R1 = R1_final,
    cluster.index = cl_kept,
    removed_snp = if (length(to_drop_pos)) snp_order[to_drop_pos] else character(0)
  )
}
