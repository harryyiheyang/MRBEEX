#' Build block-diagonal LD matrix by SNP clustering with shrinkage
#'
#' @param SNP_DF data.frame with columns: SNP, CHR, BP, P.
#' @param window_size numeric, window (bp) used to define independent SNP centers (default 1.5e6).
#' @param R numeric LD/correlation matrix; row/colnames must be SNP IDs; square.
#' @param shrink_grid numeric matrix with 2 columns:
#'   - col 1: size cutoffs (e.g., c(50, 10))
#'   - col 2: shrinkage values s in [0,1] (e.g., c(0.05, 0.01))
#'   The matrix is internally sorted by col 1 (descending).
#'   After sorting, col 2 must be nonincreasing (largest cutoff gets >= shrink of the next).
#'   If no rule matches a block, s=0 is used.
#'
#' @return list with:
#'   - R1: sparse block-diagonal matrix (Matrix::dgCMatrix)
#'   - cluster.index: data.frame with SNP, CHR, BP, cluster
#' @export
build_blockdiag_ld <- function(SNP_DF,
                               window_size = 1.5e6,
                               R,
                               shrink_grid = cbind(c(50, 10), c(0.05, 0.01))) {
  ## ---- 0) input checks ----
  req_cols <- c("SNP","CHR","BP","P")
  miss <- setdiff(req_cols, names(SNP_DF))
  if (length(miss)) stop("SNP_DF is missing columns: ", paste(miss, collapse=", "))
  if (!is.matrix(R)) stop("R must be a numeric matrix.")
  if (nrow(R) != ncol(R)) stop("R must be square.")
  if (is.null(rownames(R)) || is.null(colnames(R))) stop("R must have row/col names (SNP IDs).")
  if (!all(rownames(R) == colnames(R))) stop("R row/col names must be identical and in the same order.")
  if (!is.matrix(shrink_grid) || ncol(shrink_grid) != 2) stop("shrink_grid must be a 2-column matrix.")
  if (any(is.na(rownames(R))) || any(is.na(colnames(R)))) stop("R dimnames cannot contain NA.")

  ## ---- 1) keep SNPs present in R ----
  SNP_DF <- SNP_DF[order(SNP_DF$CHR, SNP_DF$BP), ]
  keep <- intersect(SNP_DF$SNP, colnames(R))
  if (!length(keep)) stop("No SNP overlap between SNP_DF and R.")
  SNP_DF <- SNP_DF[match(keep, SNP_DF$SNP), , drop = FALSE]

  ## ---- 2) independent centers by window clustering (your function) ----
  # cluster_snps() is expected to return the center SNPs (smallest P) per window on each chr
  IndependentSNP <- cluster_snps(
    df = SNP_DF[, c("SNP","CHR","BP","P")],
    window_size = window_size
  )
  # Make sure columns exist even if your helper returns subset:
  IndependentSNP <- IndependentSNP[, c("SNP","CHR","BP","P"), drop = FALSE]

  ## ---- 3) assign every SNP to nearest center (your function) ----
  # clump_cluster(df1 = centers, df2 = all SNPs) must add 'cluster' column to df2
  cluster_df <- clump_cluster(
    df1 = IndependentSNP[, c("SNP","CHR","BP","P")],
    df2 = SNP_DF[,  c("SNP","CHR","BP","P")]
  )
  # Normalize shape:
  cluster.index <- cluster_df[, c("SNP","CHR","BP","cluster")]
  rownames(cluster.index) <- cluster.index$SNP

  ## ---- 4) subset and order R to the clustered SNPs ----
  snps_in_R <- intersect(cluster.index$SNP, colnames(R))
  if (!length(snps_in_R)) stop("After clustering, no SNPs left in R.")
  # Keep SNP_DF order (CHR,BP), but only those present in R:
  ord <- match(SNP_DF$SNP, snps_in_R)
  ord <- which(!is.na(ord))
  snp_order <- SNP_DF$SNP[ord]
  snp_order <- snp_order[snp_order %in% snps_in_R]
  R_sub <- R[snp_order, snp_order, drop = FALSE]

  ## ---- 5) check and sort shrink_grid ----
  sg <- shrink_grid
  storage.mode(sg) <- "double"
  # sort by cutoff desc
  o <- order(sg[,1], decreasing = TRUE)
  sg <- sg[o, , drop = FALSE]
  # enforce nonincreasing shrink along this order
  if (any(diff(sg[,2]) > 1e-12)) {
    stop("After sorting shrink_grid by column 1 (descending), column 2 must be nonincreasing.")
  }
  # helper to get shrink s given block size k
  pick_shrink <- function(k) {
    # first rule whose cutoff is strictly less than k? or <= ?
    # your example used (>50). We'll use '>' to match that.
    idx <- which(k > sg[,1])[1]
    if (length(idx) && !is.na(idx)) return(sg[idx,2])
    0
  }

  ## ---- 6) build block-diagonal R1 with shrink by block size ----
  # Start with all zeros; fill block by block
  R1 <- matrix(0, nrow = nrow(R_sub), ncol = ncol(R_sub),
               dimnames = dimnames(R_sub))

  # group SNPs by 'cluster' index (based on snp_order)
  cl_map <- split(snp_order, cluster.index[snp_order, "cluster"])

  for (cid in names(cl_map)) {
    snpj <- cl_map[[cid]]
    idx  <- match(snpj, rownames(R_sub))
    idx  <- idx[!is.na(idx)]
    if (!length(idx)) next
    # copy block
    block <- R_sub[idx, idx, drop = FALSE]
    # apply shrink rule if size exceeds thresholds
    k <- length(idx)
    s <- pick_shrink(k)  # s in [0,1]
    if (s > 0) {
      block <- (1 - s) * block + s * diag(k)
    }
    R1[idx, idx] <- block
  }

  list(R1 = R1, cluster.index = cluster.index[rownames(R1), , drop = FALSE])
}
