#' Cluster SNPs (greedy by P) and assign every SNP to exactly one cluster
#'
#' @description
#' One-pass pipeline for C+T style clustering and membership assignment.
#' For each chromosome, repeatedly take the SNP with the smallest P as a center,
#' group all SNPs within \code{±window_size} bp into that cluster, remove them,
#' and continue until all SNPs are clustered. Then assign every SNP in \code{df}
#' to exactly one cluster by non-equi interval join with explicit tie-breaking
#' in overlapping regions: smallest distance to center, then smaller center P,
#' then smaller \code{cluster_id}.
#'
#' @param df A data.frame/data.table with columns: \code{SNP}, \code{CHR}, \code{BP}, and \code{P}.
#' @param window_size Integer, window radius in base pairs (e.g., \code{1e6} for ±1 Mb).
#'
#' @return A list with two data.tables:
#' \itemize{
#'   \item \code{centers}: cluster centers (like df1), columns \code{SNP, CHR, BP, P, ClusterSize, left, right, cluster_id}.
#'   \item \code{assigned}: original \code{df} annotated with \code{cluster} (ID), distance to center \code{dist_bp},
#'         and center metadata \code{center_SNP, center_BP, center_P}.
#' }
#'
#' @import data.table
#' @export
cluster_snps_and_assign <- function(df, window_size = 1e6) {
# ----- input checks -----
df <- data.table::as.data.table(data.table::copy(df))
stopifnot(all(c("SNP","CHR","BP","P") %in% names(df)))
if (!is.numeric(window_size) || length(window_size) != 1L || window_size < 0) {
  stop("window_size must be a non-negative length-1 numeric.")
}

# normalize types / order
df[, `:=`(CHR = as.integer(CHR), BP = as.integer(BP))]
data.table::setorder(df, CHR, BP)

# ----- Step 1: greedy C+T clustering per chromosome (build centers) -----
centers_list <- vector("list", 0L); j <- 1L

for (chr in unique(df$CHR)) {
  d <- data.table::copy(df[CHR == chr])
  if (nrow(d) == 0L) next

  # loop: pick min-P center, take ±window, remove, repeat
  while (nrow(d) > 0L) {
    ctr <- d[which.min(P)]
    left  <- ctr$BP - as.integer(window_size)
    right <- ctr$BP + as.integer(window_size)
    in_win <- d[BP >= left & BP <= right]

    centers_list[[j]] <- data.table::data.table(
      SNP = ctr$SNP, CHR = ctr$CHR, BP = ctr$BP, P = ctr$P,
      ClusterSize = nrow(in_win),
      left = left, right = right
    )
    j <- j + 1L

    # remove clustered SNPs in this window
    d <- d[!(BP >= left & BP <= right)]
  }
}

if (length(centers_list) == 0L) {
  # no centers (empty df)
  return(list(
    centers  = data.table::data.table(),
    assigned = data.table::data.table()
  ))
}

centers <- data.table::rbindlist(centers_list, use.names = TRUE)
data.table::setorder(centers, CHR, BP)
centers[, cluster_id := .I]  # stable cluster index across all chromosomes

# ----- Step 2: non-equi interval join → candidates (may have multiple hits per SNP) -----
cand <- centers[df,
                on = .(CHR, left <= BP, right >= BP),
                .(SNP = i.SNP, CHR = i.CHR, BP = i.BP, P = i.P,
                  cluster = cluster_id,
                  center_SNP = x.SNP, center_BP = x.BP, center_P = x.P),
                nomatch = 0L]

if (nrow(cand) < nrow(df)) {
  message("Note: some SNPs were not covered by any window; they will be absent from `assigned`.")
}
# ----- Step 3: remove overlapped SNP -----
cand[, dist_bp := abs(BP - center_BP)]
data.table::setorder(cand, CHR, SNP, BP, dist_bp, center_P, cluster)
assigned <- cand[, .SD[1L], by = .(CHR, SNP, BP)]
data.table::setorder(assigned, CHR, BP)
list(centers = centers[], assigned = assigned[])
}
