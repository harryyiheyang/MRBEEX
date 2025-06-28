#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List construct_sparse_blockwise_LD_cpp(const arma::sp_mat& LD,
                                       const IntegerVector& cluster_index,
                                       const IntegerVector& cluster_sampling,
                                       double admm_rho) {
  if (LD.n_rows != LD.n_cols) {
    stop("LD matrix must be square");
  }
  if (LD.n_rows != cluster_index.size()) {
    stop("LD matrix dimensions must match cluster_index length");
  }
  if (admm_rho < 0) {
    warning("admm_rho should be positive");
  }

  int n_vars = cluster_index.size();
  int n_clusters = cluster_sampling.size();

  std::unordered_map<int, std::vector<int>> cluster_map;
  for (int i = 0; i < n_vars; ++i) {
    int cluster_id = cluster_index[i];
    cluster_map[cluster_id].push_back(i);
  }

  std::vector<int> block_sizes(n_clusters);
  std::vector<std::vector<int>> block_indices(n_clusters);
  int total_size = 0;

  for (int i = 0; i < n_clusters; ++i) {
    int cluster_id = cluster_sampling[i];
    auto it = cluster_map.find(cluster_id);
    if (it == cluster_map.end()) {
      stop("Cluster ID " + std::to_string(cluster_id) + " not found in cluster_index");
    }
    block_indices[i] = it->second;
    block_sizes[i] = it->second.size();
    total_size += block_sizes[i];
  }

  if (total_size == 0) {
    stop("Total size is zero");
  }

  arma::sp_mat Thetaj(total_size, total_size);
  arma::sp_mat Thetarhoj(total_size, total_size);
  arma::sp_mat TCj(total_size, total_size);
  arma::sp_mat LDj(total_size, total_size);

  std::vector<int> indj_vec;
  indj_vec.reserve(total_size);

  int current_pos = 0;

  arma::mat LD_dense_all(LD);  // convert full matrix to dense once

  for (int i = 0; i < n_clusters; ++i) {
    const std::vector<int>& vars = block_indices[i];
    int block_size = vars.size();

    if (block_size == 0) {
      continue;
    }

    try {
      arma::uvec var_indices = arma::conv_to<arma::uvec>::from(vars);
      arma::mat LD_dense = LD_dense_all.submat(var_indices, var_indices);

      arma::mat L_test;
      if (!arma::chol(L_test, LD_dense, "lower")) {
        stop("Block " + std::to_string(i + 1) + " is not positive definite");
      }

      arma::mat I = arma::eye(block_size, block_size);
      arma::mat Theta_block = arma::solve(LD_dense, I, arma::solve_opts::likely_sympd);
      arma::mat Thetarho_block = arma::solve(LD_dense + admm_rho * I, I, arma::solve_opts::likely_sympd);
      arma::mat TC_block;
      if (!arma::chol(TC_block, LD_dense, "upper")) {
        stop("Cholesky decomposition failed for block " + std::to_string(i + 1));
      }

      int end_pos = current_pos + block_size - 1;

      LDj.submat(current_pos, current_pos, end_pos, end_pos) = sp_mat(LD_dense);
      Thetaj.submat(current_pos, current_pos, end_pos, end_pos) = sp_mat(Theta_block);
      Thetarhoj.submat(current_pos, current_pos, end_pos, end_pos) = sp_mat(Thetarho_block);
      TCj.submat(current_pos, current_pos, end_pos, end_pos) = sp_mat(TC_block);

      for (int var_idx : vars) {
        indj_vec.push_back(var_idx + 1);
      }

      current_pos += block_size;

    } catch (const std::exception& e) {
      stop("Error processing block " + std::to_string(i + 1) + ": " + std::string(e.what()));
    }
  }

  IntegerVector indj = wrap(indj_vec);

  return List::create(
    Named("indj") = indj,
    Named("LDj") = LDj,
    Named("Thetaj") = Thetaj,
    Named("Thetarhoj") = Thetarhoj,
    Named("TCj") = TCj
  );
}

// [[Rcpp::export]]
List get_matrix_info_cpp(const arma::sp_mat& mat) {
  return List::create(
    Named("n_rows") = mat.n_rows,
    Named("n_cols") = mat.n_cols,
    Named("n_nonzero") = mat.n_nonzero,
    Named("density") = (double)mat.n_nonzero / (mat.n_rows * mat.n_cols)
  );
}
