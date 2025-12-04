#' @title Generate data from a Mixture of Normal Distributions (Variance Contribution Model)
#'
#' @description This function generates samples from a mixture of zero-mean normal distributions.
#'              The 'var' parameter specifies the total variance contribution of each
#'              component (i.e., pi_i * sigma_i^2), not the component variance itself.
#'
#' @param m The total number of random samples to generate.
#' @param prop A numeric vector of mixing proportions (weights). These are normalized. Default: c(0.9, 0.09, 0.009, 0.0009).
#' @param var A numeric vector where var[i] represents the variance contribution of component i: V_i = prop[i] * sigma_i^2. Default: c(0.35, 0.4, 0.2, 0.05).
#'
#' @return A numeric vector of length 'm' containing the random variates.
mix_normal_generator <- function(m, prop = c(0.9, 0.09, 0.009, 0.0009), var = c(0.35, 0.4, 0.2, 0.05)) {

  # 1. Input Validation and Normalization
  if (length(prop) != length(var)) {
    stop("Error: The length of 'prop' and 'var' vectors must be equal!")
  }

  sum_prop <- sum(prop)
  prop_norm <- prop / sum_prop # Normalized proportions

  if (any(var < 0)) {
    stop("Error: All elements in 'var' (variance contributions) must be non-negative!")
  }

  n_components <- length(prop_norm)

  if (m < n_components) {
    stop(paste0("Error: Total sample size (m=", m, ") must be at least the number of components (", n_components, ")."))
  }

  # --- CALCULATE COMPONENT STANDARD DEVIATIONS (New Logic) ---
  # Component Variance (sigma_i^2) = (Variance Contribution V_i) / (Mixing Proportion pi_i)
  # Component Standard Deviation (sigma_i) = sqrt(V_i / pi_i)

  # Handle the case where a component has a non-zero contribution (V_i > 0) but a zero proportion (pi_i = 0)
  # Note: If prop_norm[i] is effectively zero, the component is not truly part of the mixture.

  component_var <- var / prop_norm
  # If prop_norm[i] is 0, the division results in Inf or NaN. We must handle this:
  component_var[prop_norm == 0 & var > 0] <- Inf # Impossible case: variance is infinite
  component_var[prop_norm == 0 & var == 0] <- 0   # Safe case: zero contribution, zero variance

  component_sd <- sqrt(component_var)
  # --- END NEW LOGIC ---

  # 2. Determine the number of samples needed for each component (using the robust allocation logic)
  m_counts_initial <- round(m * prop_norm)

  zero_count_indices <- which(m_counts_initial == 0)

  m_counts <- m_counts_initial
  m_counts[zero_count_indices] <- 1

  total_allocated_after_fix <- sum(m_counts)
  deficit_needed <- total_allocated_after_fix - m

  valid_subtraction_indices <- which(m_counts > 1)

  if (deficit_needed > 0) {
    sorted_indices <- valid_subtraction_indices[order(m_counts[valid_subtraction_indices], decreasing = TRUE)]

    for (idx in sorted_indices) {
      can_subtract <- m_counts[idx] - 1
      subtraction_amount <- min(deficit_needed, can_subtract)

      m_counts[idx] <- m_counts[idx] - subtraction_amount
      deficit_needed <- deficit_needed - subtraction_amount

      if (deficit_needed == 0) break
    }
  }

  m_counts[n_components] <- m_counts[n_components] + (m - sum(m_counts)) # Final adjustment

  # 3. Generate data from each component N(mu=0, sd=sigma_i)
  all_data <- numeric(0)

  for (i in 1:n_components) {
    count_i <- m_counts[i]
    sd_i <- component_sd[i]

    if (count_i > 0 && is.finite(sd_i)) {
      # Draw 'count_i' samples from the Normal distribution N(mean=0, sd=sd_i)
      component_data <- stats::rnorm(n = count_i, mean = 0, sd = sd_i)
      all_data <- c(all_data, component_data)
    } else if (count_i > 0 && is.infinite(sd_i)) {
      warning("Warning: Component ", i, " has infinite variance, skipping sample generation for this component.")
    }
  }

  # 4. Shuffle the data
  all_data_shuffled <- sample(all_data)

  return(all_data_shuffled)
}
