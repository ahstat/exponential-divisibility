## Characteristic function for cutting an exponential distribution
# such that sum(alpha_i X_i) = Z, with Z following Exp(1).
.info_from_alpha = function(alpha = c(1, 2, 3)) {
  if(length(alpha) == 0) {
    stop("`alpha` should have positive length")
  }
  if(any(alpha <= 0)) {
    stop("All elements of alpha should be positive")
  }

  alpha = sort(alpha)
  N = length(alpha)
  is_trivial = length(unique(alpha)) == 1

  # Get `p`
  rle_alpha = rle(sort(alpha, decreasing = TRUE))
  p = rle_alpha[[1]][1] # multiplicity of alpha[N]

  return(list(alpha = alpha, N = N, p = p, is_trivial = is_trivial))
}

.get_D_from_t = function(t, alpha_N_minus_p, alpha_N) {
  D = rep(NA, length(t))
  idx_in = which(abs(t) <= alpha_N)
  D[idx_in] = 0
  idx_out = which(abs(t) > alpha_N)
  D[idx_out] = ceiling((log(abs(t[idx_out])) - log(alpha_N)) / (log(alpha_N) - log(alpha_N_minus_p)))
  return(D)
}

.log_phi_main_func = function(t, alpha = c(1, 2), K = 100) {
  if(max(abs(t)) > max(alpha)) {
    stop("Out of the radius of convergence")
  }

  k = 1:K

  # vector of \sum_{l=1}^{N} \alpha_l^k over k, of length K
  bottom_term = sapply(k, function(k_idx){sum(alpha^k_idx)})

  # vector of \frac{i^k}{k \sum_{l=1}^{N} \alpha_l^k}
  terms_without_t = (1i)^k / (k * bottom_term)

  output = sapply(t, function(t_val) {
    # na.rm to prevent division by infinite
    sum(t_val^k * terms_without_t, na.rm = TRUE)
  })
  return(output)
}

.log_phi_compound_func = function(t, alpha_j_k_list, alpha_N, alpha, D, p, K) {
  my_terms = list()

  # First term
  first_term = -(1/p) * log(1 - (1/alpha_N)*1i*t)
  my_terms[[1]] = first_term

  # Middle terms
  if(D > 1L) {
    for(d in 2:D) {
      coeffs = apply(alpha_j_k_list[[d-1]], 1, prod) / alpha_N^d
      middle_terms = sapply(coeffs, function(coeff) {log(1 - coeff*1i*t)})
      if(length(t) > 1) {
        my_terms[[d]] = (-1)^d/p^d * apply(middle_terms, 1, sum)
      } else {
        my_terms[[d]] = (-1)^d/p^d * middle_terms
      }
    }
  }

  # Right term
  coeffs = apply(alpha_j_k_list[[D]], 1, prod) / alpha_N^D
  right_terms = sapply(coeffs, function(coeff) {.log_phi_main_func(coeff*t, alpha, K)})
  if(length(t) > 1) {
    right_term = (-1)^D/p^D * apply(right_terms, 1, sum)
  } else {
    right_term = (-1)^D/p^D * right_terms
  }
  my_terms[[D+1]] = right_term

  return(Reduce("+", my_terms))
}

#' Log of the characteristic function, see the main function \code{phi_func}.
#' @param t Real value where the characteristic function is evaluated
#' @param alpha Vector of positive numbers of length at least 1
#' @param K Number of terms of the sum on the initial radius of convergence
#' @param verbose Print additional information if TRUE
#' @return Evaluation of the log of the characteristic function at point t
#' @export
log_phi_func = function(t, alpha = c(1, 2, 3), K = 100, verbose = TRUE) {
  # `alpha` is checked and sorted
  # N is the length of `alpha`
  # p is the multiplicity of alpha[N]
  # is_trivial corresponds to the case where all elements of alpha are equal
  list_info = .info_from_alpha(alpha)
  alpha = list_info$alpha
  N = list_info$N
  p = list_info$p
  is_trivial = list_info$is_trivial

  ## If all elements of `alpha` are equal (in particular also true when N = 1)
  if(is_trivial) {
    k = 1/N # shape of the gamma
    theta = 1/alpha[1] # scale of the gamma
    # logarithm of the characteristic function of a gamma(k, theta)
    return(-k * log(1 - 1i * t * theta))
  }

  ## Other case
  D_vals = .get_D_from_t(t, alpha[N-p], alpha[N])
  # prevent slow convergence when close to the radius of convergence:
  # (see tests of `.log_phi_main_func` for example)
  D_vals = D_vals + 1
  if(verbose) {
    print(paste0(t[1], " (length of vector is ", length(t), ")"))
    print(paste0("Max D is ", max(D_vals)))
  }

  log_phi_output = rep(NA, length(t))
  for(D in sort(unique(D_vals))) {
    if(verbose) {
      print(D)
    }
    idx = which(D_vals == D) # t for which it is the current D
    if(D == 0L) {
      log_phi_output[idx] = .log_phi_main_func(t[idx], alpha, K)
    } else {
      alpha_j_k_list = list()
      for(d in 1:D) {
        alpha_j_k_list[[d]] = expand.grid(replicate(d, alpha[1:(N-p)], simplify = FALSE))
      }
      log_phi_output[idx] = .log_phi_compound_func(t[idx], alpha_j_k_list, alpha[N], alpha, D, p, K)
    }
  }

  return(log_phi_output)
}

#' Characteristic function of the variable X_1 such that, given i.i.d. copies
#' X_1, ... X_n, the linear combination alpha_1 X_1 + ... + alpha_n X_n follows
#' an exponential distribution.
#' @param t Real value where the characteristic function is evaluated
#' @param alpha Vector of positive numbers of length at least 1
#' @param K Number of terms of the sum on the initial radius of convergence
#' @param verbose Print additional information if TRUE
#' @return Evaluation of the characteristic function at point t
#' @export
phi_func = function(t, alpha = c(1, 2, 3), K = 100, verbose = TRUE) {
  return(exp(log_phi_func(t, alpha, K, verbose)))
}
