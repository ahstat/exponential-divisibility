## Preparation for the inversion ----

#' Create a log-scale of length.out points to be evaluated
#' @param length.out Number of evaluated values t (a valuated log-scale is used)
#' @param logTmin The smallest positive value evaluated is given by exp(logTmin)
#' @param logTmax The extremes evaluated are given by -exp(logTmax) and +exp(logTmax), at most 709
#' @export
precalculated_points = function(length.out = 1000, logTmin = -100, logTmax = 100) {
  ts = exp(seq(from = logTmin, to = logTmax, length.out = length.out))
  ts = c(-rev(ts), 0, ts)
  return(ts)
}

#' Create a pre-calculated characteristic function corresponding to phi_func.
#' Given a point t, instead of computing the result from scratch, it will output
#' the value of the closest t.
#' @param alpha Vector of positive numbers of length at least 1
#' @param ts Vector of points t on which the function is evaluated
#' @param K Number of terms of the sum on the initial radius of convergence
#' @return Function of the pre-calculated characteristic function the ts points
#' @export
precalculated_phi_func = function(alpha = c(1, 2), ts = precalculated_points(), K = 10000) {
  phi_output = exp(log_phi_func(ts, alpha, K, verbose = FALSE))

  cf_out = function(t) {
    idx = which.min(abs(t - ts))
    phi_output[idx]
  }
  cf_out2 = function(t) {
    sapply(t, cf_out)
  }
  return(cf_out2)
}

#' Create a not pre-calculated characteristic function corresponding to phi_func,
#' which will give much slower results.
#' @param alpha Vector of positive numbers of length at least 1
#' @param K Number of terms of the sum on the initial radius of convergence
#' @return Function of the not pre-calculated characteristic function
#' @export
not_precalculated_phi_func = function(alpha = c(1, 2), K = 10000) {
  cf_out2 = function(t) {
    if(length(t) == 1) {
      if(t == 1e+300) { # point evaluated in CharFunToolR::cf2DistGP for the inversion
        return(0)
      }
    }
    return(phi_func(t, alpha, K, verbose = FALSE))
  }
  return(cf_out2)
}

## Characteristic function inversion ----

# Different methods exist to retrieve the cumulative distribution
# from the characteristic function: the Lévy inversion formula for positive
# variables, or the Gil-Pelaez (1951) inversion. We use the Gil-Palaez inversion.
#
# We use the CharFunToolR package developed by the Viktor Witkovský's team
#
# The reference is (2016 Witkovský)
# "Numerical inversion of a characteristic function:
#  An alternative tool to form the probability distribution
#  of output quantity in linear measurement models"
#                      https://arxiv.org/pdf/1701.08299.pdf
#
# The first version of the package is https://github.com/Simkova/CharFun
# with last release on 24/05/2017
# install_github("Simkova/CharFun")
# library(CharFun)
#
# The second version of the package is https://github.com/gajdosandrej/CharFunToolR
# with last release on 28/12/2018
# install_github("gajdosandrej/CharFunToolR") # need V8: sudo apt-get install libv8-dev
#
# We use this second version library(CharFunToolR)

#' Rule of thumb for the parameter T in the fourier_inversion function
#' @param cf The characteristic function to inverse
#' @param rel.tol Accepted tolerance (to select as small as possible)
#' @param range_of_search_the_tol Range of search (default should be ok)
#' @return Rule of thumb for selecting the parameter T
#' @export
rule_of_thumb_T = function(cf, rel.tol = 1e-4, range_of_search_the_tol = c(1e-3, 1e10)) {
  # Next we choose T such that, for t > T: abs(cf(t)/t) <= rel.tol
  # We thus solve: abs(cf(t)/t) - rel.tol
  h = function(t) {
    abs(cf(t)/t) - rel.tol
  }

  # t = 10^(log(range_of_search_the_tol, 10)[1]:log(range_of_search_the_tol, 10)[2])
  # plot(t, abs(cf(t)/t), type = "l", log = "y")
  if(h(range_of_search_the_tol[1]) < 0) {
    # T = 1e-3 is working, but it's dubious.
    stop("rel.tol should be small, as close as possible to 0")
  }
  if(h(range_of_search_the_tol[2]) > 0) {
    # The correct T is larger that 1e10...
    stop("rel.tol should be larger")
  }

  T_selected = round(uniroot(h, range_of_search_the_tol, tol = .Machine$double.eps)$root)
  # h(T_selected) # ~ 0
  return(T_selected)
}

#' Rule of thumb for the parameter dt in the fourier_inversion function
#' @param xMean Mean of the random variable of interest
#' @param xStd Standard deviation of the random variable of interest
#' @param k Number of sigma of distance, 6 by default
#' @param T_selected Selected T (that can also be computed with the rule of thumb)
#' @return Rule of thumb for selecting the parameter dt
#' @export
rule_of_thumb_N = function(xMean, xStd, k = 6, T_selected) {
  # Rule of thumb of the documentation is detailed in the arxiv page 7
  # k-sigma with k=6 for example.
  # dt = (2*pi)/(xMax - xMin)
  # with: xMax = xMean + k*xStd, and xMin = xMean - 6*xStd
  # Here xMean and xStd are the mean and sd of the random variable of interest
  # In our case we know the cumulants so:
  xMax = xMean + k*xStd
  xMin = xMean - k*xStd
  dt = (2*pi)/(xMax - xMin)
  N_selected = round(T_selected/dt)
  return(dt)
}

#' Inversion of the probabilistic characteristic function to retrieve the proba
#' @param x Interval where the cumulative distribution function needs to be calculated
#' @param cf Characteristic function to inverse
#' @param params_T Selected T (function rule_of_thumb_T could be used)
#' @param params_N Selected N (function rule_of_thumb_N could be used)
#' @param verbose Additional information
#' @return Perform the inversion of the characteristic function following CharFunToolR::cf2DistGP
#' @export
fourier_inversion = function(x = seq(from = -3, to = 3, length.out = 101),
                             cf,
                             T_selected = 100,
                             N_selected = 100,
                             verbose = TRUE) {
  result = CharFunToolR::cf2DistGP(cf = cf,
                                   x = x,
                                   options = list(isPlot = verbose,
                                                  N = N_selected,
                                                  T = T_selected))
  return(result)
}
