library(zeallot)
test_that(".info_from_alpha answers correct information for alpha, N, p, and is_trivial", {
  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(1, 2, 3))
  expect_equal(alpha, c(1, 2, 3))
  expect_equal(N, 3)
  expect_equal(p, 1)
  expect_equal(is_trivial, FALSE)

  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(3, 1, 2))
  expect_equal(alpha, c(1, 2, 3))
  expect_equal(N, 3)
  expect_equal(p, 1)
  expect_equal(is_trivial, FALSE)

  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(3, 3))
  expect_equal(alpha, c(3, 3))
  expect_equal(N, 2)
  expect_equal(p, 2)
  expect_equal(is_trivial, TRUE)

  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(8))
  expect_equal(alpha, c(8))
  expect_equal(N, 1)
  expect_equal(p, 1)
  expect_equal(is_trivial, TRUE)

  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(9, 9, 7, 7, 3, 1, 3, 3))
  expect_equal(alpha, c(1, 3, 3, 3, 7, 7, 9, 9))
  expect_equal(N, 8)
  expect_equal(p, 2)
  expect_equal(is_trivial, FALSE)

  c(alpha, N, p, is_trivial) %<-% .info_from_alpha(c(1, 4, 4, 4))
  expect_equal(alpha, c(1, 4, 4, 4))
  expect_equal(N, 4)
  expect_equal(p, 3)
  expect_equal(is_trivial, FALSE)

  expect_error(.info_from_alpha(c()))
  expect_error(.info_from_alpha(c(-1)))
  expect_error(.info_from_alpha(c(-1, 1)))
})

test_that(".get_D_from_t verifies the inequality", {
  ## Plotting of D as a function of t
  # library(zeallot)
  t = seq(from = -100, to = 100, length.out = 1000)
  alpha = c(1, 2) # c(0.2, 0.5, 2)
  c(alpha, N, p, is_trivial) %<-% expdiv:::.info_from_alpha(alpha)
  D = expdiv:::.get_D_from_t(t, alpha[N-p], alpha[N])
  # plot(t, D, type = "l")

  ## Common test that the condition on D is verified
  t = seq(from = -100, to = 100, length.out = 1000)
  alpha = c(1, 2)
  c(alpha, N, p, is_trivial) %<-% expdiv:::.info_from_alpha(alpha)
  D = expdiv:::.get_D_from_t(t, alpha[N-p], alpha[N])
  expect_equal(length(D), length(t))
  expect_true(all((alpha[N-p]/alpha[N])^D * abs(t) <= alpha[N]))
  # plot((alpha[N-p]/alpha[N])^D * abs(t) - alpha[N]) # never reach 0 and quite optimal selection of D

  ## Check when no value outside [-alpha[N], alpha[N]]
  t = seq(from = -1, to = 1, length.out = 1000)
  alpha = c(1, 2)
  c(alpha, N, p, is_trivial) %<-% expdiv:::.info_from_alpha(alpha)
  D = expdiv:::.get_D_from_t(t, alpha[N-p], alpha[N])
  expect_equal(length(D), length(t))
  expect_true(all((alpha[N-p]/alpha[N])^D * abs(t) <= alpha[N]))
  expect_true(all(D == 0))

  ## Check when no value inside [-alpha[N], alpha[N]]
  t = seq(from = 5, to = 10, length.out = 1000)
  alpha = c(1, 2)
  c(alpha, N, p, is_trivial) %<-% expdiv:::.info_from_alpha(alpha)
  D = expdiv:::.get_D_from_t(t, alpha[N-p], alpha[N])
  expect_equal(length(D), length(t))
  expect_true(all((alpha[N-p]/alpha[N])^D * abs(t) <= alpha[N]))
  expect_true(all(D > 0))

  ## .get_D_from_t is not tested for N = 1 because it is not needed in this case
})

test_that(".log_phi_main_func verifies the equation in [-1, 1]", {
  # 1 is the radius of convergence for the characteristic function of
  # and exponential
  t = seq(from = -1, 1, length.out = 1000)
  alpha = c(1, 2)
  K = 1e3
  # Compute the product of the characteristic, which should be an exponential
  pred = apply(sapply(alpha, function(alpha_j) {expdiv:::.log_phi_main_func(alpha_j*t, alpha, K)}), 1, sum)
  true = -log(1-1i*t)
  # plot(t, abs(pred - true))
  expect_equal(pred, true, tolerance = 1e-5)
  # The convergence is not good close to -1 and 1, and this is not improving
  # by increasing K.
  # In the `log_phi_func` function, it is solved by taking D one more than the
  # ideal one, so that the sum close to 1 are converted into sum close to 1/2 for
  # example.

  ## Test with other alpha
  K = 1e3
  t = seq(from = -1, 1, length.out = 1000)
  true = -log(1-1i*t)
  alphas = list(c(1, 3), c(1, 1, 2), c(1, 2, 2), c(1, 2, 3, 4), c(1, 1, 2, 2))
  for(alpha in alphas) {
    pred = apply(sapply(alpha, function(alpha_j) {expdiv:::.log_phi_main_func(alpha_j*t, alpha, K)}), 1, sum)
    # plot(t, abs(pred - true), main = paste(alpha, collapse = " "))
    expect_equal(pred, true, tolerance = 1e-4)
  }
})

test_that("log_phi_func verifies the equation with good accuracy in [-1, 1]", {
  t = seq(from = -1, 1, length.out = 1000)
  alpha = c(1, 2)
  K = 1e3
  pred = apply(sapply(alpha, function(alpha_j) {
    log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
    }), 1, sum)
  true = -log(1-1i*t)
  # plot(t, abs(pred - true))
  expect_equal(pred, true, tolerance = 1e-15) # much better convergence

  ## Test with other alpha
  K = 1e3
  t = seq(from = -1, 1, length.out = 1000)
  true = -log(1-1i*t)
  alphas = list(c(1, 3), c(1, 1, 2), c(1, 2, 2), c(1, 2, 3, 4), c(1, 1, 2, 2), c(1), c(2), c(3, 3))
  for(alpha in alphas) {
    pred = apply(sapply(alpha, function(alpha_j) {
      log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
      }), 1, sum)
    # plot(t, abs(pred - true), main = paste(alpha, collapse = " "))
    expect_equal(pred, true, tolerance = 1e-15)
  }
})

test_that("log_phi_func verifies the equation with good accuracy on R", {
  t = seq(from = -10, 10, length.out = 1000)
  alpha = c(1, 2)
  K = 1e3
  pred = apply(sapply(alpha, function(alpha_j) {
    log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
  }), 1, sum)
  true = -log(1-1i*t)
  # plot(t, abs(pred - true))
  expect_equal(pred, true, tolerance = 1e-15) # much better convergence

  ## Test with other alpha
  K = 1e3
  t = seq(from = -3, 3, length.out = 100)
  true = -log(1-1i*t)
  alphas = list(c(1, 3), c(1, 1, 2), c(1, 2, 2), c(1, 2, 3, 4), c(1, 1, 2, 2), c(1), c(2), c(3, 3))
  for(alpha in alphas) {
    pred = apply(sapply(alpha, function(alpha_j) {
      log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
      }), 1, sum)
    # plot(t, abs(pred - true), main = paste(alpha, collapse = " "))
    expect_equal(pred, true, tolerance = 1e-15)
  }
})

test_that("log_phi_func verifies the equation with good accuracy on the radius", {
  t = seq(from = 0.9, 1.1, length.out = 1001)
  t[501] = 1L

  alpha = c(1, 2)
  K = 1e3
  pred = apply(sapply(alpha, function(alpha_j) {
    log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
  }), 1, sum)
  true = -log(1-1i*t)
  # plot(t, abs(pred - true))
  expect_equal(pred, true, tolerance = 1e-15) # much better convergence

  ## Test with other alpha
  K = 1e3
  t = seq(from = 0.9, 1.1, length.out = 1001)
  true = -log(1-1i*t)
  alphas = list(c(1, 3), c(1, 1, 2), c(1, 2, 2), c(1, 2, 3, 4), c(1, 1, 2, 2), c(1), c(2), c(3, 3))
  for(alpha in alphas) {
    pred = apply(sapply(alpha, function(alpha_j) {
      log_phi_func(alpha_j*t, alpha, K, verbose = FALSE)
    }), 1, sum)
    # plot(t, abs(pred - true), main = paste(alpha, collapse = " "))
    expect_equal(pred, true, tolerance = 1e-15)
  }
})
