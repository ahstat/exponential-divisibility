library(dplyr)

test_that("fourier_inversion works for the normal and exponential distributions", {
  cf <- function(t) exp(-t^2/2) # the standard normal characteristic function
  T_selected = 100
  N_selected = 100
  verbose = FALSE
  x = seq(from = -3, to = 3, length.out = 101)
  F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
  expect_equal(F_out$cdf, pnorm(x), 2e-4)
  expect_equal(F_out$pdf, dnorm(x), 2e-3)

  cf <- function(t) 1/(1 - 1i*t)
  T_selected = 1000
  N_selected = 10000
  verbose = FALSE
  x = seq(from = -3, to = 3, length.out = 1001)
  F_out = fourier_inversion(x, cf, T_selected, N_selected, TRUE)
  expect_equal(F_out$cdf, pexp(x), 3e-6)
  # Problem of effect around the discontinuity in 0 for the pdf

  x = seq(from = 0.05, to = 3, length.out = 1001)
  F_out = fourier_inversion(x, cf, T_selected, N_selected, TRUE)
  expect_equal(F_out$cdf, pexp(x), 3e-6)
  expect_equal(F_out$pdf, dexp(x), 1e-3)
  # ok for positive values
})

test_that("check absolute error for the normal and exponential distributions as function of T/N", {
  # norm ----
  cf <- function(t) exp(-t^2/2) # the standard normal characteristic function
  verbose = FALSE
  x = seq(from = -3, to = 3, length.out = 101)
  list_df = list()

  Ts = c(1, 10, 100, 1000)
  Ns = c(1, 10, 100, 1000)
  for(i in 1:length(Ts)) {
    for(j in 1:length(Ns)) {
      T_selected = Ts[i]
      N_selected = Ns[j]
      F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
      diff_cdf = max(abs(F_out$cdf - pnorm(x)))
      diff_pdf = max(abs(F_out$pdf - dnorm(x)))
      list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                  N_selected = N_selected,
                                                  diff_cdf = diff_cdf,
                                                  diff_pdf = diff_pdf)
    }
  }
  df = bind_rows(list_df)
  # df %>% filter(T_selected == 1)
  # df %>% filter(T_selected == 10)
  # df %>% filter(T_selected == 100)
  # OK when N = 10*T
  out = df %>% filter(T_selected == 10, N_selected == 100)
  expect_true(out$diff_cdf < 1e-15)
  expect_true(out$diff_pdf < 1e-15)

  # exp ----
  cf <- function(t) 1/(1 - 1i*t)
  verbose = FALSE
  x = seq(from = 0.05, to = 3, length.out = 101)
  list_df = list()

  Ts = c(1, 10, 100, 1e3, 1e4, 1e5)
  for(i in 1:length(Ts)) {
    T_selected = Ts[i]
    N_selected = 10*T_selected
    F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
    diff_cdf = max(abs(F_out$cdf - pexp(x)))
    diff_pdf = max(abs(F_out$pdf - dexp(x)))
    list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                N_selected = N_selected,
                                                diff_cdf = diff_cdf,
                                                diff_pdf = diff_pdf)
  }
  df = bind_rows(list_df)
  # df %>% filter(T_selected == 1)
  # df %>% filter(T_selected == 10)
  # df %>% filter(T_selected == 100)
  # df %>% filter(T_selected == 1000)
  # df %>% filter(T_selected == 1e4)
  df %>% filter(T_selected == 1e5)
  out = df %>% filter(T_selected == 1e5)
  expect_true(out$diff_cdf < 1e-9)
  expect_true(out$diff_pdf < 1e-5)
})

test_that("check absolute error for alpha = c(1, 1) for which closed-form is known", {
  ## (1,1) cdf ----
  alpha = c(1, 1)
  length.out = 10000 # 100000 also tested, but no slower without any improvement
  logTmin = -10
  logTmax = 10
  K = 100
  ts = precalculated_points(length.out, logTmin, logTmax) # complete interval
  cf = precalculated_phi_func(alpha, ts, K)
  verbose = FALSE

  x = seq(from = -3, to = 3, length.out = 101)
  list_df = list()
  cdf_alpha11 = pgamma(x, shape = 1/2, scale = 1/1)
  pdf_alpha11 = dgamma(x, shape = 1/2, scale = 1/1)
  Ts = c(1, 10, 100, 1e3, 1e4)
  for(i in 1:length(Ts)) {
    T_selected = Ts[i]
    # print(T_selected)
    N_selected = 10*T_selected
    F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
    diff_cdf = max(abs(F_out$cdf - cdf_alpha11))
    diff_pdf = max(abs(F_out$pdf - pdf_alpha11))
    list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                N_selected = N_selected,
                                                diff_cdf = diff_cdf,
                                                diff_pdf = diff_pdf
    )
  }
  df = bind_rows(list_df)
  # plot(log(df$T_selected), log(df$diff_cdf), type = "o")
  expect_true(df$diff_cdf[5] < 5e-3)
  # df$diff_pdf issue in 0 with the discontinuity

  ## (1,1) pdf ----
  x = seq(from = 0.05, to = 3, length.out = 101)
  list_df = list()
  cdf_alpha11 = pgamma(x, shape = 1/2, scale = 1/1)
  pdf_alpha11 = dgamma(x, shape = 1/2, scale = 1/1)
  Ts = c(1, 10, 100, 1e3, 1e4)
  for(i in 1:length(Ts)) {
    T_selected = Ts[i]
    # print(T_selected)
    N_selected = 10*T_selected
    F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
    diff_cdf = max(abs(F_out$cdf - cdf_alpha11))
    diff_pdf = max(abs(F_out$pdf - pdf_alpha11))
    list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                N_selected = N_selected,
                                                diff_cdf = diff_cdf,
                                                diff_pdf = diff_pdf
    )
  }
  df = bind_rows(list_df)
  # plot(log(df$T_selected), log(df$diff_cdf), type = "o")
  expect_true(df$diff_cdf[5] < 2e-5) # better without taking 0 inside

  # plot(log(df$T_selected), log(df$diff_pdf), type = "o")
  expect_true(df$diff_pdf[5] < 3e-2) # a little better without taking 0 inside
})

test_that("check absolute error for alpha = c(2, 2) for which closed-form is known", {
  ## (2,2) cdf ----
  alpha = c(2, 2)
  length.out = 10000
  logTmin = -10
  logTmax = 10
  K = 100
  ts = precalculated_points(length.out, logTmin, logTmax) # complete interval
  cf = precalculated_phi_func(alpha, ts, K)

  verbose = FALSE
  x = seq(from = -3, to = 3, length.out = 101)
  list_df = list()

  cdf_alpha22 = pgamma(x, shape = 1/2, scale = 1/2)
  pdf_alpha22 = dgamma(x, shape = 1/2, scale = 1/2)

  Ts = c(1, 10, 100, 1e3, 1e4)
  for(i in 1:length(Ts)) {
    T_selected = Ts[i]
    # print(T_selected)
    N_selected = 10*T_selected
    F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
    diff_cdf = max(abs(F_out$cdf - cdf_alpha22))
    diff_pdf = max(abs(F_out$pdf - pdf_alpha22))
    list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                N_selected = N_selected,
                                                diff_cdf = diff_cdf,
                                                diff_pdf = diff_pdf
    )
  }
  df = bind_rows(list_df)
  plot(log(df$T_selected), log(df$diff_cdf), type = "o")
  expect_true(df$diff_cdf[5] < 7e-3)
  # df$diff_pdf issue in 0 with the discontinuity

  ## (2,2) pdf ----
  x = seq(from = 0.05, to = 3, length.out = 101)
  cdf_alpha22 = pgamma(x, shape = 1/2, scale = 1/2)
  pdf_alpha22 = dgamma(x, shape = 1/2, scale = 1/2)
  list_df = list()

  Ts = c(1, 10, 100, 1e3, 1e4)
  for(i in 1:length(Ts)) {
    T_selected = Ts[i]
    # print(T_selected)
    N_selected = 10*T_selected
    F_out = fourier_inversion(x, cf, T_selected, N_selected, verbose)
    diff_cdf = max(abs(F_out$cdf - cdf_alpha22))
    diff_pdf = max(abs(F_out$pdf - pdf_alpha22))
    list_df[[length(list_df) + 1]] = data.frame(T_selected = T_selected,
                                                N_selected = N_selected,
                                                diff_cdf = diff_cdf,
                                                diff_pdf = diff_pdf
    )
  }
  df = bind_rows(list_df)
  # plot(log(df$T_selected), log(df$diff_cdf), type = "o")
  expect_true(df$diff_cdf[5] < 2e-5) # better without taking 0 inside

  # plot(log(df$T_selected), log(df$diff_pdf), type = "o")
  expect_true(df$diff_pdf[5] < 4e-2) # a little better without taking 0 inside
})
