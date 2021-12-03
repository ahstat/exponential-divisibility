rm(list = ls())
library(expdiv)

# In unit tests (before with alpha)
# Step 2: Check for alpha = c(1, 2)

alpha = c(1, 2)
length.out = 10000 # 100000 also feasible
logTmin = -10
logTmax = 10
K = 100
# Key is to compute well on the most important part on [-exp(logTmax), exp(logTmax)]
# which is also the part which is quick, and to discard on outside.
# Also, we further focus around zero with logTmin.

ts = precalculated_points(length.out, logTmin, logTmax) # complete interval
#ts = seq(from = 0.75, to = 22, length.out = 1000) # specific interval 1
cf_main = precalculated_phi_func(alpha, ts, K)
cf_equalMIN = precalculated_phi_func(rep(min(alpha), length(alpha)), ts, K)
cf_equalMAX = precalculated_phi_func(rep(max(alpha), length(alpha)), ts, K)
cf_equalMEAN = precalculated_phi_func(rep(mean(alpha), length(alpha)), ts, K)
cf_equalGEOM = precalculated_phi_func(rep(prod(alpha)^(1/length(alpha)), length(alpha)), ts, K) # very good for t large, equivalent for t large?

# t = seq(from = -1000, to = 1000, length.out = 10000)
t = ts

# Complex characteristic function ----
plot(cf_main(t), type = "l")
lines(cf_equalMIN(t), type = "l", col = "blue")
lines(cf_equalMAX(t), type = "l", col = "red")
lines(cf_equalGEOM(t), type = "l", col = "green")
lines(cf_equalMEAN(t), type = "l", col = "orange")
text(Re(cf_main(t)), Im(cf_main(t)), round(t, 2))

# Real part of the characteristic function ----
plot(t, Re(cf_main(t)), type = "l")
lines(t, Re(cf_equalMIN(t)), type = "l", col = "blue")
lines(t, Re(cf_equalMAX(t)), type = "l", col = "red")
lines(t, Re(cf_equalGEOM(t)), type = "l", col = "green") # perfect for t large
lines(t, Re(cf_equalMEAN(t)), type = "l", col = "orange")
# lines(t, Re(exp(log_phi_func(t, c(1.55, 1.55), K, verbose = FALSE))), type = "l", col = "orange") # perfect for t small

# Imaginary part of the characteristic function ----
plot(t, Im(cf_main(t)), type = "l")
lines(t, Im(cf_equalMIN(t)), type = "l", col = "blue")
lines(t, Im(cf_equalMAX(t)), type = "l", col = "red")
lines(t, Im(cf_equalGEOM(t)), type = "l", col = "green") # perfect for t large
lines(t, Im(cf_equalMEAN(t)), type = "l", col = "orange") # perfect for t small
# lines(t, Im(exp(log_phi_func(t, c(1.55, 1.55), K, verbose = FALSE))), type = "l", col = "orange")

# Distance to 0 of the characteristic function ----
plot(t, abs(cf_main(t)), type = "l")
lines(t, abs(cf_equalMIN(t)), type = "l", col = "blue")
lines(t, abs(cf_equalMAX(t)), type = "l", col = "red")
lines(t, abs(cf_equalGEOM(t)), type = "l", col = "green") # perfect for t large
lines(t, abs(cf_equalMEAN(t)), type = "l", col = "orange")
# lines(t, abs(exp(log_phi_func(t, c(1.59, 1.59), K, verbose = FALSE))), type = "l", col = "orange") # perfect for t small

# Checking of the characteristic function for large t ----
t = seq(from = -30000, to = -29000, length.out = 10)
out = apply(sapply(alpha, function(alpha_j) {log_phi_func(alpha_j*t, alpha, K)}), 1, sum)
plot(out, type = "l")
true = -log(1-1i*t)
plot(t, abs(out - true))

x = 1
alpha = c(1, 2)
K = 100
verbose = FALSE
rel.tol = 1e-4

start.time <- Sys.time()
print("to do a computation")
print(Sys.time() - start.time)



## Next ----
alpha = c(1, 2)
K = 10000
x = seq(from = -3, to = 3, length.out = 101)

cf = not_precalculated_phi_func(alpha, K)
F_out = fourier_inversion(x, cf, T_selected = 100, N_selected = 100, verbose = TRUE)

x = seq(from = -10, to = 10, length.out = 1000)
cdf_alpha22 = pgamma(x, shape = 1/2, scale = 1/2)
cdf_alpha11 = pgamma(x, shape = 1/2, scale = 1/1)
lines(x, cdf_alpha11, col = "black")
lines(x, cdf_alpha22, col = "green")

plot(result$x, result$pdf, type = "l")

# ok seems good

lines(result$x, pgamma(result$x, shape = 1/2, scale = 1), type = "l", col = "blue")
lines(result$x, pgamma(result$x, shape = 1/2, scale = 1/2), type = "l", col = "blue")
# X + Y = Z --> Gamma(1/2, 1)

https://stats.stackexchange.com/questions/358524/sampling-from-characteristic-moment-generating-function
https://stats.stackexchange.com/questions/88697/sample-from-a-custom-continuous-distribution-in-r

N = 1e5
hist(rexp(N, 2), breaks = 1000, probability = TRUE)
hist(rexp(N, 1)/2, breaks = 1000, probability = TRUE)

## To recycle 1 ----

out_fakemiddle_cdf_true = pgamma(x, shape = 1/2, scale = 1/sqrt(prod(alpha))) # 0.907391
plot(out_fakemiddle$cdf - out_fakemiddle_cdf_true, type = "l")
plot(abs(out_fakemiddle$cdf - out_fakemiddle_cdf_true), type = "l", log = "y")

out_fakemiddle2$cdf # 0.9139234
pgamma(x, shape = 1/2, scale = 1/mean(alpha)) # 0.9167355

out_1$cdf # 0.8266215
pgamma(x, shape = 1/2, scale = 1) # 0.8427008

out_2$cdf # 0.9539599
pgamma(x, shape = 1/2, scale = 1/2) # 0.9544997



plot(out$x, out$cdf, type = "l")
lines(out_1$x, out_1$cdf, type = "l", col = "red")
lines(out_2$x, out_2$cdf, type = "l", col = "blue")

plot(out_fakemiddle$x, out_fakemiddle$cdf, type = "l", col = "purple")
lines(x, pgamma(x, shape = 1/2, scale = 1/sqrt(prod(alpha))), col = "blue")
# Mean = k*theta = 1/(2*sqrt(2))
# Sd = sqrt(k)*theta = (1/sqrt(2))*(1/sqrt(2)) = 1/2

plot(out_fakemiddle2$x, out_fakemiddle2$cdf, type = "l", col = "purple")
lines(x, pgamma(x, shape = 1/2, scale = 1/mean(alpha)), col = "blue")
# Mean = k*theta = 1/(2*(3/2)) = 1/3, same as the true one
# Sd = sqrt(k)*theta = (1/sqrt(2))*(sqrt(2)/sqrt(3)) = 1/sqrt(3) = sqrt(3)/3

# But for the true, it's 1/sqrt(5)

# It's like one approximation is good for small t and the other for large t
# true?
plot(out$x, out$cdf, type = "l")
lines(x, pgamma(x, shape = 1/2, scale = 1/sqrt(prod(alpha))), col = "blue")
lines(x, pgamma(x, shape = 1/2, scale = 1/mean(alpha)), col = "blue")

plot(out$x, out$pdf, type = "l", log = "y", xlim = c(0.4, 1), ylim = c(0.1, 1))
lines(x, dgamma(x, shape = 1/2, scale = 1/sqrt(prod(alpha))), col = "blue")
lines(x, dgamma(x, shape = 1/2, scale = 1/mean(alpha)), col = "blue")
