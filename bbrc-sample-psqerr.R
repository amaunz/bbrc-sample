# # # Probabilistic squared error
# Integrates over squaredError(x) * p(x), with X~Poisson(lambda)
# Integral approximation: rectangular in limits from, to
# 
# Params: Integral bounds (from, to)
#         Poisson parameter (lambda)
#         Expected value of x (spx)

# Rectangular integration with step width 1
rectIntegrate = function(f, from, to, ...) { 
  n=seq(floor(from), ceiling(to))
  sum(do.call(f, c(list(x=n, ...))))
}

# Inner integrand; weighted squared error values
# Params: point to estimate (x)
wSquaredErr = function(x, lambda, spx) {
  dpois(x,lambda)*squaredErr(x,spx)
}

# Squared error (Yates correction)
squaredErr = function(x, spx) {
  (x-spx-0.5)^2 / spx
}
