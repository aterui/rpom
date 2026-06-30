# Generate test data -------------------------------------------------------

## Food web
fwb <- matrix(1, 2, 2)
diag(fwb) <- 0

## Ecosystem parameters
rl <- runif(1, 100, 200)
lambda_b <- runif(1, 0.5, 1)

h <- 1
delta <- runif(1, 0, 0.1)

## Resource and propagule supply
r0 <- runif(1, 0.5, 1)
g <- runif(1, 1, 10)
b <- (1 - r0) / rl

## Extinction rates
mu0 <- runif(1, 0, 0.1)
mu_p <- runif(1, 5, 10)
mu_c <- runif(1, 0.1, 1)

## No spatial synchrony
nu <- 0

## Numerical solution (includes predation)
cout <- npom(
  w = fwb,
  size = rl,
  lambda = lambda_b,
  h = h,
  delta = delta,
  r0 = r0,
  b = b,
  g = g,
  mu0 = mu0,
  mu_p = mu_p,
  mu_c = mu_c,
  nu = nu,
  threshold = 1e-5,
  n_timestep = 200
)

y <- unname(cout[nrow(cout), -1])

## Analytical solution (ignores predation)
y_base <- p_base(
  size = rl,
  lambda = lambda_b,
  h = h,
  delta = delta,
  r0 = r0,
  b = b,
  g = g,
  mu = mu0,
  nu = nu,
  exact = TRUE
)

y_cnsm <- p_cnsm(
  size = rl,
  lambda = lambda_b,
  h = h,
  delta = delta,
  prey = y_base,
  max_prey = sum(fwb[, 2]),
  g = g,
  mu = c(mu0, mu_p),
  nu = nu,
  exact = TRUE
)

y0 <- c(y_base, y_cnsm)

# Tests --------------------------------------------------------------------

test_that("predation reduces basal occupancy", {
  expect_lt(y[1], y0[1])
})
