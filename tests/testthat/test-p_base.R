# Generate test data -------------------------------------------------------

## Food web
fwb <- matrix(0, 1, 1)

## Ecosystem parameters
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)

h <- 1
delta <- 1

## Resource and propagule supply
r0 <- runif(1, 0, 1)
g <- runif(1, 1, 10)
b <- (1 - r0) / rl

## Extinction rate
mu0 <- runif(1, 0, 0.1)

## Spatial synchrony
## Set nu = 1 / rl to minimize spatial effects and isolate local dynamics.
nu <- 1 / rl
rho0 <- runif(1)

## Numerical solution
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
  mu_c = 0,
  rho0 = rho0,
  nu = nu,
  n_timestep = 300
)

(y0 <- unname(cout[nrow(cout), -1]))

## Analytical solution
(y <- p_base(
  size = rl,
  lambda = lambda_b,
  h = h,
  delta = delta,
  r0 = r0,
  b = b,
  g = g,
  mu = mu0,
  rho0 = rho0,
  nu = nu,
  exact = TRUE
))

# Tests --------------------------------------------------------------------

test_that("p_base() matches the equilibrium occupancy from npom()", {
  expect_equal(y, y0, tolerance = 1e-5)
})
