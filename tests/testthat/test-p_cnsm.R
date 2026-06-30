# Generate test data -------------------------------------------------------

## Food web
fwb <- matrix(1, 2, 2)
diag(fwb) <- 0

## Ecosystem parameters
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)

h <- 1
delta <- runif(1, 0, 0.1)

## Resource and propagule supply
r0 <- runif(1, 0, 1)
g <- runif(1, 1, 10)
b <- (1 - r0) / rl

## Extinction rates
mu0 <- runif(1, 0, 0.1)
mu_p <- runif(1, 0, 0.1)

## Spatial synchrony
## Set nu = 1 / rl to minimize spatial effects and isolate food-web dynamics.
nu <- 1 / rl

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
  mu_p = mu_p,
  mu_c = 0,
  nu = nu
)

y0 <- unname(cout[nrow(cout), -1])

## Analytical solution
p1 <- p_base(
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

y <- p_cnsm(
  size = rl,
  lambda = lambda_b,
  h = h,
  delta = delta,
  prey = p1,
  max_prey = 1,
  g = g,
  mu = c(mu0, mu_p),
  nu = nu,
  exact = TRUE
)

# Tests --------------------------------------------------------------------

test_that("p_cnsm() matches the equilibrium occupancy from npom()", {
  expect_equal(y, y0[2], tolerance = 1e-5)
})
