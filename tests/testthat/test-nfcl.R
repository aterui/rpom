# Generate test data -------------------------------------------------------

## Food web
n_sp <- round(runif(1, 1, 10))

fwb <- matrix(0, n_sp, n_sp)
for (i in seq_len(n_sp - 1)) {
  fwb[i, i + 1] <- 1
}
fwb <- fwb + t(fwb)

## Ecosystem parameters
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)

h <- 1
delta <- runif(1, 0, 0.1)

## Resource and propagule supply
r0 <- runif(1, 0, 1)
g <- runif(1, 1, 10)
b <- (1 - r0) / rl

## Numerical solution
y <- nfcl(
  w = fwb,
  lambda = lambda_b,
  size = rl,
  h = h,
  delta = delta,
  r0 = r0,
  b = b,
  g = g,
  mu0 = 0,
  mu_p = 0,
  mu_c = 0,
  rho0 = 0,
  nu = 0,
  exact = TRUE
)

# Tests --------------------------------------------------------------------

test_that("nfcl() returns the correct food chain length for a linear food web", {
  expect_equal(as.numeric(y), n_sp)
})
