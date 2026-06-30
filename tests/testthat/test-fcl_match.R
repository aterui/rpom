# Generate test data -------------------------------------------------------

## Food webs
n_sp <- round(runif(1, 1, 10))

fw_linear <- matrix(0, n_sp, n_sp)
for (i in seq_len(n_sp - 1)) {
  fw_linear[i, i + 1] <- 1
}
fw_linear <- fw_linear + t(fw_linear)

fw_omnivory <- matrix(1, 3, 3)
diag(fw_omnivory) <- 0

## Ecosystem parameters
rl <- runif(1, 100, 1000)
lambda_b <- runif(1, 0.3, 1)

h <- 1
delta <- 1

## Resource and propagule supply
r0 <- runif(1, 0.5, 1)
g <- runif(1, 10, 100)
b <- (1 - r0) / rl

## Extinction rates
mu0 <- runif(1, 0, 1.5)
mu_p <- runif(1, 0, 1.5)

## Synchrony
nu <- 1 / rl
rho0 <- runif(1)

check_fcl <- function(w) {

  analytical <- fcl(
    w = w,
    lambda = lambda_b,
    size = rl,
    h = h,
    delta = delta,
    r0 = r0,
    b = b,
    g = g,
    mu0 = mu0,
    mu_p = mu_p,
    rho0 = rho0,
    nu = nu,
    weight = TRUE,
    exact = TRUE
  )

  numerical <- nfcl(
    w = w,
    lambda = lambda_b,
    size = rl,
    h = h,
    delta = delta,
    r0 = r0,
    b = b,
    g = g,
    mu0 = mu0,
    mu_p = mu_p,
    mu_c = 0,
    rho0 = rho0,
    nu = nu,
    n_timestep = 250,
    threshold = 1e-5,
    weight = TRUE
  )

  print(list(analytical, numerical))
  attributes(analytical) <- attributes(numerical) <- NULL

  expect_equal(analytical, numerical, tolerance = 1e-4)
}

# Tests --------------------------------------------------------------------

test_that("fcl() matches nfcl() for a linear food web", {
  check_fcl(fw_linear)
})

test_that("fcl() matches nfcl() for an omnivorous food web", {
  check_fcl(fw_omnivory)
})
