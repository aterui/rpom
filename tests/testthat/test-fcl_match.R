
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
n_sp <- round(runif(1, 1, 10))
fwbl <- matrix(0, n_sp, n_sp)
for (i in seq_len(nrow(fwbl) - 1)) {
  fwbl[i, i + 1] <- 1
}

fwbl <- fwbl + t(fwbl)

fwb <- matrix(1, 3, 3)
diag(fwb) <- 0

## - ecosystem structure
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 0, 1)
g <- runif(1, 1, 10)

## - extinction rates
mu0 <- runif(1, 0, 1.5)
mu_p <- runif(1, 0, 1.5)

# test --------------------------------------------------------------------

test_that("fcl() = nfcl() with linear food web", {

  ## analytical
  y1 <- fcl(foodweb = fwbl,
            lambda = lambda_b,
            size = rl,
            h = h,
            delta = delta,
            rsrc = rsrc,
            g = g,
            mu0 = mu0,
            mu_p = mu_p,
            rho = 0,
            weight = TRUE)

  ## numerical
  y2 <- nfcl(foodweb = fwbl,
             lambda = lambda_b,
             size = rl,
             h = h,
             delta = delta,
             rsrc = rsrc,
             g = g,
             mu0 = mu0,
             mu_p = mu_p,
             mu_c = 0,
             rho = 0,
             n_timestep = 250,
             weight = TRUE,
             threshold = 1e-5)

  expect_equal(c(round(y1, 5)),
               c(round(y2, 5)))
})

test_that("fcl() = nfcl() with omnivory", {

  ## analytical
  y1 <- fcl(foodweb = fwb,
            lambda = lambda_b,
            size = rl,
            h = h,
            delta = delta,
            rsrc = rsrc,
            g = g,
            mu0 = mu0,
            mu_p = mu_p,
            rho = 0,
            weight = TRUE)

  ## numerical
  y2 <- nfcl(foodweb = fwb,
             lambda = lambda_b,
             size = rl,
             h = h,
             delta = delta,
             rsrc = rsrc,
             g = g,
             mu0 = mu0,
             mu_p = mu_p,
             mu_c = 0,
             rho = 0,
             n_timestep = 250,
             weight = TRUE,
             threshold = 1e-5)

  expect_equal(c(round(y1, 5)),
               c(round(y2, 5)))
})
