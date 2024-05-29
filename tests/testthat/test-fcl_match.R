
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
n_sp <- round(runif(1, 1, 10))
fwb <- matrix(0, n_sp, n_sp)
for(i in seq_len(nrow(fwb) - 1)) {
  fwb[i, i + 1] <- 1
}

fwb <- fwb + t(fwb)

## - ecosystem structure
L <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 1, 10)
g <- runif(1, 1, 10)

## - extinction rates
mu0 <- runif(1, 0, 1.5)
mu_p <- runif(1, 0, 1.5)

## analytical
y1 <- fcl(foodweb = fwb,
          lambda = lambda_b,
          size = L,
          h = h,
          delta = delta,
          rsrc = rsrc,
          g = g,
          mu0 = mu0,
          mu_p = mu_p,
          mu_s = 0,
          rho = 0)

## numerical
y2 <- nfcl(foodweb = fwb,
           lambda = lambda_b,
           size = L,
           h = h,
           delta = delta,
           rsrc = rsrc,
           g = g,
           mu0 = mu0,
           mu_p = mu_p,
           mu_c = 0,
           mu_s = 0,
           rho = 0)

# test --------------------------------------------------------------------

test_that("fcl() = nfcl()", {
  expect_equal(c(y1), c(y2))
})
