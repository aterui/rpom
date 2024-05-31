
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
n_sp <- round(runif(1, 1, 10))
fwb <- matrix(0, n_sp, n_sp)
for (i in seq_len(nrow(fwb) - 1)) {
  fwb[i, i + 1] <- 1
}

fwb <- fwb + t(fwb)

## - ecosystem structure
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 1, 10)
g <- runif(1, 1, 10)

## analytical
y <- fcl(foodweb = fwb,
         lambda = lambda_b,
         size = rl,
         h = h,
         delta = delta,
         rsrc = rsrc,
         g = g,
         mu0 = 0,
         mu_p = 0,
         rho = 0)

# test --------------------------------------------------------------------

test_that("fcl() works", {
  expect_equal(c(y), n_sp)
})
