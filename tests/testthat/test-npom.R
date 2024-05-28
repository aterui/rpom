
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
fwb <- matrix(1, 2, 2)
diag(fwb) <- 0

## - ecosystem structure
L <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 1, 10)
g <- runif(1, 1, 10)

## - extinction rates
mu0 <- runif(1, 0, 1)
mu_p <- runif(1, 0, 0)
mu_c <- runif(1, 0, mu0 * 0.1)
mu_s <- runif(1, 0, 1)

## - synchrony prob.
## - set 0 to remove spatial effect; focus on food webs
rho <- 0

## numerical
cout <- npom(foodweb = fwb,
             size = L,
             lambda = lambda_b,
             h = h,
             delta = delta,
             rsrc = rsrc,
             g = g,
             mu0 = mu0,
             mu_p = mu_p,
             mu_c = mu_c,
             mu_s = mu_s,
             rho = rho)

print(y0 <- c(cout[nrow(cout), -1]))
names(y0) <- NULL

## analytical
y <- p_base(size = L,
            lambda = lambda_b,
            h = h,
            delta = delta,
            rsrc = rsrc,
            g = g,
            mu = c(mu0, mu_s),
            rho = rho)


# test --------------------------------------------------------------------

test_that("p_base is smaller when predation included", {
  expect_lt(y0[1], y)
})

test_that("basal species has greater occupancy than consumers", {
  expect_lt(y0[2], y0[1])
})
