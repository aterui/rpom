
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
fwb <- matrix(1, 2, 2)
diag(fwb) <- 0

## - ecosystem structure
rl <- runif(1, 10, 100)
lambda_b <- runif(1, 0.1, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 1, 10)
g <- runif(1, 1, 10)

## - extinction rates
mu0 <- runif(1, 0, 0.1)
mu_p <- runif(1, 0, 0.1)
mu_s <- runif(1, 0, 0.1)

## - synchrony prob.
## - set 0 to remove spatial effect; focus on food webs
rho <- runif(1, 0, 0.5)

## numerical
cout <- npom(foodweb = fwb,
             size = rl,
             lambda = lambda_b,
             h = h,
             delta = delta,
             rsrc = rsrc,
             g = g,
             mu0 = mu0,
             mu_p = mu_p,
             mu_c = 0,
             mu_s = mu_s,
             rho = rho)

y0 <- c(cout[nrow(cout), -1])
names(y0) <- NULL

## analytical
p1 <- p_base(size = rl,
             lambda = lambda_b,
             h = h,
             delta = delta,
             rsrc = rsrc,
             g = g,
             mu = c(mu0, mu_s),
             rho = rho)

y <- p_cnsm(size = rl,
            lambda = lambda_b,
            h = h,
            delta = delta,
            prey = p1,
            max_prey = 1,
            g = g,
            mu = c(mu0, mu_p, mu_s),
            rho = rho)


# test --------------------------------------------------------------------

test_that("p_base = npom output", {
  expect_equal(round(y, 5), round(y0[2], 5))
})
