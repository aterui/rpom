
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
r0 <- runif(1, 0, 1)
g <- runif(1, 1, 10)
b <- (1 / rl) * (1 - r0)

## - extinction rates
mu0 <- runif(1, 0, 0.1)
mu_p <- runif(1, 0, 0.1)

## - synchrony prob.
## - set 0 to remove spatial effect; focus on food webs
nu <- 1 / rl

## numerical
cout <- npom(w = fwb,
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
             nu = nu)

y0 <- c(cout[nrow(cout), -1])
names(y0) <- NULL

## analytical
p1 <- p_base(size = rl,
             lambda = lambda_b,
             h = h,
             delta = delta,
             r0 = r0,
             b = b,
             g = g,
             mu = mu0,
             nu = nu,
             exact = TRUE)

y <- p_cnsm(size = rl,
            lambda = lambda_b,
            h = h,
            delta = delta,
            prey = p1,
            max_prey = 1,
            g = g,
            mu = c(mu0, mu_p),
            nu = nu,
            exact = TRUE)


# test --------------------------------------------------------------------

test_that("p_cnsm = npom output", {
  expect_equal(round(y, 5), round(y0[2], 5))
})
