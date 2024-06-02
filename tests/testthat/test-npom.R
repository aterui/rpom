
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
fwb <- matrix(1, 2, 2)
diag(fwb) <- 0

## - ecosystem structure
rl <- runif(1, 100, 200)
lambda_b <- runif(1, 0.5, 1)
h <- delta <- 1

## - resource and propagules
rsrc <- runif(1, 0.5, 1)
g <- runif(1, 1, 10)

## - extinction rates
mu0 <- runif(1, 0, 0.1)
mu_p <- runif(1, 5, 10)
mu_c <- runif(1, 0.1, 1)

## - synchrony prob.
## - set 0 to remove spatial effect; focus on food webs
rho <- 0

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
             mu_c = mu_c,
             rho = rho,
             threshold = 1e-05,
             n_timestep = 200)

(y <- c(cout[nrow(cout), -1]))
names(y) <- NULL

## analytical
y_base <- p_base(size = rl,
                 lambda = lambda_b,
                 h = h,
                 delta = delta,
                 rsrc = rsrc,
                 g = g,
                 mu = mu0,
                 rho = rho)

y_cnsm <- p_cnsm(size = rl,
                 lambda = lambda_b,
                 h = h,
                 delta = delta,
                 prey = y_base,
                 max_prey = sum(fwb[,2]),
                 g = g,
                 mu = c(mu0, mu_p),
                 rho = rho)

(y0 <- c(y_base, y_cnsm))

# test --------------------------------------------------------------------

test_that("p_base is smaller when predation included", {
  expect_lt(y[1], y0[1])
})

test_that("basal species has greater/less occupancy than consumers", {

  ## extinction to colonization ratio for consumers
  clnz <- rsrc * g * (1 - exp(-delta * h))
  extn <- mu0 + mu_p * (1 - y_base)
  z <- extn / clnz

  if (z < y_base * (1 - y_base)) {
    # when z < y_base * (1 - y_base), basal p < consumer p
    expect_lt(y0[1], y0[2])
  } else {
    # when z > y_base * (1 - y_base), basal p > consumer p
    expect_gt(y0[1], y0[2])
  }

})
