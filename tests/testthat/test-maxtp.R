
# data generation ---------------------------------------------------------

## set parameters
## - foodweb
n_sp <- 3
fwb <- matrix(1, n_sp, n_sp)
diag(fwb) <- 0
p_hat <- runif(n_sp, 0, 1)

## manual calculation
wtp0 <- 1 * p_hat[1] / sum(p_hat[1:2]) + 2 * p_hat[2] / sum(p_hat[1:2]) + 1

wtp <- maxtp(fwb, p_hat, weight = TRUE)
uwtp <- maxtp(fwb, p_hat, weight = FALSE)


# test --------------------------------------------------------------------

test_that("maxtp() works", {
  expect_equal(c(wtp), wtp0)
  expect_equal(c(uwtp), 2.5)
})
