# Generate test data -------------------------------------------------------

## Food web
n_sp <- 3
fwb <- matrix(1, n_sp, n_sp)
diag(fwb) <- 0

## Species occupancies
p_hat <- runif(n_sp)

## Expected weighted maximum trophic position
wtp_expected <-
  1 +
  p_hat[1] / sum(p_hat[1:2]) +
  2 * p_hat[2] / sum(p_hat[1:2])

## Computed values
wtp <- maxtp(fwb, p_hat, weight = TRUE)
uwtp <- maxtp(fwb, p_hat, weight = FALSE)

# Tests --------------------------------------------------------------------

test_that("maxtp() computes weighted and unweighted trophic positions", {
  expect_equal(as.numeric(wtp), wtp_expected)
  expect_equal(as.numeric(uwtp), 2.5)
})
