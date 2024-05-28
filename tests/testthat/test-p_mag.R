test_that("p_mag works", {
  ## output values are compared with values in Shreve 1967
  expect_equal(p_mag(1, 1), 1)

  expect_equal(round(p_mag(1, 100), 5),
               0.50251)

  expect_equal(round(p_mag(8, 200), 5),
               0.01336)

  expect_equal(p_mag(1, Inf), 0.5)
})
