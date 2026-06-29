test_that("p_mag() matches published values from Shreve (1967)", {
  ## Values compared against Shreve (1967)

  expect_equal(p_mag(1, 1), 1)

  expect_equal(
    p_mag(1, 100),
    0.50251,
    tolerance = 1e-5
  )

  expect_equal(
    p_mag(8, 200),
    0.01336,
    tolerance = 1e-5
  )
})
