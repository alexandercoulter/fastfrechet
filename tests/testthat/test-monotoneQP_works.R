test_that("monotoneQP output satisfies constraints", {
  # Set parameters
  lower <- 0.5
  upper <- 1.5

  # Generate example vector
  m <- 100
  set.seed(31)
  y <- rnorm(m, 2 * seq(0, 1, len = m), 0.1)

  # Calculate monotone, box-constrained projection
  output <- monotoneQP(y, lower = lower, upper = upper)

  # Check that the output is non-decreasing
  monotonicity_check <- round(min(diff(output$Solution[1, ])), digits = 5)
  expect_true(
    monotonicity_check >= 0,
    sprintf("The Solution is not non-decreasing: min(diff) = %.5f", monotonicity_check)
  )

  # Check lower bound
  lower_bound_diff <- min(output$Solution) - lower
  expect_true(
    lower_bound_diff >= 0,
    sprintf("The Solution violates the lower bound: min(Solution) - lower = %.5f", lower_bound_diff)
  )

  # Check upper bound
  upper_bound_diff <- max(output$Solution) - upper
  expect_true(
    upper_bound_diff <= 0,
    sprintf("The Solution violates the upper bound: max(Solution) - upper = %.5f", upper_bound_diff)
  )
})


test_that("monotoneQP works with random active constraint estimate C_init", {
  # Set parameters
  lower <- 0.5
  upper <- 1.5

  # Generate example vector
  m <- 100
  set.seed(31)
  y <- rnorm(m, 2 * seq(0, 1, len = m), 0.1)

  # Generate random C_init matrix:
  C_init <- matrix(rbinom(m + 1, 1, 0.5), 1, m + 1)

  # Calculate monotone, box-constrained projection
  output1 <- monotoneQP(y, lower = lower, upper = upper, C_init = C_init)
  output2 <- monotoneQP(y, lower = lower, upper = upper, C_init = NULL)

  predicted_QF_diff <- output1$Solution - output2$Solution

  expect_true(max(abs(predicted_QF_diff)) <= 1e-10,
    info = "QF values differ by C_init"
  )
})


test_that("monotoneQP works when given impossible full constraint matrix C_init = 1", {
  # Set parameters
  lower <- 0.5
  upper <- 1.5

  # Generate example vector
  m <- 100
  set.seed(31)
  y <- rnorm(m, 2 * seq(0, 1, len = m), 0.1)

  # Generate random C_init matrix:
  C_init <- matrix(1, 1, m + 1)

  # Calculate monotone, box-constrained projection
  output1 <- monotoneQP(y, lower = lower, upper = upper, C_init = C_init)
  output2 <- monotoneQP(y, lower = lower, upper = upper, C_init = NULL)

  predicted_QF_diff <- output1$Solution - output2$Solution

  expect_true(max(abs(predicted_QF_diff)) <= 1e-10,
    info = "QF values differ by C_init"
  )
})
