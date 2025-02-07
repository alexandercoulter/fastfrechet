test_that("compatibility checks", {
  expect_error(check_wholenumber(2.2),
    regexp = ".*", # Replace with specific error message if applicable
    info = "Function check_wholenumber should throw an error for input 2.2"
  )

  expect_error(check_wholenumber(-2),
    regexp = ".*", # Replace with specific error message if applicable
    info = "Function check_wholenumber should throw an error for input 2.2"
  )
  expect_error(check_numeric("a"),
    regexp = ".*", # Replace with specific error message if applicable
    info = "Function check_numeruc should throw an error"
  )

  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- mydata$X
  Y <- mydata$Y

  Y[Y > upper] <- upper

  # Estimate conditional QFs
  expect_error(
    frechetreg_univar2wass(
      X = X,
      Y = Y,
      Z = X,
      C_init = NULL,
      lambda = NULL,
      lower = lower,
      upper = lower - 1
    ),
    regexp = ".*", # Replace with specific error message if applicable
    info = "upper is less than lower"
  )

  expect_error(monotoneQP(Y[1, ], lower = lower, upper = lower - 2))
})
