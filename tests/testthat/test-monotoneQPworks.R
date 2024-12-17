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
  expect_true(monotonicity_check >= 0, 
              sprintf("The Solution is not non-decreasing: min(diff) = %.5f", monotonicity_check))
  
  # Check lower bound
  lower_bound_diff <- min(output$Solution) - lower
  expect_true(lower_bound_diff >= 0, 
              sprintf("The Solution violates the lower bound: min(Solution) - lower = %.5f", lower_bound_diff))
  
  # Check upper bound
  upper_bound_diff <- max(output$Solution) - upper
  expect_true(upper_bound_diff <= 0, 
              sprintf("The Solution violates the upper bound: max(Solution) - upper = %.5f", upper_bound_diff))
})


test_that("predicted QFs are monotonic and within [0, 40]", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] = upper
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = NULL,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})


test_that("works with large p", {
  # Parameters
  lower <- 0   # Lower bound
  upper <- 40  # Upper bound
  
  # Generate data
  n <- 50
  p <- 70
  m <- 30
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)
  
  X <- mydata$X
  Y <- mydata$Y
  
  Y[Y > upper] = upper
  
  # Estimate conditional QFs
  output <- frechetreg_univar2wass(X = X,
                                   Y = Y,
                                   Z = X,
                                   C_init = NULL,
                                   lambda = NULL,
                                   lower = lower,
                                   upper = upper)
  
  # Extract predicted QFs
  predicted_QFs <- output$Qhat  # Assuming this contains (n x m) matrix of QFs
  
  # Check monotonicity for each row (QFs should be non-decreasing)
  apply(predicted_QFs, 1, function(row) {
    expect_true(all(diff(row) >= -1e-5), info = "QF row is not monotonic (non-decreasing)")
  })
  
  # Check lower bound for each element
  expect_true(all(predicted_QFs >= lower - 1e-5), 
              info = "QF values violate the lower bound")
  
  # Check upper bound for each element
  expect_true(all(predicted_QFs <= upper + 1e-5), 
              info = "QF values violate the upper bound of 40")
})

