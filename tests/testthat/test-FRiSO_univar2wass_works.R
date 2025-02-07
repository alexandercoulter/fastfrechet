test_that("works with sequence of tau's", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- mydata$X
  Y <- mydata$Y

  Y[Y > upper] <- upper

  tauseq <- 1:10

  # Calculate lambda-hat values
  output <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = NULL,
    eps = 1e-5,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 1
  )

  # Check all lambda values are non-negative:
  expect_true(all(output >= 0),
    info = "lambda values violate non-negativity"
  )

  # Check lambda vectors add to tauseq values:
  expect_true(all(max(abs(colSums(output) - tauseq)) <= 1e-10),
    info = "lambda values do not add to tau"
  )
})


test_that("works with sequence of tau's and large p", {
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

  tauseq <- 1:10

  # Calculate lambda-hat values
  output <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = NULL,
    eps = 1e-5,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 1
  )

  # Check all lambda values are non-negative:
  expect_true(all(output >= 0),
    info = "lambda values violate non-negativity"
  )

  # Check lambda vectors add to tauseq values:
  expect_true(all(max(abs(colSums(output) - tauseq)) <= 1e-10),
    info = "lambda values do not add to tau"
  )
})


test_that("works with sequence of tau's and momentum (impulse < 1)", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- mydata$X
  Y <- mydata$Y

  Y[Y > upper] <- upper

  tauseq <- 1:10


  # Calculate lambda-hat values with momentum
  output1 <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = NULL,
    eps = 1e-6,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 0.8
  )
  # Calculate lambda-hat values without momentum
  output2 <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = NULL,
    eps = 1e-6,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 1
  )

  # Check all lambda values are non-negative:
  expect_true(all(output1 >= 0),
    info = "lambda values violate non-negativity"
  )

  # Check lambda vectors add to tauseq values:
  expect_true(all(max(abs(colSums(output1) - tauseq)) <= 1e-10),
    info = "lambda values do not add to tau"
  )

  # Check fitted lambda values are approximately equal:
  expect_true(all(output1 >= 0),
    info = "lambda values violate non-negativity"
  )

  # Check lambda vectors do not differ too much based on hyperparameters:
  expect_true(max(abs(output1 - output2)) <= 1e-5,
    info = "lambda values change based on hyperparameters"
  )
})


test_that("works with sequence of tau's and a given lambda_init", {
  # Parameters
  lower <- 0 # Lower bound
  upper <- 40 # Upper bound

  # Generate data
  n <- 100
  p <- 10
  m <- 100
  set.seed(31)
  mydata <- generate_zinbinom_qf(n = n, p = p, m = m)

  X <- mydata$X
  Y <- mydata$Y

  Y[Y > upper] <- upper

  tauseq <- 1:10

  # Calculate lambda-hat values with given lambda_init
  output1 <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = rexp(p),
    eps = 1e-6,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 1
  )
  # Calculate lambda-hat values with no lambda_init
  output2 <- FRiSO_univar2wass(
    X = X,
    Y = Y,
    lower = lower,
    upper = upper,
    tauseq = tauseq,
    lambda_init = NULL,
    eps = 1e-6,
    nudge = 0.001,
    alpha = 0.9,
    max_iter = 1000,
    max_theta = pi / 4,
    impulse = 1
  )

  # Check all lambda values are non-negative:
  expect_true(all(output1 >= 0),
    info = "lambda values violate non-negativity"
  )

  # Check lambda vectors add to tauseq values:
  expect_true(all(max(abs(colSums(output1) - tauseq)) <= 1e-10),
    info = "lambda values do not add to tau"
  )

  # Check lambda vectors do not differ too much based on hyperparameters:
  expect_true(max(abs(output1 - output2)) <= 1e-5,
    info = "lambda values change based on hyperparameters"
  )
})
