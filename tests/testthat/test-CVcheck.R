test_that("Support matched lambda in CV", {
  # Generate data for X and Y inputs by using the output of `generate_zinbinom_qf`
  n = 100  # number of samples - nrow(X) and nrow(Y).
  p = 10   # number of covariates - ncol(X).
  m = 30   # EQF grid density - ncol(Y).
  lower = 0
  upper = Inf
  
  set.seed(31)
  mydata = generate_zinbinom_qf(n = n, p = p, m = m)
  X = mydata$X  # (n x p) matrix of covariates
  Y = mydata$Y  # (n x m) matrix of EQFs, stored row-wise
  
  # Set cross-validation parameters
  K = 5
  thresh = 0.0001
  tauseq = seq(0.1, 20, 0.5)
  eps = 0.01
  
  # Run complementary pairs stability selection
  cv = FRiSO_CV_univar2wass(X = X,
                            Y = Y,
                            K = K,
                            thresh = thresh,
                            lower = lower,
                            upper = upper,
                            tauseq = tauseq,
                            eps = eps)
  # Extract results
  opt_lambda <- cv$opt_lambda  # Vector of lambda values
  opt_selected <- cv$opt_selected  # Indices selected
  
  # Find positions where opt_lambda > thresh
  valid_positions <- which(opt_lambda > thresh)
  
  # Check that opt_selected matches the valid positions
  expect_equal(opt_selected, valid_positions, 
               info = "opt_selected indices do not match positions where opt_lambda > thresh")
}) 
  
