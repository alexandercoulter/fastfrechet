test_that("Shah and Samworth threshold is bounded and monotone in theta", {
  # Set parameters
  p <- 30
  q <- 1:29
  B <- 50
  E_thr <- 1
  
  # Calculate Shah and Samworth thresholds:
  output <- Shah_Samworth_thresholds(p = p, q = q, B = B, E_thr = E_thr)
  
  # Check that given thresholds are monotone in relative model size q / p:
  monotonicity_check <- round(min(diff(output$pointwise_thresholds[order(output$relative_model_size)])), 5)
  expect_true(monotonicity_check >= 0,
              sprintf("The Thresholds are not non-decreasing in theta: min(diff) = %.5f", monotonicity_check))
  
  # Check that Shah and Samworth thresholds are bounded in (0, 1]:
  bounded_check <- range(output$pointwise_thresholds)
  expect_true(bounded_check[1] > 0 & bounded_check[2] <= 1,
              sprintf("The Thresholds are not bounded in (0, 1]: min(thr) = %.5f, max(thr) = %.5f", bounded_check[1], bounded_check[2]))
  
})