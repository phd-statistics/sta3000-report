circulant_symmetric_mat <- function(p) {
  # 1. Initialize a vector of zeros of length p
  v <- rep(0, p)
  
  # 2. Set the diagonal entry (j = k)
  v[1] <- 1
  
  # 3. Set the first set of off-diagonals: k in {j+1, ..., j+5}
  # For the first row (j=1), this is indices 2, 3, 4, 5, 6
  v[2:6] <- 0.1
  
  # 4. Set the wrap-around off-diagonals: k in {j+p-5, ..., j+p-1}
  # For the first row (j=1), this is indices p-4, p-3, p-2, p-1, p
  v[(p-4):p] <- 0.1
  
  # 5. Use the toeplitz function to expand this row into a symmetric matrix
  Sigma <- toeplitz(v)
  
  return(Sigma)
}

compute_metrics <- function(thetahat, ci_lower, ci_upper, pvals, theta0, S, alpha) {
  p <- length(theta0)
  Sc <- setdiff(1:p, S) # Inactive set (S^c)
  s0 <- length(S)
  
  # Interval Lengths (Eq 57 & 58)
  lengths <- ci_upper - ci_lower
  avg_len <- mean(lengths)
  avg_len_S <- mean(lengths[S])
  avg_len_Sc <- mean(lengths[Sc])
  
  # Coverage Flags (1 if true theta is inside CI, 0 otherwise)
  covered <- (theta0 >= ci_lower) & (theta0 <= ci_upper)
  
  # Average Coverage (Eq 59, 60, 61)
  dCov <- mean(covered)
  dCov_S <- mean(covered[S])
  dCov_Sc <- mean(covered[Sc])
  
  # False Positive and True Positive Rates
  rejected <- (pvals < alpha) 
  
  TPR <- sum(rejected[S]) / s0              # Power
  FPR <- sum(rejected[Sc]) / (p - s0)       # Type I Error
  
  # Estimation Metrics (Using thetahat)
  # Handle cases where a method might not return a point estimate
  if (!is.null(thetahat)) {
    bias <- mean(thetahat - theta0)
    mse <- mean((thetahat - theta0)^2)
  } else {
    bias <- NA
    mse <- NA
  }
  
  return(list(
    AvgLength = avg_len,
    AvgLength_S = avg_len_S,
    AvgLength_Sc = avg_len_Sc,
    dCov = dCov,
    dCov_S = dCov_S,
    dCov_Sc = dCov_Sc,
    TPR = TPR,
    FPR = FPR,
    Bias = bias,
    MSE = mse,
    Thetahat = thetahat # Returns the actual vector if you need it
  ))
}