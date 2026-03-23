library(mvnfast)
library(hdi)
source("lasso_inference.R")
source("setup.R")

configurations <- list(
  c(n=1000, p=600, s0=10, b=0.5),
  c(n=1000, p=600, s0=10, b=0.25),
  c(n=1000, p=600, s0=10, b=0.1),
  c(n=1000, p=600, s0=30, b=0.5),
  c(n=1000, p=600, s0=30, b=0.25),
  c(n=1000, p=600, s0=30, b=0.1),
  c(n=2000, p=1500, s0=50, b=0.5),
  c(n=2000, p=1500, s0=50, b=0.25),
  c(n=2000, p=1500, s0=50, b=0.1),
  c(n=2000, p=1500, s0=25, b=0.1),
  c(n=2000, p=1500, s0=25, b=0.1),
  c(n=2000, p=1500, s0=25, b=0.1)
)

config_no <- length(configurations)
res <- list(list(NULL), length(configurations)) #Pre-allocated List Vector

n_sims <- 20 # Number of independent realizations
res <- list() # Final results list

for(i in 1:config_no){
  set.seed(i)
  # 1. Extract Configuration Parameters
  dat <- configurations[[i]] 
  n <- dat["n"]
  p <- dat["p"]
  s0 <- dat["s0"]
  b <- dat["b"]
  alpha <- 0.05
  
  Sigma <- circulant_symmetric_mat(p)
  S <- sample(1:p, size = s0, replace = FALSE)
  theta0 <- rep(0, p) 
  theta0[S] <- b 
  
  # Temporary lists to hold the 20 runs for this specific config
  sim_sslasso <- list()
  sim_multi <- list()
  sim_ridge <- list()
  sim_lasso_proj <- list()
  
  # --- INNER LOOP: 20 Realizations of Noise ---
  for(sim in 1:n_sims){
    
    # Generate fresh data for this realization
    X <- rmvn(n, mu=rep(0, p), sigma=Sigma) 
    W <- rnorm(n) # NEW measurement noise!
    y <- X %*% theta0 + W # NEW response!
    
    # Fit Models
    fit_sslasso <- SSLasso(X, y, alpha=alpha)
    fit_multi_splitting <- hdi(X, y, method="multi.split", alpha=alpha)
    fit_ridge_proj <- hdi(X, y, method="ridge.proj", alpha=alpha) 
    fit_lasso_proj <- lasso.proj(X, y, alpha=alpha) 
    
    # Compute Metrics (Assume compute_metrics now returns a 1-row data.frame or named numeric vector)
    # We drop thetahat from the return list for easier averaging, or handle it separately
    sim_sslasso[[sim]] <- as.data.frame(compute_metrics(fit_sslasso$unb.coef, fit_sslasso$low.lim, fit_sslasso$up.lim, fit_sslasso$pvals, theta0, S, alpha))
    sim_multi[[sim]] <- as.data.frame(compute_metrics(fit_multi_splitting$estimate, fit_multi_splitting$ci.lower, fit_multi_splitting$ci.upper, fit_multi_splitting$pval, theta0, S, alpha))
    sim_ridge[[sim]] <- as.data.frame(compute_metrics(fit_ridge_proj$estimate, fit_ridge_proj$ci.lower, fit_ridge_proj$ci.upper, fit_ridge_proj$pval, theta0, S, alpha))
    sim_lasso_proj[[sim]] <- as.data.frame(compute_metrics(fit_lasso_proj$estimate, fit_lasso_proj$ci.lower, fit_lasso_proj$ci.upper, fit_lasso_proj$pval, theta0, S, alpha))
    
    cat(sprintf("  Config %d: Finished simulation %d / %d\n", i, sim, n_sims))
  }
  
  # --- AVERAGE THE 20 REALIZATIONS ---
  # Bind the 20 rows together and calculate the column means
  avg_sslasso <- colMeans(do.call(rbind, sim_sslasso), na.rm = TRUE)
  avg_multi <- colMeans(do.call(rbind, sim_multi), na.rm = TRUE)
  avg_ridge <- colMeans(do.call(rbind, sim_ridge), na.rm = TRUE)
  avg_lasso_proj <- colMeans(do.call(rbind, sim_lasso_proj), na.rm = TRUE)
  
  # Compile final averaged results for this configuration
  sim_res <- list(
    params = list(n = n, p = p, s0 = s0, b = b, alpha = alpha),
    metrics = list(
      sslasso = avg_sslasso,
      multi_split = avg_multi,
      ridge_proj = avg_ridge,
      lasso_proj = avg_lasso_proj
    )
  )
  
  res[[i]] <- sim_res
  cat(sprintf("COMPLETED CONFIGURATION %d / %d\n", i, config_no))
}

saveRDS(res, "res.rds")
