library(mvnfast)
library(hdi)
library(parallel)
library(doParallel)
library(foreach)

source("lasso_inference.R")
source("setup.R")

# 1. SETUP CLUSTER
cores_to_use <-2
cl <- makeCluster(cores_to_use, outfile = "") 
registerDoParallel(cl)

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
res <- vector("list", config_no) 
n_sims <- 20
alpha <- 0.05
# Critical value for 95% CI construction for the Projection methods
z_crit <- qnorm(1 - alpha/2)

for(i in 1:config_no){
  dat <- configurations[[i]] 
  n <- dat["n"]; p <- dat["p"]; s0 <- dat["s0"]; b <- dat["b"]
  
  Sigma <- circulant_symmetric_mat(p)
  S <- sample(1:p, size = s0, replace = FALSE)
  theta0 <- rep(0, p) 
  theta0[S] <- b 
  
  cat(sprintf("\n--- Starting Config %d (n=%d, p=%d) ---\n", i, n, p))
  
  sim_results <- foreach(sim = 1:n_sims, 
                         .packages = c("mvnfast", "hdi"),
                         .export = c("SSLasso", "compute_metrics")) %dopar% {
                           
                           source("lasso_inference.R", local = TRUE)
                           source("setup.R", local = TRUE)
                           set.seed(i * 1000 + sim)
                           
                           X <- rmvn(n, mu = rep(0, p), sigma = Sigma)
                           W <- rnorm(n) 
                           y <- as.numeric(X %*% theta0 + W) 
                           
                           fit_ms <- multi.split(X, y)
                           
                           compute_metrics(NULL, fit_ms$lci, fit_ms$uci, fit_ms$pval.corr, theta0, S, alpha)
                         }
  
  res[[i]] <- list(
    params = list(n = n, p = p, s0 = s0, b = b),
    metrics = list(
      multi_split = colMeans(do.call(rbind, sim_results), na.rm = TRUE)
    )
  )
  
  saveRDS(res, "res_commute_backup_ms.rds")
}

stopCluster(cl)
saveRDS(res, "res_ms.rds")