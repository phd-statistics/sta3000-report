library(mvnfast)
library(hdi)
source("lasso_inference.R")
source("setup.R")

data("riboflavin")
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)

X<- riboflavin$x
y<- riboflavin$y
fit.multi <- multi.split(X,y)
fit.ridge.proj <- ridge.proj(X, y)
fit.sslasso <- SSLasso(X,y)
fit.lasso.proj <- lasso.proj(X,y)

# 1. SS Lasso (Uses unb.coef, low.lim, up.lim, pvals)
m_ss <- compute_metrics(fit.sslasso$unb.coef, fit.sslasso$low.lim, fit.sslasso$up.lim, fit.sslasso$pvals, theta0, S, alpha)

# 2. Multi-Split (Uses estimate [if exists], lci, uci, pval)
# Note: multi.split often doesn't return a point estimate. We pass NULL if missing.
m_ms <- compute_metrics(NULL, fit.multi$lci, fit.multi$uci, fit.multi$pval, theta0, S, alpha)

# 3. Ridge Proj (Uses bhat, pval, and SE to build intervals)
ci_low_rp <- fit.ridge.proj$bhat - z_crit * fit_rp$se
ci_up_rp  <- fit.ridge.proj$bhat + z_crit * fit_rp$se
m_rp <- compute_metrics(fit.ridge.proj$bhat, ci_low_rp, ci_up_rp, fit.ridge.proj$pval, theta0, S, alpha)

# 4. Lasso Proj (Uses bhat, pval, and SE to build intervals)
ci_low_lp <- fit.lasso.proj$bhat - z_crit * fit_lp$se
ci_up_lp  <- fit.lasso.proj$bhat + z_crit * fit_lp$se
m_lp <- compute_metrics(fit.lasso.proj$bhat, ci_low_lp, ci_up_lp, fit.lasso.proj$pval, theta0, S, alpha)

list(sslasso = m_ss, multi = m_ms, ridge = m_rp, lasso_proj = m_lp)