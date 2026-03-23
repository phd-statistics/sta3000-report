set.seed(3000)
library(mvnfast)
library(hdi)
library(ggplot2)
source("lasso_inference.R")
source("setup.R")
# A replication of the simulation study used in the paper

####################
# 5.1 Reproduction #
####################

#########
# Setup #
#########

# Setting up Dimensions
n <- 1000
p <- 600
alpha <- 0.05
Sigma <- circulant_symmetric_mat(p) # Covariance Matrix
mu <- rep(0,p)
# True Model Specifications
s0 <- 10 #Sparsity 
b <- 1 # Signal Strength

S <- sample(1:p, size = s0, replace = FALSE) # supp(theta0)
theta0 <- matrix(rep(0, p)) # Initialize theta0 as a vector of zeros
theta0[S,] <- b # Assign the value of b to the indices in supp(theta0)

# Data
X<- mvnfast::rmvn(n,mu=mu, sigma=Sigma) #Design Matrix
W <- rnorm(n) # Noise
y <- X%*%theta0 + W #Response

#################
# Model Fitting #
#################

# Paper Shows these

# SSLASSO (Has LASSO in it)
fit_sslasso<-SSLasso(X, y)
# Multisample Splitting
fit_multi_splitting <- hdi(X, y, method= "multi.split")
# Ridge Projection Estimator
fit_ridge_proj <- ridge.proj(X,y)

# Paper DOES NOT Show these (too similar?)
fit_lasso_proj <- lasso.proj(X,y)
fit_lasso <- boot.lasso.proj(X,y)

####################################
# Visualize the Covariance Matrix  #
####################################

# Convert matrix to a long data frame
Sigma_df <- as.data.frame(as.table(Sigma))
colnames(Sigma_df) <- c("Row", "Col", "Value")
Sigma_df$Row <- as.numeric(Sigma_df$Row)
Sigma_df$Col <- as.numeric(Sigma_df$Col)

ggplot(Sigma_df, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  # Using 'magma' or 'viridis' helps differentiate small values from zero
  scale_fill_viridis_c(option = "magma", direction = -1) + 
  scale_y_reverse() + # Standard matrix view (row 1 at top)
  coord_fixed() +
  theme_minimal() +
  labs(title = "Circulant Symmetric Matrix Sigma",
       subtitle = "Highlighting bands at 0.1 and diagonal at 1",
       fill = "Value")

################################
# Visualize SSLASSO Simulation #
################################

plot(fit_sslasso$coef,ylim=c(-1.5*b,1.5*b), main='Confidence Intervals based on de-biased LASSO', ylab='', xlab = 'Coefficients');
points(fit_sslasso$unb.coef,col="blue");
points(theta0,col="green");
lines(fit_sslasso$up.lim,col="red");
lines(fit_sslasso$low.lim,col="red");
legend('topright', legend=c('LASSO','de-biased LASSO','Ground-truth','Confidence Intervals'), col=c('black', 'blue','green','red'), pch=c(1,1,1,NA_integer_), lty = c(0,0,0,1))


print(paste("Estimated Noise Sdev = ",fit_sslasso$noise.sd));
print(paste("Estimated Norm0 = ",fit_sslasso$norm0));
nc <- sum(as.numeric(fit_sslasso$low.lim>theta0)+as.numeric(fit_sslasso$up.lim<theta0));
print(paste("Coverage = ",(p-nc)/p));
print("P-values:")
print(fit_sslasso$pvals)
