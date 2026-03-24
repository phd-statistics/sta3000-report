# Simulation and Data Analysis Codes for "Confidence Intervals and Hypothesis Testing for High-Dimensional Regression: A Review and Discussion"

This repository contains the replication codes for the report "Confidence Intervals and Hypothesis Testing for High-Dimensional Regression: A Review and Discussion" by Benjamin Smith. The code aims to replicate the the results of Javanmard and Montari's JMLR publication ["Confidence Intervals and Hypothesis Testing For High-Dimensional Regression" (2014)](https://web.stanford.edu/~montanar/sslasso/) based on their description, but relaxes the specifications related to specification of the regularization parameter and optimization constraint. 

To make the replication experience seamless, please ensure that you have the install the following packages with the following line:

```
# Use the pak package to ensure installation is seemless
# install.packages("pak")
pak::pkg_install(c("Matrix", "glmnet", "expm", "flare","hdi", "mvnfast", "ggplot2"))
```
# Simulation Study

To run the simulation study, ensure that you