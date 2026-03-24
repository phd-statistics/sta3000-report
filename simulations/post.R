library(xtable)

res_main <- readRDS("./res_commute_backup.rds")
res_ms <- readRDS("./res_commute_backup.rds") 


tab1 <- do.call(rbind, lapply(res_main[1:6], function(x) {
  # Extract params (handling them as vectors)
  p_str <- paste0("(", x$params[["n"]], ",", x$params[["p"]], ",", 
                  x$params[["s0"]], ",", x$params[["b"]], ")")
  
  # Extract metrics using [ ] since sslasso is a named vector
  res <- x$metrics$sslasso
  
  data.frame(
    configuration = p_str,
    ell     = res["AvgLength"],
    ell_S   = res["AvgLength_S"],
    ell_Sc  = res["AvgLength_Sc"],
    dCov    = res["dCov"],
    dCov_S  = res["dCov_S"],
    dCov_Sc = res["dCov_Sc"],
    row.names = NULL # Prevents messy row names in the final table
  )
}))