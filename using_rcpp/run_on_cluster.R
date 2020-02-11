t1 <- Sys.time()
source("functions.R")
library(parallel)

args <- commandArgs(TRUE)

beta <- as.numeric(args[1])
pf1 <- as.numeric(args[2])
delta <- as.numeric(args[3])
r <- as.numeric(args[4])

# unchanged parameters
N <- 300
eta <- 1

dat_list <- lapply(1:N, function(x)
  {generate_censored_data(r, pf1, beta, eta)})
output <- mclapply(dat_list, function(x) {
  prediction_five_methods(x, qweibull(pf1+delta, beta, eta), beta, eta)
}, mc.cores = 3)

saveRDS(output, file = paste("b",beta,"er", r,"p",pf1,"d",delta,".rds",sep=""))
Sys.time() -t1