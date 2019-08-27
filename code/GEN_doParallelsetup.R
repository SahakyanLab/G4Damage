library(doParallel)
if(NCPU > 1){
  registerDoParallel(cores = NCPU)
  `%op%` <- `%dopar%`
  print(paste("Running with ", NCPU, " core(s).", sep=""), quote=F)
} else {
              `%op%` <- `%do%`
}
