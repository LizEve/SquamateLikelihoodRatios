library(ggplot2)
mainDir <- "/Users/ChatNoir/Projects/Squam/BayesFactors/"
setwd(mainDir)

bfs <- read.csv("BFs.out", stringsAsFactors = FALSE, header=T)
hist(bfs$lnBF, breaks=20)
abline(v=median(simulated_inference_pps_example$mean_rf), col="green", lty=2)
abline(v=empirical_inference_pps_example$mean_rf, col="red", lty=2)


mainDir <- "/Users/ChatNoir/Projects/Squam/BayesFactors/threegGeneTrees/uce-314/results_GTRG/"
setwd(mainDir)
empirical_inference_pps_example <- read.csv("empirical_inference_uce-314.csv", header=TRUE)
simulated_inference_pps_example <- read.csv("simulated_inference_uce-314.csv", header=TRUE)
hist(simulated_inference_pps_example$mean_rf, breaks=20)
abline(v=median(simulated_inference_pps_example$mean_rf), col="green", lty=2)
abline(v=empirical_inference_pps_example$mean_rf, col="red", lty=2)
