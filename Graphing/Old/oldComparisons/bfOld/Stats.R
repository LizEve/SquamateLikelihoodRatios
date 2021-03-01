library(ggplot2)


histogram <- function(sim,emp,xlab,main,breaks=25){
  low=min(min(sim),emp)
  high=max(max(sim),emp)
  hist(sim, breaks=25, xlim=c(low,high), xlab = xlab, main = main)
  abline(v=median(sim), col="orange", lty=1, lwd=3)
  abline(v=emp, col="purple", lty=3, lwd=3)
}

setwd(getwd())
getwd()
empIn <- "uce-314/results_GTRG/empirical_inference_uce-314.csv"
simIn <- "uce-314/results_GTRG/simulated_inference_uce-314.csv"
#empIn <- "uce-12/results_GTRG/empirical_inference_uce-12.csv"
#simIn <- "uce-12/results_GTRG/simulated_inference_uce-12.csv"
#empIn <- "AHE-L13/results_GTRG/empirical_inference_AHE-L13.csv"
#simIn <- "AHE-L13/results_GTRG/simulated_inference_AHE-L13.csv"

empirical_inference_pps_example <- read.csv(empIn, header=TRUE)
simulated_inference_pps_example <- read.csv(simIn, header=TRUE)
jpeg("uce-314.jpg",width = 600, height = 600)
#jpeg("uce-12.jpg",width = 600, height = 600)
#jpeg("AHE-L13.jpg",width = 600, height = 600)
par(mfrow=c(3,3))
# RF mean
histogram(simulated_inference_pps_example$mean_rf,empirical_inference_pps_example$mean_rf, xlab="Mean RF", main = "Mean RF")
# Tree length mean
histogram(simulated_inference_pps_example$mean_tl,empirical_inference_pps_example$mean_tl, xlab="Mean Tree Length", main = "Mean Tree Length")
# Tree length varience 
histogram(simulated_inference_pps_example$var_tl,empirical_inference_pps_example$var_tl, xlab="Variance Tree Length", main = "Variance Tree Length")
# Entropy
histogram(simulated_inference_pps_example$entropy,empirical_inference_pps_example$entropy,xlab="Entropy",main = "Entropy")
# Quantile 25
histogram(simulated_inference_pps_example$quantile25,empirical_inference_pps_example$quantile25, xlab="25th Quantile", main = "25th Quantile")
# Quantile 50
histogram(simulated_inference_pps_example$quantile50,empirical_inference_pps_example$quantile50, xlab="50th Quantile", main = "50th Quantile")
# Quantile 75
histogram(simulated_inference_pps_example$quantile75,empirical_inference_pps_example$quantile75, xlab="75th Quantile", main = "75th Quantile")
# Quantile 99
histogram(simulated_inference_pps_example$quantile99,empirical_inference_pps_example$quantile99, xlab="99th Quantile", main = "99th Quantile")
# Quantile 999
histogram(simulated_inference_pps_example$quantile999,empirical_inference_pps_example$quantile999, xlab="999th Quantile", main = "999th Quantile")
dev.off()

