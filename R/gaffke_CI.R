#' Gaffke bound as described in Learned-Miller and Thomas https://arxiv.org/pdf/1905.06208.pdf
#' @export
#' @param x A vector of iid random samples from some distribution F with mean mu
#' @param alpha The level; the output is a (1-alpha) confidence bound
#' @param B the number of Monte Carlo iterations to estimate the interval
#' @param side the side of the confidence interval, either "upper" or "lower"
#' @return An upper or lower level alpha confidence bound on the population mean mu
gaffke_CI <- function(x, alpha = .05, B = 10000, side = "upper"){
  n <- length(x)
  if(side == "lower"){
    x <- 1 - x
  }
  z <- sort(x, decreasing = FALSE)
  ms <- rep(NA, B)
  s <- c(diff(z), 1 - z[n])
  for(i in 1:B){
    u <- sort(runif(n), decreasing = FALSE)
    ms[i] <- 1 - sum(u * s)
  }
  ms_alpha <- quantile(ms, 1 - alpha)
  if(side == "lower"){
    1 - ms_alpha
  } else if(side == "upper"){
    ms_alpha
  }
}