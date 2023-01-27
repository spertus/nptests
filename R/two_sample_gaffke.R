#' Two-sample test of mean dominance
#' @export
#' @description A function to run a one-sided, two-sample nonparametric test based on the Gaffke test.
#' In particular, given IID samples from two distributions with means \eqn{\mu_1} and \eqn{\mu_2} the hypothesis H_0: \eqn{\mu_2 <= \mu_1} is tested. 
#' The test assumes that the distributions share common bounds \[a,b\], where a is the smallest and b is the largest value that samples from the populations can take.
#' If the distribution is contained within these bounds, the test is valid (P-values are dominated by the uniform distribution) no matter the shape of the population distributions.
#' Standard two-sample tests (e.g., the t-test) are not valid over all bounded distributions in finite-samples, and may fail to control the significance level.
#' @param sample_1 a numeric vector, a random sample from a population with unknown mean \eqn{\mu_1}
#' @param sample_2 a numeric vector, a random sample from a second population with unknown mean \eqn{\mu_2}
#' @param bounds a length-2 vector of the (shared) upper and lower bounds on both populations, \[a,b\]
#' @param B a positive integer, the number of Monte Carlo iterations to be used in running the test
#' @param method a string, the method of combining two one-sample tests into a two-sample test of equality
#' * Sidak: using Sidak's correction, compute a 1-sqrt(1-alpha) upper and lower bound (respectively) and see if they overlap
#' * Fisher: compute a combined p-value using Fisher's combining function of whether both means are equal to a particular mu_0, then maximize this over possible values of mu_0
#' * Liptak: compute a combined p-value using Liptak's combining function
#' @return A P-value indicating for the hypothesis test H_0: \eqn{\mu_2 <= \mu_1}
#' @examples 
#' two_sample_gaffke(sample_1 = rbeta(30, shape1 = 5, shape2 = 0.5), sample_2 = rbeta(30, shape1 = 100, shape2 = 10))
two_sample_gaffke <- function(sample_1, sample_2, B = 1000, method = "fisher", bounds = c(0,1)){
  if(anyNA(sample_1) | anyNA(sample_2)){stop("Samples have NA values")}
  
  sample_1 <- (sample_1 - bounds[1]) / diff(bounds)
  sample_2 <- (sample_2 - bounds[1]) / diff(bounds)
  if(method == "sidak"){
    upper_1 <- gaffke_CI(x = sample_1, alpha = 1-sqrt(1-alpha), B = B, side = "upper")
    lower_2 <- gaffke_CI(x = sample_2, alpha = 1-sqrt(1-alpha), B = B, side = "lower")
    reject <- ifelse(upper_1 < lower_2, TRUE, FALSE)
    reject
  } else if(method %in% c("fisher","liptak","tippett")){
    n_1 <- length(sample_1)
    n_2 <- length(sample_2)
    x_1 <- sample_1
    x_2 <- sample_2
    ms_1 <- rep(NA, B)
    ms_2 <- rep(NA, B)
    for(b in 1:B){
      Z_1 <- rexp(n_1 + 1)
      D_1 <- Z_1 / sum(Z_1)
      Z_2 <- rexp(n_2 + 1)
      D_2 <- Z_2 / sum(Z_2)
      ms_1[b] <- sum(D_1 * c(x_1, 1))
      ms_2[b] <- sum(D_2 * c(x_2, 0))
    }
    if(method == "fisher"){
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- -2 * (log(p_1) + log(p_2))
        p_val <- pchisq(q = combined_test_stat, df = 4, lower.tail = FALSE)
        p_val
      }
    } else if(method == "liptak"){
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- qnorm(1 - p_1) + qnorm(1 - p_2)
        p_val <- pnorm(q = combined_test_stat, sd = 2, lower.tail = FALSE)
        p_val
      }
    } else{
      combined_p <- function(mu_0){
        p_1 <- mean(c(ms_1 >= mu_0, TRUE))
        p_2 <- mean(c(ms_2 <= mu_0, TRUE))
        combined_test_stat <- min(p_1, p_2)
        p_val <- pbeta(q = combined_test_stat, shape1 = 1, shape2 = 2, lower.tail = TRUE)
        p_val
      }
    }
    
    max_p_val <- optimize(combined_p, lower = min(c(ms_1, ms_2)), upper = max(c(ms_1, ms_2)), maximum = TRUE)$objective
    max_p_val
  } else{
    stop("Input a valid method: sidak, fisher, liptak, or tippett")
  }
}
