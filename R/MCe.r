#' @title  compute Monte Carlo estimate of the Beta(3,3) cdf.
#' @description Use a function to compute Monte Carlo estimate of the Beta(3,3) cdf.
#' @param x vector of quantiles
#' @param n number of random sampling
#' @return a vector of size x\code{n}
#' @examples
#' \dontrun{
#' x <- seq(0.1, 0.9, 0.1)
#' n <- 10000
#' res <- MCe(x, n)
#' }
#' @export
MCe <- function(x, n)
{
  y <- runif(n)
  cdf <- numeric(length(x))

  for (i in 1:length(x))
  {
      g <- (x[i])^3*(y*(1-y*x[i]))^2
      cdf[i] <- mean(g)/beta(3,3)
  }
   return(cdf)
}
