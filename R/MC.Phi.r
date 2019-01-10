#' @title generate samples from a Rayleigh distribution.
#' @description Implement a function to generate samples from a Rayleigh distribution.
#' @param x vector of quantiles
#' @param sigma  a parameter
#' @param r number of sample
#' @param antithetic logical value, the default is true
#' @return a vector of size x\code{n}
#' @examples
#' \dontrun{
#' x <- seq(1, 10, 1)
#' sigma <- 2
#' r <- 10000
#' res <- MC.Phi(x, sigma, r)
#' }
#' @export
MC.Phi <- function(x, sigma, r, antithetic = TRUE)
{
  y <- runif(r)
  if(!antithetic) v <- runif(r)
  else
    v <- 1-y
  y <- c(y, v)
  cdf <- numeric(length(x))
  for(i in 1:length(x))
  {
    g <- (y*x[i]^2)*exp(-(y*x[i]^2)/(2*(sigma^2)))
    cdf[i] <-  1/(sigma^2)*mean(g)
  }
  cdf
}
