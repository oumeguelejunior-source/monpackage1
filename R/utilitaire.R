#' Density Function of the Laplace Distribution
#'
#' Evaluates the probability density function (PDF) of the Laplace
#' (double exponential) distribution at specified values.
#'
#' @param x Numeric vector of quantiles.
#' @param mu Location parameter of the Laplace distribution.
#' @param b Scale parameter of the Laplace distribution (must be positive).
#'
#' @return A numeric vector giving the value of the density function at \code{x}.
#'
#' @details
#' The Laplace density is defined as
#' \deqn{
#' f(x \mid \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right).
#' }
#'
#' @export

dlaplace <- function(x, mu, b) {
  if (b <= 0) stop("b doit être positif")
  1/(2*b) * exp(-abs(x - mu)/b)
}

#' Distribution Function of the Laplace Distribution
#'
#' Computes the cumulative distribution function (CDF) of the Laplace
#' distribution.
#'
#' @param q Numeric vector of quantiles.
#' @param mu Location parameter of the Laplace distribution.
#' @param b Scale parameter of the Laplace distribution (must be positive).
#'
#' @return A numeric vector giving the cumulative probabilities
#' \eqn{P(X \le q)}.
#'
#' @details
#' The cumulative distribution function is given by
#' \deqn{
#' F(q \mid \mu, b) =
#' \begin{cases}
#' \frac{1}{2} \exp\left(\frac{q - \mu}{b}\right), & q < \mu, \\
#' 1 - \frac{1}{2} \exp\left(-\frac{q - \mu}{b}\right), & q \ge \mu.
#' \end{cases}
#' }
#'
#' @export


plaplace <- function(q, mu, b) {
  if (b <= 0) stop("b doit être positif")
  0.5 * (1 + sign(q - mu) * (1 - exp(-abs(q - mu)/b)))
}

#' Random Generation from the Laplace Distribution
#'
#' Generates random deviates from a Laplace (double exponential)
#' distribution using the inverse transform method.
#'
#' @param n Number of observations to generate.
#' @param mu Location parameter of the Laplace distribution.
#' @param b Scale parameter of the Laplace distribution (must be positive).
#'
#' @return A numeric vector of length \code{n} containing Laplace-distributed
#' random values.
#'
#' @details
#' Random values are generated using the inverse cumulative distribution
#' function applied to uniform random variables.
#'
#' @export


rlaplace <- function(n, mu, b) {
  if (b <= 0) stop("b doit être positif")
  u <- runif(n) - 0.5
  mu - b * sign(u) * log(1 - 2*abs(u))
}

#' Predictive Interval for a New Measurement
#'
#' Computes a predictive confidence interval for a future measurement
#' based on a fitted Laplace measurement error model.
#'
#' @param object An object of class \code{"laplace_measure"}.
#' @param level Confidence level for the predictive interval
#'   (default is 0.95).
#'
#' @return A list containing:
#' \itemize{
#'   \item \strong{level}: Confidence level of the interval.
#'   \item \strong{interval}: Lower and upper bounds of the predictive interval.
#'   \item \strong{mu}: Estimated location parameter.
#'   \item \strong{b}: Estimated scale parameter.
#' }
#'
#' @details
#' The predictive interval is derived from the quantiles of the Laplace
#' distribution:
#' \deqn{
#' P\left( \mu + b \log(\alpha) \le X \le
#' \mu - b \log(\alpha) \right) = 1 - \alpha,
#' }
#' where \eqn{\alpha = 1 - \text{level}}.
#'
#' @export


predict_interval <- function(object, level = 0.95) {
  alpha <- 1 - level
  mu <- object$parameters$mu
  b <- object$parameters$b

  lower <- mu + b * log(2 * alpha/2)
  upper <- mu - b * log(2 * alpha/2)

  list(
    level = level,
    interval = c(lower, upper),
    mu = mu,
    b = b
  )
}
