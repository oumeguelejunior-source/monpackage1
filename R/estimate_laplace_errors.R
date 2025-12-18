#' Estimation of Measurement Error Parameters under a Laplace Model
#'
#' Estimates the parameters of a Laplace-distributed measurement error model.
#' The function supports two usage modes: (i) computing errors from observed
#' measurements and a known true value, or (ii) directly providing precomputed
#' errors.
#'
#' Under the Laplace assumption, the location parameter \eqn{\mu} is estimated
#' by the sample median, while the scale parameter \eqn{b} is estimated by the
#' mean absolute deviation from the median.
#'
#' @param measurements Numeric vector of observed measurements.
#' @param true_value Known true value of the measured quantity (optional).
#' @param errors Numeric vector of precomputed measurement errors
#'   (alternative to \code{measurements} and \code{true_value}).
#' @param alpha Significance level for confidence intervals
#'   (default is 0.05, corresponding to a 95\% confidence level).
#'
#' @return An object of class \code{"laplace_measure"} containing:
#' \itemize{
#'   \item \strong{parameters}: A list with estimated parameters \eqn{\mu} and \eqn{b}.
#'   \item \strong{variance}: Estimated variance of the Laplace distribution
#'     (\eqn{2b^2}).
#'   \item \strong{confidence_intervals}: Confidence intervals for \eqn{\mu} and \eqn{b}.
#'   \item \strong{data}: List containing the data used in the estimation.
#'   \item \strong{log_likelihood}: Value of the log-likelihood at the estimated parameters.
#'   \item \strong{alpha}: Significance level used for inference.
#' }
#'
#' @details
#' Let the measurement error \eqn{X} follow a Laplace distribution with density
#' \deqn{
#' f(x \mid \mu, b) = \frac{1}{2b} \exp\left(-\frac{|x - \mu|}{b}\right).
#' }
#'
#' The maximum likelihood estimator of \eqn{\mu} is the sample median, while
#' the estimator of \eqn{b} is
#' \deqn{
#' \hat{b} = \frac{1}{n} \sum_{i=1}^n |X_i - \hat{\mu}|.
#' }
#'
#' The confidence interval for \eqn{\mu} is based on a normal approximation,
#' while the confidence interval for \eqn{b} is obtained using the likelihood
#' ratio method.
#'
#' @export
#'
#' @examples
#' # Example 1: Known true value
#' measurements <- c(10.1, 10.2, 9.8, 10.0, 10.3)
#' result <- estimate_laplace_errors(
#'   measurements = measurements,
#'   true_value = 10
#' )
#'
#' # Example 2: Precomputed errors
#' errors <- c(0.1, 0.2, -0.2, 0.0, 0.3)
#' result <- estimate_laplace_errors(errors = errors)
#'
estimate_laplace_errors <- function(measurements = NULL,
                                    true_value = NULL,
                                    errors = NULL,
                                    alpha = 0.05) {

  # Input validation
  if (is.null(errors)) {
    if (is.null(measurements) || is.null(true_value)) {
      stop("Either 'errors' or both 'measurements' and 'true_value' must be provided")
    }
    errors <- measurements - true_value
  }

  if (!is.numeric(errors) || length(errors) < 2) {
    stop("errors must be a numeric vector of length at least 2")
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1)")
  }

  n <- length(errors)

  # Parameter estimation
  mu_hat <- median(errors)
  b_hat  <- mean(abs(errors - mu_hat))

  # Variance
  variance <- 2 * b_hat^2

  # Log-likelihood
  log_lik <- -n * log(2 * b_hat) - sum(abs(errors - mu_hat)) / b_hat

  # Confidence interval for mu (normal approximation)
  se_mu <- sqrt(variance / n)
  ci_mu <- mu_hat + c(-1, 1) * qnorm(1 - alpha / 2) * se_mu

  # Confidence interval for b (likelihood ratio)
  ci_b <- estimate_b_ci(errors, mu_hat, b_hat, alpha)

  # Result object
  result <- list(
    parameters = list(mu = mu_hat, b = b_hat),
    variance = variance,
    confidence_intervals = list(mu = ci_mu, b = ci_b),
    data = list(
      errors = errors,
      n = n,
      measurements = measurements,
      true_value = true_value
    ),
    log_likelihood = log_lik,
    alpha = alpha
  )

  class(result) <- "laplace_measure"
  result
}
