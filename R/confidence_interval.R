#' Likelihood Ratio Confidence Interval for the Scale Parameter b
#'
#' Computes a confidence interval for the scale parameter \eqn{b} of a
#' Laplace distribution using the likelihood ratio method.
#'
#' @param errors Numeric vector of observed errors.
#' @param mu_hat Estimated location parameter \eqn{\mu}.
#' @param b_hat Estimated scale parameter \eqn{b}.
#' @param alpha Significance level (e.g., 0.05 for a 95\% confidence interval).
#'
#' @return A numeric vector of length two containing the lower and upper
#' bounds of the confidence interval for \eqn{b}.
#'
#' @details
#' The confidence interval is obtained by solving
#' \deqn{2(\ell(\hat{b}) - \ell(b)) = \chi^2_{1, 1-\alpha},}
#' where \eqn{\ell(b)} is the profile log-likelihood.
#'
#' @keywords internal
#'
estimate_b_ci <- function(errors, mu_hat, b_hat, alpha) {
  n <- length(errors)

  # Profile log-likelihood for b
  profile_loglik <- function(b) {
    -n * log(2 * b) - sum(abs(errors - mu_hat)) / b
  }

  # Likelihood ratio statistic
  D <- function(b) {
    2 * (profile_loglik(b_hat) - profile_loglik(b))
  }

  # Root-finding functions
  f_lower <- function(b) D(b) - qchisq(1 - alpha, 1)
  f_upper <- function(b) D(b) - qchisq(1 - alpha, 1)

  # Numerical search for bounds
  lower_bound <- tryCatch({
    uniroot(f_lower, c(b_hat * 0.01, b_hat))$root
  }, error = function(e) NA)

  upper_bound <- tryCatch({
    uniroot(f_upper, c(b_hat, b_hat * 10))$root
  }, error = function(e) NA)

  c(lower_bound, upper_bound)
}
