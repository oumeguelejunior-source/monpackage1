#' Print Method for Laplace Measurement Error Objects
#'
#' Displays a formatted summary of the results obtained from a
#' Laplace measurement error model fitted with
#' \code{\link{estimate_laplace_errors}}.
#'
#' The output includes the number of observations, estimated parameters,
#' variance, standard deviation, confidence intervals, and the value of the
#' log-likelihood.
#'
#' @param x An object of class \code{"laplace_measure"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} is returned invisibly.
#'
#' @export
#'
#' @method print laplace_measure


print.laplace_measure <- function(x, ...) {
  cat("Modèle d'erreur de mesure - Loi de Laplace\n")
  cat("==========================================\n")
  cat(sprintf("Nombre d'observations: %d\n", x$data$n))
  cat("\nParamètres estimés:\n")
  cat(sprintf("  μ (erreur systématique): %.4f\n", x$parameters$mu))
  cat(sprintf("  b (paramètre d'échelle): %.4f\n", x$parameters$b))
  cat(sprintf("  Variance: %.4f\n", x$variance))
  cat(sprintf("  Écart-type: %.4f\n", sqrt(x$variance)))
  cat("\nIntervalles de confiance à %.1f%%:\n", (1-x$alpha)*100)
  cat(sprintf("  μ: [%.4f, %.4f]\n", x$confidence_intervals$mu[1], x$confidence_intervals$mu[2]))
  cat(sprintf("  b: [%.4f, %.4f]\n", x$confidence_intervals$b[1], x$confidence_intervals$b[2]))
  cat(sprintf("\nLog-vraisemblance: %.4f\n", x$log_likelihood))
}
#' Print Method for Laplace Measurement Error Objects
#'
#' Displays a formatted summary of the results obtained from a
#' Laplace measurement error model fitted with
#' \code{\link{estimate_laplace_errors}}.
#'
#' The output includes the number of observations, estimated parameters,
#' variance, standard deviation, confidence intervals, and the value of the
#' log-likelihood.
#'
#' @param x An object of class \code{"laplace_measure"}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The object \code{x} is returned invisibly.
#'
#' @export
#'
#' @method print laplace_measure


summary.laplace_measure <- function(object, ...) {
  list(
    parameters = object$parameters,
    variance = object$variance,
    confidence_intervals = object$confidence_intervals,
    n = object$data$n,
    log_likelihood = object$log_likelihood,
    AIC = -2 * object$log_likelihood + 4,
    BIC = -2 * object$log_likelihood + 2 * log(object$data$n)
  )
}

#' Diagnostic Plots for Laplace Measurement Error Models
#'
#' Produces diagnostic plots for a fitted Laplace measurement error model.
#'
#' The following plots are displayed:
#' \enumerate{
#'   \item Histogram of the errors with the fitted Laplace density.
#'   \item Normal Q--Q plot of residuals against theoretical Laplace quantiles.
#'   \item Errors plotted against observation index.
#'   \item Empirical cumulative distribution function (ECDF) with fitted
#'   Laplace cumulative distribution function.
#' }
#'
#' These plots are intended to assess the adequacy of the Laplace model
#' and to identify potential deviations such as outliers or model misspecification.
#'
#' @param x An object of class \code{"laplace_measure"}.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
#'
#' @method plot laplace_measure


plot.laplace_measure <- function(x, ...) {
  old_par <- par(mfrow = c(2, 2))
  on.exit(par(old_par))

  errors <- x$data$errors
  mu <- x$parameters$mu
  b <- x$parameters$b

  # Histogramme avec densité théorique
  hist(errors, prob = TRUE, main = "Distribution des erreurs",
       xlab = "Erreurs", ylab = "Densité", col = "lightblue")
  curve(1/(2*b) * exp(-abs(x - mu)/b), add = TRUE, col = "red", lwd = 2)
  legend("topright", legend = c("Données", "Loi de Laplace"),
         col = c("lightblue", "red"), lwd = c(NA, 2), fill = c("lightblue", NA))

  # QQ-plot
  theoretical_quantiles <- mu + b * sign(runif(length(errors)) - 0.5) *
    log(1 - 2*abs(runif(length(errors)) - 0.5))
  qqnorm(errors - theoretical_quantiles,
         main = "QQ-plot des résidus")
  qqline(errors - theoretical_quantiles, col = "red")

  # Tracé des erreurs dans l'ordre
  plot(errors, type = "b", main = "Erreurs par observation",
       xlab = "Index", ylab = "Erreur", pch = 19)
  abline(h = mu, col = "red", lwd = 2)
  abline(h = mu + c(-2, 2)*sqrt(x$variance), col = "red", lty = 2)

  # Fonction de densité cumulative
  plot(ecdf(errors), main = "Fonction de répartition",
       xlab = "Erreurs", ylab = "Probabilité cumulative")
  curve(0.5 * (1 + sign(x - mu) * (1 - exp(-abs(x - mu)/b))),
        add = TRUE, col = "red", lwd = 2)
}
