# LaplaceMeasure 📊

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/oumeguelejunior-source/LaplaceMeasure/workflows/R-CMD-check/badge.svg)](https://github.com/oumeguelejunior-source/LaplaceMeasure/actions)
[![Codecov test coverage](https://codecov.io/gh/oumeguelejunior-source/LaplaceMeasure/branch/main/graph/badge.svg)](https://app.codecov.io/gh/oumeguelejunior-source/LaplaceMeasure)

> **Robust Measurement Error Modeling with Laplace Distribution** for instrument calibration, uncertainty quantification, and precision analysis.

## ✨ Features

- **Maximum Likelihood Estimation** for Laplace distribution parameters
- **Robust uncertainty quantification** using likelihood ratio methods
- **Comprehensive diagnostic plots** for model validation
- **Comparison tools** for instrument evaluation
- **Production-ready** with full S3 method support

## 📦 Installation

```r
# Install from GitHub
devtools::install_github("oumeguelejunior-source/laplace1")

# Or install from local source
devtools::install_local("path/to/LaplaceMeasure")
```

## 🚀 Quick Start

```r
library(LaplaceMeasure)

# Analyze measurement errors from an instrument
measurements <- c(10.1, 10.2, 9.8, 10.0, 10.3)
true_value <- 10.0

result <- estimate_laplace_errors(
  measurements = measurements,
  true_value = true_value,
  alpha = 0.05
)

# View results
print(result)
plot(result)
```

## 📊 Example: Instrument Comparison

```r
# Simulate two temperature sensors
set.seed(123)
n <- 30
true_temp <- 37.0

# Sensor A: Good precision, slight bias
sensor_A <- true_temp + rlaplace(n, mu = 0.08, b = 0.12)
# Sensor B: No bias, less precise
sensor_B <- true_temp + rlaplace(n, mu = 0.02, b = 0.18)

# Analyze both sensors
result_A <- estimate_laplace_errors(sensor_A, true_temp)
result_B <- estimate_laplace_errors(sensor_B, true_temp)

# Compare performance
compare_sensors(result_A, result_B)

# Generate comprehensive report
generate_report(result_A, result_B, true_temp)
```

## 🔧 Core Functions

### Parameter Estimation
```r
estimate_laplace_errors()  # Main estimation function
dlaplace()                 # Laplace density function
plaplace()                 # Laplace distribution function
rlaplace()                 # Random generation from Laplace
```

### Analysis & Diagnostics
```r
print()                    # Concise results summary
summary()                  # Detailed statistics (AIC, BIC)
plot()                     # Diagnostic plots (4 plots)
predict_interval()         # Prediction intervals
```

### Instrument Evaluation
```r
compare_sensors()          # Statistical comparison of instruments
analyze_robustness()       # Outlier impact analysis
generate_report()          # Comprehensive evaluation report
```

## 📈 Visual Diagnostics

The package provides comprehensive diagnostic plots:

```r
plot(result)  # Generates 4-panel diagnostic plot
```

1. **Histogram with fitted density** - Check distribution fit
2. **Q-Q plot** - Assess distributional assumptions
3. **Time series of errors** - Detect trends or outliers
4. **Empirical vs theoretical CDF** - Goodness-of-fit test

## 🎯 Real-World Applications

### 1. **Instrument Calibration**
```r
# Determine systematic bias for recalibration
calibration_data <- read.csv("calibration_data.csv")
cal_result <- estimate_laplace_errors(
  measurements = calibration_data$readings,
  true_value = calibration_data$reference
)
cat(sprintf("Recalibration offset: %.3f units\n", cal_result$parameters$mu))
```

### 2. **Quality Control**
```r
# Monitor measurement precision over time
qc_results <- lapply(daily_measurements, function(day_data) {
  estimate_laplace_errors(day_data, reference_value)
})
precision_trend <- sapply(qc_results, function(r) r$parameters$b)
plot(precision_trend, type = "b", main = "Precision Trend Over Time")
```

### 3. **Method Validation**
```r
# Compare new vs reference measurement method
new_method <- c(25.1, 25.3, 24.9, 25.2, 25.0)
reference_method <- c(25.0, 25.0, 25.0, 25.0, 25.0)

errors <- new_method - reference_method
validation_result <- estimate_laplace_errors(errors = errors)

if(validation_result$confidence_intervals$mu[2] < 0.1) {
  cat("✓ Method bias within acceptable limits\n")
}
```

## 📋 Statistical Model

The package implements the Laplace (double exponential) distribution for modeling measurement errors:

\[
f(x|\mu,b) = \frac{1}{2b} \exp\left(-\frac{|x-\mu|}{b}\right)
\]

Where:
- **μ** (mu): Systematic error (bias) - estimated by median
- **b**: Scale parameter (precision) - estimated by mean absolute deviation
- **Variance**: \(2b^2\)

### Parameter Estimation
- **MLE for μ**: \(\hat{\mu} = \text{median}(e_i)\)
- **MLE for b**: \(\hat{b} = \frac{1}{n}\sum_{i=1}^n |e_i - \hat{\mu}|\)

### Confidence Intervals
- **μ**: Based on asymptotic normality
- **b**: Using likelihood ratio method for better accuracy

## 📚 Advanced Usage

### Custom Confidence Levels
```r
# 99% confidence intervals
result_99 <- estimate_laplace_errors(measurements, true_value, alpha = 0.01)

# Prediction intervals for new measurements
pred_95 <- predict_interval(result, level = 0.95)
pred_99 <- predict_interval(result, level = 0.99)
```

### Batch Processing
```r
# Analyze multiple instruments
instruments <- list("Sensor_A" = sensor_A, "Sensor_B" = sensor_B, "Sensor_C" = sensor_C)
results <- lapply(instruments, function(x) estimate_laplace_errors(x, true_value))

# Create comparison table
comparison_table <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    Bias = r$parameters$mu,
    Precision = r$parameters$b,
    Variance = r$variance,
    AIC = summary(r)$AIC
  )
}))
```

## 🧪 Testing

```r
# Run test suite
library(testthat)
test_package("LaplaceMeasure")

# Example test
test_that("Parameter estimation works", {
  errors <- rlaplace(100, mu = 0.1, b = 0.5)
  result <- estimate_laplace_errors(errors = errors)
  expect_equal(result$parameters$mu, median(errors), tolerance = 0.1)
})
```

## 📖 Documentation

Full documentation is available:

```r
# View function documentation
?estimate_laplace_errors
?plot.laplace_measure
?compare_sensors

# Browse vignettes
browseVignettes("LaplaceMeasure")
```

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## 📄 License

This package is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## 🔗 Related Packages

- `fitdistrplus` - General distribution fitting
- `MASS` - Robust statistical methods
- `metRology` - Measurement uncertainty

## 📧 Support

For questions, issues, or feature requests:
- 📝 [Open an issue](https://github.com/oumeguelejunior-source/LaplaceMeasure/issues)
- 💬 Start a discussion on GitHub

---

**LaplaceMeasure** - Making measurement uncertainty analysis robust, accessible, and statistically sound. 📐🔬

*"Precision in measurement is the foundation of scientific truth."*
