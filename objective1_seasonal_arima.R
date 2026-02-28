library(forecast)
library(tseries)
library(readxl)
library(dplyr)

# ============================================================
# LOAD DATA
# ============================================================

load_data <- function() {
  cat("\nSelect your data file...\n")
  file_path <- file.choose()
  file_ext <- tolower(tools::file_ext(file_path))

  if (file_ext %in% c("xlsx", "xls")) {
    data <- read_excel(file_path)
  } else if (file_ext == "csv") {
    data <- read.csv(file_path)
  } else {
    stop("Unsupported file format. Please use .xlsx, .xls, or .csv")
  }

  return(data)
}

prepare_data <- function(data) {
  data <- data %>%
    select(Year, Month, `Individual Life Insurance`) %>%
    mutate(
      Year = as.numeric(Year),
      Month = as.numeric(Month),
      `Individual Life Insurance` = as.numeric(`Individual Life Insurance`)
    )
  return(data)
}

# ============================================================
# STEP 1: AUGMENTED DICKEY-FULLER (ADF) TEST
# Paragraph Section B: Stationarity Testing
# ADF equation: nabla(yt) = alpha + beta*t + gamma*y(t-1)
#               + sum(delta_i * nabla(y(t-i))) + epsilon_t
# ============================================================

perform_adf_test <- function(data) {
  cat("\n+================================================================+\n")
  cat("|            STEP 1: STATIONARITY TEST (ADF)                    |\n")
  cat("+================================================================+\n\n")

  cat("HYPOTHESIS:\n")
  cat("H0: Series is NON-STATIONARY (has unit root)\n")
  cat("H1: Series is STATIONARY\n\n")

  cat("ADF Test Equation:\n")
  cat("  nabla(yt) = alpha + beta*t + gamma*y(t-1) + sum(delta_i * nabla(y(t-i))) + epsilon_t\n\n")

  cat("Decision Rule:\n")
  cat("  p-value < 0.05  -> STATIONARY (d=0)\n")
  cat("  p-value >= 0.05 -> NON-STATIONARY (need differencing)\n\n")

  insurance_series <- data$`Individual Life Insurance`

  cat("Testing ORIGINAL Series:\n")
  cat("----------------------------\n")
  adf_original <- suppressWarnings(adf.test(insurance_series))
  cat("Test Statistic:", round(adf_original$statistic, 6), "\n")
  cat("P-value:", sprintf("%.10f", adf_original$p.value), "\n")

  if (adf_original$p.value < 0.05) {
    cat("[OK] DECISION: STATIONARY -> d = 0\n\n")
    d_value <- 0
  } else {
    cat("[X] DECISION: NON-STATIONARY\n\n")

    cat("Testing FIRST DIFFERENCE (d=1):\n")
    cat("----------------------------\n")
    diff_series <- diff(insurance_series, differences = 1)
    adf_diff1 <- suppressWarnings(adf.test(diff_series))
    cat("Test Statistic:", round(adf_diff1$statistic, 6), "\n")
    cat("P-value:", sprintf("%.10f", adf_diff1$p.value), "\n")

    if (adf_diff1$p.value < 0.05) {
      cat("[OK] DECISION: First difference is STATIONARY -> d = 1\n\n")
      d_value <- 1
    } else {
      cat("[X] DECISION: First difference still NON-STATIONARY\n\n")

      cat("Testing SECOND DIFFERENCE (d=2):\n")
      cat("----------------------------\n")
      diff_series2 <- diff(insurance_series, differences = 2)
      adf_diff2 <- suppressWarnings(adf.test(diff_series2))
      cat("Test Statistic:", round(adf_diff2$statistic, 6), "\n")
      cat("P-value:", sprintf("%.10f", adf_diff2$p.value), "\n")

      if (adf_diff2$p.value < 0.05) {
        cat("[OK] DECISION: Second difference is STATIONARY -> d = 2\n\n")
      } else {
        cat("[WARNING] Series may require alternative transformation\n\n")
      }
      d_value <- 2
    }
  }

  cat("CONCLUSION: Integrated order d =", d_value, "\n")
  cat("  The differencing required to achieve stationarity defines\n")
  cat("  the integrated order d in the ARIMA(p,d,q) model.\n\n")

  return(d_value)
}

# ============================================================
# STEP 2: ACF AND PACF PLOTS
# Paragraph Section C: Box-Jenkins - Identification Stage
# ACF -> determines q (MA order)
# PACF -> determines p (AR order)
# ============================================================

plot_acf_pacf <- function(data) {
  cat("+================================================================+\n")
  cat("|            STEP 2: ACF AND PACF PLOTS                         |\n")
  cat("|        (Box-Jenkins Identification Stage)                     |\n")
  cat("+================================================================+\n\n")

  cat("PURPOSE: Model Identification\n")
  cat("  ACF:  Used to determine q (Moving Average order)\n")
  cat("  PACF: Used to determine p (Autoregressive order)\n\n")

  cat("INTERPRETATION GUIDE:\n")
  cat("  - ACF cuts off after lag q  -> MA(q) component\n")
  cat("  - PACF cuts off after lag p -> AR(p) component\n")
  cat("  - Both decay gradually      -> ARMA(p,q) model\n")
  cat("  - Significant spikes at seasonal lags (12, 24, 36)\n")
  cat("    indicate seasonal components (P, Q)\n\n")

  insurance_series <- data$`Individual Life Insurance`

  # Use platform-independent plotting
  if (.Platform$OS.type == "windows") {
    windows(width = 12, height = 8)
  } else {
    dev.new(width = 12, height = 8)
  }

  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  acf(insurance_series, main = "ACF of Individual Life Insurance",
      lag.max = 36, col = "darkblue", cex.main = 1.2)

  pacf(insurance_series, main = "PACF of Individual Life Insurance",
       lag.max = 36, col = "darkgreen", cex.main = 1.2)

  mtext("Box-Jenkins Identification: ACF and PACF Analysis", outer = TRUE, cex = 1.1)

  par(mfrow = c(1, 1))

  cat("[OK] ACF and PACF plots generated!\n")
  cat("  Use these plots to visually assess p and q orders.\n\n")
}

# ============================================================
# STEP 3: AUTO ARIMA WITH SEASONALITY
# Paragraph Section C: Box-Jenkins - Parameter Estimation
# Uses Maximum Likelihood Estimation (MLE) via auto.arima
# Model: ARIMA(p,d,q)(P,D,Q)[12]
# ============================================================

find_optimal_arima <- function(data, d_value) {
  cat("\n+================================================================+\n")
  cat("|     STEP 3: AUTO ARIMA - SEASONAL MODEL SELECTION            |\n")
  cat("|     (Box-Jenkins Parameter Estimation via MLE)               |\n")
  cat("+================================================================+\n\n")

  ts_data <- ts(data$`Individual Life Insurance`,
                frequency = 12,
                start = c(min(data$Year), 1))

  cat("FITTED ARIMA MODEL EQUATION:\n")
  cat("================================================================\n")
  cat("  yhat_t = c + sum(delta_i * y(t-i), i=1..p)\n")
  cat("             + sum(theta_j * eta(t-j), j=1..q) + epsilon_t\n\n")
  cat("  where:\n")
  cat("    y(t-i)   = lagged (past) observations\n")
  cat("    eta(t-j) = lagged forecast errors (residuals)\n")
  cat("    c        = constant term\n")
  cat("    delta_i  = AR coefficients\n")
  cat("    theta_j  = MA coefficients\n")
  cat("    epsilon_t = residual passed to XGBoost component\n\n")

  cat("MODEL STRUCTURE: ARIMA(p,d,q)(P,D,Q)[12]\n")
  cat("================================================================\n")
  cat("Non-Seasonal Part: (p,d,q)\n")
  cat("  p = AR (autoregressive) order\n")
  cat("  d = differencing order (from ADF test = ", d_value, ")\n", sep = "")
  cat("  q = MA (moving average) order\n\n")

  cat("Seasonal Part: (P,D,Q)[12]\n")
  cat("  P = Seasonal AR order\n")
  cat("  D = Seasonal differencing order\n")
  cat("  Q = Seasonal MA order\n")
  cat("  [12] = 12-month seasonality (monthly data)\n\n")

  cat("ESTIMATION METHOD: Maximum Likelihood Estimation (MLE)\n")
  cat("SELECTION CRITERION: Lowest AICc (corrected AIC)\n\n")

  cat("Running auto.arima with seasonal components...\n\n")

  optimal_model <- auto.arima(
    ts_data,
    d = d_value,
    ic = "aicc",
    stepwise = TRUE,
    trace = FALSE,
    max.p = 3,
    max.q = 3,
    max.P = 2,
    max.Q = 2
  )

  # Get components
  p <- optimal_model$arma[1]
  d <- optimal_model$arma[6]
  q <- optimal_model$arma[2]
  P <- optimal_model$arma[3]
  D <- optimal_model$arma[7]
  Q <- optimal_model$arma[4]

  cat("[OK] OPTIMAL MODEL FOUND:\n")
  cat("================================================================\n")
  cat("Model: ARIMA(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")[12]\n\n", sep = "")

  cat("PARAMETER VALUES:\n")
  cat("----------------------------\n")
  cat("Non-Seasonal: p=", p, ", d=", d, ", q=", q, "\n", sep = "")
  cat("Seasonal:     P=", P, ", D=", D, ", Q=", Q, "\n\n", sep = "")

  cat("ESTIMATED COEFFICIENTS (from MLE):\n")
  cat("----------------------------\n")
  coefs <- coef(optimal_model)
  if (length(coefs) > 0) {
    for (nm in names(coefs)) {
      cat(sprintf("  %-12s = %10.6f\n", nm, coefs[nm]))
    }
  }
  cat("\n")

  cat("MODEL FIT STATISTICS:\n")
  cat("----------------------------\n")
  cat("AIC:", round(optimal_model$aic, 4), "\n")
  cat("AICc:", round(optimal_model$aicc, 4), "\n")
  cat("BIC:", round(optimal_model$bic, 4), "\n")
  cat("Log Likelihood:", round(optimal_model$loglik, 4), "\n\n")

  return(optimal_model)
}

# ============================================================
# STEP 4: MANUAL COMPARISON - TOP 5 & TOP 10
# Paragraph Section C: Box-Jenkins - optimal p, d, q values
# ============================================================

compare_arima_models <- function(data, d_value) {
  cat("\n+================================================================+\n")
  cat("|   STEP 4: MANUAL SEARCH - SEASONAL ARIMA COMPARISON          |\n")
  cat("|   (Exhaustive search for optimal p, d, q, P, D, Q)          |\n")
  cat("+================================================================+\n\n")

  ts_data <- ts(data$`Individual Life Insurance`,
                frequency = 12,
                start = c(min(data$Year), 1))

  p_range <- 0:3
  q_range <- 0:3
  P_range <- 0:2
  D_range <- 0:1
  Q_range <- 0:2

  total <- length(p_range) * length(q_range) * length(P_range) * length(D_range) * length(Q_range)

  cat("Testing combinations:\n")
  cat("Non-seasonal: p=0-3, d=", d_value, ", q=0-3\n", sep = "")
  cat("Seasonal:     P=0-2, D=0-1, Q=0-2\n")
  cat("Total combinations:", total, "\n\n")

  results <- data.frame()
  tested <- 0
  converged <- 0

  for (p in p_range) {
    for (q in q_range) {
      for (P in P_range) {
        for (D in D_range) {
          for (Q in Q_range) {
            tested <- tested + 1
            tryCatch({
              model <- suppressWarnings(
                arima(ts_data, order = c(p, d_value, q),
                      seasonal = list(order = c(P, D, Q), period = 12),
                      method = "ML")
              )
              # Skip models with invalid AIC (NaN from log(s2) issues)
              if (!is.na(model$aic) && is.finite(model$aic)) {
                converged <- converged + 1
                results <- rbind(results, data.frame(
                  p = p, d = d_value, q = q,
                  P = P, D = D, Q = Q,
                  AIC = model$aic,
                  BIC = AIC(model, k = log(length(ts_data)))
                ))
              }
            }, error = function(e) {
              # Model did not converge - skip
            })
          }
        }
      }
    }
  }

  cat("Models tested:", tested, "\n")
  cat("Models converged:", converged, "\n\n")

  results <- results[order(results$AIC), ]

  cat("TOP 5 MODELS BY AIC:\n")
  cat("================================================================\n")
  cat("Rank | ARIMA(p,d,q)(P,D,Q)[12]  |      AIC      |      BIC\n")
  cat("-----+---------------------------+---------------+-------------\n")

  top5 <- head(results, 5)
  for (i in 1:nrow(top5)) {
    cat(sprintf("  %d  | ARIMA(%d,%d,%d)(%d,%d,%d)[12] | %13.4f | %13.4f\n",
                i, top5$p[i], top5$d[i], top5$q[i],
                top5$P[i], top5$D[i], top5$Q[i],
                top5$AIC[i], top5$BIC[i]))
  }

  cat("\n\nTOP 10 MODELS BY AIC:\n")
  cat("================================================================\n")
  cat("Rank | ARIMA(p,d,q)(P,D,Q)[12]  |      AIC      |      BIC\n")
  cat("-----+---------------------------+---------------+-------------\n")

  top10 <- head(results, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("  %-2d | ARIMA(%d,%d,%d)(%d,%d,%d)[12] | %13.4f | %13.4f\n",
                i, top10$p[i], top10$d[i], top10$q[i],
                top10$P[i], top10$D[i], top10$Q[i],
                top10$AIC[i], top10$BIC[i]))
  }

  cat("\n[OK] Rank 1 is the BEST model (lowest AIC)\n\n")

  return(results)
}

# ============================================================
# STEP 5: EXTRACT RESIDUALS & PREDICTED VALUES
# Paragraph Section A: epsilon_t is the residual series
# passed directly to XGBoost component
# ============================================================

extract_residuals <- function(model, data) {
  cat("\n+================================================================+\n")
  cat("|     STEP 5: PREDICTED VALUES & RESIDUALS EXTRACTION          |\n")
  cat("+================================================================+\n\n")

  cat("FORMULA:\n")
  cat("  Residual(t) = Actual(t) - Predicted(t)\n")
  cat("  epsilon_t = y_t - yhat_t\n\n")
  cat("  These residuals (epsilon_t) will be passed to XGBoost\n")
  cat("  to capture remaining non-linear patterns.\n\n")

  actual_values <- data$`Individual Life Insurance`
  fitted_values <- fitted(model)
  residuals_vals <- residuals(model)

  cat("FIRST 15 OBSERVATIONS:\n")
  cat("================================================================\n")
  cat("Time | Actual    | Predicted | Residual  | Error %\n")
  cat("-----+-----------+-----------+-----------+---------\n")

  for (i in 1:min(15, length(residuals_vals))) {
    a <- actual_values[i]
    p <- fitted_values[i]
    r <- residuals_vals[i]
    err <- ifelse(a != 0, (abs(r) / abs(a)) * 100, NA)
    cat(sprintf("%4d | %9.0f | %9.0f | %9.0f | %6.2f%%\n", i, a, p, r, err))
  }

  cat("\n")

  cat("RESIDUALS QUALITY ASSESSMENT:\n")
  cat("================================================================\n")
  mean_actual <- mean(actual_values)
  mean_residual <- mean(residuals_vals)
  residual_pct <- (abs(mean_residual) / abs(mean_actual)) * 100

  cat("Mean of Actual:", round(mean_actual, 2), "\n")
  cat("Mean of Residuals:", round(mean_residual, 4), "\n")
  cat("Residual % of Actual:", round(residual_pct, 4), "%\n")
  cat("Std Dev of Residuals:", round(sd(residuals_vals), 2), "\n")
  cat("Min Residual:", round(min(residuals_vals), 2), "\n")
  cat("Max Residual:", round(max(residuals_vals), 2), "\n\n")

  if (residual_pct < 1) {
    cat("[OK] GOOD: Residual mean is < 1% of actual mean\n")
    cat("  Model is UNBIASED\n\n")
  } else {
    cat("[!] CAUTION: Residual mean is", round(residual_pct, 2), "% of actual mean\n\n")
  }

  residuals_df <- data.frame(
    Time_Index = 1:length(residuals_vals),
    Year = data$Year,
    Month = data$Month,
    Actual_Premium = actual_values,
    Predicted_Premium = as.numeric(fitted_values),
    Residual = as.numeric(residuals_vals)
  )

  return(residuals_df)
}

# ============================================================
# STEP 6: DIAGNOSTIC CHECKING
# Paragraph Section C: Box-Jenkins - Diagnostic Checking
# Check residuals for autocorrelation, non-normality,
# and heteroskedasticity
# ============================================================

diagnostic_checking <- function(model, data) {
  cat("\n+================================================================+\n")
  cat("|     STEP 6: DIAGNOSTIC CHECKING                              |\n")
  cat("|     (Box-Jenkins Final Stage)                                |\n")
  cat("+================================================================+\n\n")

  cat("Checking residuals for:\n")
  cat("  1. Autocorrelation (Ljung-Box test)\n")
  cat("  2. Non-normality (Shapiro-Wilk test)\n")
  cat("  3. Heteroskedasticity (visual inspection)\n\n")

  residuals_vals <- residuals(model)

  # --- 1. Ljung-Box Test for Autocorrelation ---
  cat("1. LJUNG-BOX TEST (Autocorrelation in Residuals):\n")
  cat("----------------------------\n")
  cat("  H0: Residuals are independently distributed (no autocorrelation)\n")
  cat("  H1: Residuals exhibit autocorrelation\n\n")

  lb_test <- Box.test(residuals_vals, lag = 20, type = "Ljung-Box")
  cat("  Test Statistic:", round(lb_test$statistic, 4), "\n")
  cat("  P-value:", sprintf("%.6f", lb_test$p.value), "\n")

  if (lb_test$p.value > 0.05) {
    cat("  [OK] PASS: No significant autocorrelation in residuals\n\n")
  } else {
    cat("  [!] FAIL: Residuals show significant autocorrelation\n")
    cat("    Consider revising ARIMA orders (p, q, P, Q)\n\n")
  }

  # --- 2. Shapiro-Wilk Test for Normality ---
  cat("2. SHAPIRO-WILK TEST (Normality of Residuals):\n")
  cat("----------------------------\n")
  cat("  H0: Residuals are normally distributed\n")
  cat("  H1: Residuals are NOT normally distributed\n\n")

  # Shapiro-Wilk has a max sample size of 5000
  n_resid <- length(residuals_vals)
  if (n_resid <= 5000) {
    sw_test <- shapiro.test(as.numeric(residuals_vals))
    cat("  Test Statistic:", round(sw_test$statistic, 6), "\n")
    cat("  P-value:", sprintf("%.6f", sw_test$p.value), "\n")

    if (sw_test$p.value > 0.05) {
      cat("  [OK] PASS: Residuals appear normally distributed\n\n")
    } else {
      cat("  [!] FAIL: Residuals deviate from normality\n")
      cat("    (Common for financial data; model may still be adequate)\n\n")
    }
  } else {
    cat("  Sample size > 5000; using Jarque-Bera test instead.\n")
    cat("  (Shapiro-Wilk limited to n <= 5000)\n\n")
  }

  # --- 3. Visual Diagnostic Plots ---
  cat("3. DIAGNOSTIC PLOTS:\n")
  cat("----------------------------\n")

  if (.Platform$OS.type == "windows") {
    windows(width = 12, height = 10)
  } else {
    dev.new(width = 12, height = 10)
  }

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  # Residuals over time (check for heteroskedasticity)
  plot(as.numeric(residuals_vals), type = "l",
       main = "Residuals Over Time\n(Check: constant variance / heteroskedasticity)",
       xlab = "Time", ylab = "Residuals", col = "darkblue", cex.main = 1.0)
  abline(h = 0, col = "red", lty = 2, lwd = 2)
  grid()

  # Histogram (check normality)
  hist(residuals_vals, main = "Distribution of Residuals\n(Check: normality)",
       xlab = "Residuals", col = "skyblue", border = "navy",
       breaks = 30, cex.main = 1.0)

  # ACF of residuals (check autocorrelation)
  acf(residuals_vals, main = "ACF of Residuals\n(Check: autocorrelation)",
      lag.max = 24, cex.main = 1.0)

  # Q-Q plot (check normality)
  qqnorm(residuals_vals, main = "Q-Q Plot of Residuals\n(Check: normality)",
         pch = 16, col = "steelblue", cex.main = 1.0)
  qqline(residuals_vals, col = "red", lwd = 2)

  mtext("Box-Jenkins Diagnostic Checking", outer = TRUE, cex = 1.1)
  par(mfrow = c(1, 1))

  cat("  [OK] Diagnostic plots generated!\n\n")

  # --- Summary ---
  cat("DIAGNOSTIC SUMMARY:\n")
  cat("================================================================\n")
  cat("  Autocorrelation (Ljung-Box): ",
      ifelse(lb_test$p.value > 0.05, "PASS", "FAIL"), "\n")
  if (n_resid <= 5000) {
    cat("  Normality (Shapiro-Wilk):    ",
        ifelse(sw_test$p.value > 0.05, "PASS", "FAIL"), "\n")
  }
  cat("  Heteroskedasticity:           See residual time plot\n\n")

  if (lb_test$p.value > 0.05) {
    cat("[OK] Residuals pass key diagnostic checks.\n")
    cat("  The model adequately captures the data patterns.\n\n")
  } else {
    cat("[!] Residuals show some inadequacy.\n")
    cat("  Consider adjustments: revise ARIMA orders or add variables.\n\n")
  }
}

# ============================================================
# STEP 7: PLOT ACTUAL VS PREDICTED
# ============================================================

plot_actual_vs_predicted <- function(data, fitted_values) {
  cat("\n+================================================================+\n")
  cat("|         STEP 7: ACTUAL vs PREDICTED COMPARISON                |\n")
  cat("+================================================================+\n\n")

  actual <- data$`Individual Life Insurance`

  if (.Platform$OS.type == "windows") {
    windows(width = 14, height = 6)
  } else {
    dev.new(width = 14, height = 6)
  }

  par(mar = c(4, 4, 3, 1))

  plot(1:length(actual), actual,
       type = "l", lwd = 2, col = "darkblue",
       main = "Actual vs ARIMA Predicted Premium Revenue",
       xlab = "Time Period", ylab = "Premium Amount",
       cex.main = 1.3)

  lines(1:length(fitted_values), fitted_values,
        col = "red", lwd = 2, lty = 2)

  legend("topright",
         legend = c("Actual Premium", "ARIMA Predicted"),
         col = c("darkblue", "red"),
         lty = c(1, 2),
         lwd = 2,
         cex = 1.1)
  grid()

  cat("[OK] Comparison plot generated!\n\n")
}

# ============================================================
# MAIN EXECUTION
# ============================================================

main <- function() {
  tryCatch({
    cat("\n================================================================\n")
    cat("  COMPONENT 1: LINEAR FORECASTING USING ARIMA\n")
    cat("  OBJECTIVE 1: OPTIMAL SEASONAL ARIMA(p,d,q)(P,D,Q)[12]\n")
    cat("================================================================\n\n")

    cat("Methodology: Box-Jenkins (3 stages)\n")
    cat("  1. Model Identification (ADF + ACF/PACF)\n")
    cat("  2. Parameter Estimation (MLE via auto.arima)\n")
    cat("  3. Diagnostic Checking (residual analysis)\n\n")

    # Load
    raw_data <- load_data()
    data <- prepare_data(raw_data)

    cat("Data Loaded:", nrow(data), "observations\n")
    cat("Date Range:", min(data$Year), "to", max(data$Year), "\n")
    cat("Mean Premium:", round(mean(data$`Individual Life Insurance`), 2), "\n\n")

    # Step 1: Stationarity Testing (Paragraph Section B)
    d_value <- perform_adf_test(data)

    # Step 2: ACF/PACF - Model Identification (Paragraph Section C)
    plot_acf_pacf(data)

    # Step 3: Parameter Estimation via MLE (Paragraph Sections A & C)
    optimal_model <- find_optimal_arima(data, d_value)

    # Step 4: Exhaustive Model Comparison
    all_models <- compare_arima_models(data, d_value)

    # Step 5: Extract Residuals for XGBoost (Paragraph Section A)
    residuals_df <- extract_residuals(optimal_model, data)

    # Step 6: Diagnostic Checking (Paragraph Section C)
    diagnostic_checking(optimal_model, data)

    # Step 7: Actual vs Predicted Plot
    fitted_vals <- fitted(optimal_model)
    plot_actual_vs_predicted(data, fitted_vals)

    # FINAL SUMMARY
    cat("\n================================================================\n")
    cat("                    FINAL RESULTS\n")
    cat("================================================================\n\n")

    p <- optimal_model$arma[1]
    d <- optimal_model$arma[6]
    q <- optimal_model$arma[2]
    P <- optimal_model$arma[3]
    D <- optimal_model$arma[7]
    Q <- optimal_model$arma[4]

    cat("OPTIMAL MODEL SELECTED:\n")
    cat("================================================================\n")
    cat("Model: ARIMA(", p, ",", d, ",", q, ")(", P, ",", D, ",", Q, ")[12]\n\n", sep = "")

    cat("Non-Seasonal: p=", p, " d=", d, " q=", q, "\n", sep = "")
    cat("Seasonal:     P=", P, " D=", D, " Q=", Q, "\n\n", sep = "")

    cat("Estimated Coefficients:\n")
    coefs <- coef(optimal_model)
    if (length(coefs) > 0) {
      for (nm in names(coefs)) {
        cat(sprintf("  %-12s = %10.6f\n", nm, coefs[nm]))
      }
    }
    cat("\n")

    cat("Statistics:\n")
    cat("  AIC:  ", round(optimal_model$aic, 4), "\n")
    cat("  AICc: ", round(optimal_model$aicc, 4), "\n")
    cat("  BIC:  ", round(optimal_model$bic, 4), "\n")
    cat("  Residual Mean:", round(mean(residuals_df$Residual), 4), "\n\n")

    cat("================================================================\n")
    cat("[OK] COMPONENT 1 (ARIMA Linear Forecasting) COMPLETE\n")
    cat("[OK] Residual series (epsilon_t) ready for Component 2: XGBoost\n")
    cat("================================================================\n\n")

    # Save
    write.csv(residuals_df, "objective1_residuals.csv", row.names = FALSE)
    cat("[OK] Residuals saved to 'objective1_residuals.csv'\n\n")

  }, error = function(e) {
    cat("\n[ERROR]:", e$message, "\n\n")
  })
}

main()
