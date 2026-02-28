library(xgboost)
library(readxl)
library(dplyr)

# ============================================================
# COMPONENT 2: NON-LINEAR RESIDUAL MODELING USING XGBOOST
# ============================================================
# eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)       (Eq 3.4)
# Objective = sum(Loss) + sum(Omega(fk))            (Eq 3.5)
# Omega(fk) = gamma*T + 0.5*lambda*sum(wk^2)       (Eq 3.6)
# eta_hat_t = sum(eta_hat_t(k), k=1..K)            (Eq 3.7)
# ============================================================

# --- Helper: load file via file.choose() ---
load_file <- function(prompt_msg) {
  cat("\n", prompt_msg, "\n")
  file_path <- file.choose()
  ext <- tolower(tools::file_ext(file_path))
  if (ext %in% c("xlsx", "xls")) return(read_excel(file_path))
  if (ext == "csv") return(read.csv(file_path))
  stop("Unsupported format. Use .xlsx, .xls, or .csv")
}

# ============================================================
# STEP 1: LOAD ARIMA RESIDUALS
# ============================================================
load_residuals <- function() {
  cat("\n+================================================================+\n")
  cat("|  STEP 1: LOAD ARIMA RESIDUALS                                |\n")
  cat("+================================================================+\n\n")

  if (file.exists("objective1_residuals.csv")) {
    df <- read.csv("objective1_residuals.csv")
    cat("[OK] Loaded from 'objective1_residuals.csv'\n")
  } else {
    df <- load_file("Select the ARIMA residuals CSV file...")
  }
  cat("  Observations:", nrow(df), " | Year range:", min(df$Year), "-", max(df$Year), "\n")
  cat("  Residual mean:", round(mean(df$Residual), 4), " | SD:", round(sd(df$Residual), 2), "\n\n")
  return(df)
}

# ============================================================
# STEP 2: LOAD & PREPARE MACROECONOMIC DATA
# ============================================================
load_and_prepare_gdp <- function() {
  cat("+================================================================+\n")
  cat("|  STEP 2A: LOAD GDP (Yearly -> Monthly via Spline)             |\n")
  cat("+================================================================+\n\n")

  gdp_raw <- load_file("Select your GDP data file (yearly)...")
  col_names <- names(gdp_raw)
  col_lower <- tolower(col_names)
  year_col <- col_names[grep("year", col_lower)[1]]
  gdp_col <- col_names[grep("gdp", col_lower)[1]]
  if (is.na(year_col) || is.na(gdp_col)) { year_col <- col_names[1]; gdp_col <- col_names[2] }

  gdp <- data.frame(Year = as.numeric(gdp_raw[[year_col]]), GDP = as.numeric(gdp_raw[[gdp_col]]))
  gdp <- gdp[complete.cases(gdp), ]
  gdp <- gdp[order(gdp$Year), ]

  # Spline interpolation: yearly -> monthly
  yearly_pts <- gdp$Year + 0.5
  min_yr <- min(gdp$Year); max_yr <- max(gdp$Year)
  monthly_pts <- seq(min_yr + (1 - 0.5)/12, max_yr + (12 - 0.5)/12, by = 1/12)
  sp <- spline(yearly_pts, gdp$GDP, xout = monthly_pts, method = "natural")

  gdp_m <- data.frame(Year = rep(min_yr:max_yr, each = 12), Month = rep(1:12, times = max_yr - min_yr + 1))
  gdp_m <- gdp_m[1:length(sp$y), ]
  gdp_m$GDP <- sp$y
  cat("[OK] GDP interpolated:", nrow(gdp_m), "monthly obs\n\n")
  return(gdp_m)
}

load_and_prepare_cpi <- function() {
  cat("+================================================================+\n")
  cat("|  STEP 2B: LOAD CPI / INFLATION (Monthly)                     |\n")
  cat("+================================================================+\n\n")

  cpi_raw <- load_file("Select your CPI data file (monthly)...")
  col_names <- names(cpi_raw)
  col_lower <- tolower(col_names)
  year_col <- col_names[grep("year", col_lower)[1]]
  month_col <- col_names[grep("month", col_lower)[1]]
  cpi_col <- col_names[grep("cpi|inflation|index", col_lower)[1]]
  if (is.na(year_col) || is.na(month_col) || is.na(cpi_col)) {
    year_col <- col_names[1]; month_col <- col_names[2]; cpi_col <- col_names[3]
  }

  cpi <- data.frame(
    Year = as.numeric(cpi_raw[[year_col]]),
    Month = as.numeric(cpi_raw[[month_col]]),
    CPI = as.numeric(cpi_raw[[cpi_col]])
  )
  cpi <- cpi[complete.cases(cpi), ]
  cpi <- cpi[order(cpi$Year, cpi$Month), ]
  cat("[OK] CPI loaded:", nrow(cpi), "obs\n\n")
  return(cpi)
}

# ============================================================
# STEP 3: MERGE, FEATURE ENGINEERING, ALIGN
# ============================================================
merge_and_engineer_features <- function(residuals_df, gdp_monthly, cpi_data) {
  cat("+================================================================+\n")
  cat("|  STEP 3: MERGE & FEATURE ENGINEERING                         |\n")
  cat("+================================================================+\n\n")

  # Inner join on Year + Month
  merged <- merge(residuals_df, gdp_monthly, by = c("Year", "Month"))
  merged <- merge(merged, cpi_data, by = c("Year", "Month"))
  merged <- merged[order(merged$Year, merged$Month), ]

  cat("  After merge:", nrow(merged), "rows (common date range)\n\n")

  # --- FEATURE ENGINEERING ---
  # These derived features give XGBoost the non-linear signal
  # that raw GDP/CPI levels alone cannot provide
  n <- nrow(merged)

  # GDP growth rate (month-over-month change from interpolated GDP)
  merged$GDP_Growth <- c(NA, diff(merged$GDP) / merged$GDP[-n] * 100)

  # CPI inflation rate (month-over-month % change)
  merged$CPI_Change <- c(NA, diff(merged$CPI) / merged$CPI[-n] * 100)

  # Year-over-year CPI change (if enough data)
  merged$CPI_YoY <- NA
  if (n > 12) {
    merged$CPI_YoY[13:n] <- (merged$CPI[13:n] - merged$CPI[1:(n-12)]) / merged$CPI[1:(n-12)] * 100
  }

  # Lagged residuals (autoregressive features for XGBoost)
  merged$Resid_Lag1 <- c(NA, merged$Residual[-n])
  merged$Resid_Lag2 <- c(NA, NA, merged$Residual[-c(n-1, n)])
  merged$Resid_Lag3 <- c(NA, NA, NA, merged$Residual[-c(n-2, n-1, n)])

  # Month indicator (captures residual seasonality ARIMA may have missed)
  merged$Month_Num <- merged$Month

  # GDP level (still useful as context)
  # CPI level (still useful as context)
  # Both already in merged

  # Remove rows with NA from lagging
  merged <- merged[complete.cases(merged), ]

  feature_cols <- c("GDP", "CPI", "GDP_Growth", "CPI_Change", "CPI_YoY",
                     "Resid_Lag1", "Resid_Lag2", "Resid_Lag3", "Month_Num")

  cat("ENGINEERED FEATURES (xi):\n")
  cat("----------------------------\n")
  cat("  GDP          - Interpolated monthly GDP level\n")
  cat("  CPI          - Monthly Consumer Price Index\n")
  cat("  GDP_Growth   - Month-over-month GDP growth rate (%)\n")
  cat("  CPI_Change   - Month-over-month CPI change (%)\n")
  cat("  CPI_YoY      - Year-over-year CPI inflation (%)\n")
  cat("  Resid_Lag1   - Previous month residual\n")
  cat("  Resid_Lag2   - 2-month lagged residual\n")
  cat("  Resid_Lag3   - 3-month lagged residual\n")
  cat("  Month_Num    - Month indicator (residual seasonality)\n\n")
  cat("  Total features:", length(feature_cols), "\n")
  cat("  Usable observations:", nrow(merged), "\n\n")

  return(list(data = merged, feature_cols = feature_cols))
}

# ============================================================
# STEP 4: TRAIN/TEST SPLIT (time-based)
# ============================================================
split_data <- function(merged_data, train_ratio = 0.8) {
  cat("+================================================================+\n")
  cat("|  STEP 4: TRAIN/TEST SPLIT                                    |\n")
  cat("+================================================================+\n\n")

  n <- nrow(merged_data)
  train_size <- floor(train_ratio * n)
  train <- merged_data[1:train_size, ]
  test <- merged_data[(train_size + 1):n, ]

  cat("  Train:", train_size, "obs | Test:", n - train_size, "obs\n")
  cat("  Train:", min(train$Year), "M", min(train$Month[train$Year == min(train$Year)]),
      "to", max(train$Year), "M", max(train$Month[train$Year == max(train$Year)]), "\n")
  cat("  Test: ", min(test$Year), "M", min(test$Month[test$Year == min(test$Year)]),
      "to", max(test$Year), "M", max(test$Month[test$Year == max(test$Year)]), "\n\n")

  return(list(train = train, test = test))
}

# ============================================================
# STEP 5: XGBOOST TRAINING
# ============================================================
train_xgboost <- function(train_data, test_data, feature_cols) {
  cat("+================================================================+\n")
  cat("|  STEP 5: XGBOOST MODEL TRAINING                              |\n")
  cat("+================================================================+\n\n")

  cat("Prediction: eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)  (Eq 3.4)\n")
  cat("Objective:  L = sum(Loss) + sum(Omega(fk))              (Eq 3.5)\n")
  cat("Omega(fk) = gamma*T + 0.5*lambda*sum(wk^2)             (Eq 3.6)\n\n")

  train_X <- as.matrix(train_data[, feature_cols])
  train_y <- train_data$Residual
  test_X <- as.matrix(test_data[, feature_cols])
  test_y <- test_data$Residual

  dtrain <- xgb.DMatrix(data = train_X, label = train_y)
  dtest <- xgb.DMatrix(data = test_X, label = test_y)

  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse",
    eta = 0.05,              # Lower learning rate for better generalization
    max_depth = 3,           # Shallow trees to prevent overfitting
    min_child_weight = 5,
    subsample = 0.8,
    colsample_bytree = 0.8,
    gamma = 0.1,             # Eq 3.6: min loss reduction for split
    lambda = 1.5,            # Eq 3.6: L2 regularization
    alpha = 0.1              # L1 regularization
  )

  cat("Hyperparameters:\n")
  cat("  eta=", params$eta, " max_depth=", params$max_depth,
      " gamma=", params$gamma, " lambda=", params$lambda, "\n\n")

  # Cross-validation for optimal K (number of trees)
  set.seed(42)
  cv <- xgb.cv(params = params, data = dtrain, nrounds = 1000,
               nfold = 5, early_stopping_rounds = 50, verbose = 0)

  best_K <- cv$best_iteration
  best_cv_rmse <- cv$evaluation_log$test_rmse_mean[best_K]
  cat("  CV optimal K (trees):", best_K, " | CV RMSE:", round(best_cv_rmse, 4), "\n\n")

  # Train final model
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = best_K,
                          watchlist = list(train = dtrain, test = dtest), verbose = 0)

  cat("[OK] XGBoost trained with K =", best_K, "trees\n\n")

  return(list(model = xgb_model, dtrain = dtrain, dtest = dtest,
              train_y = train_y, test_y = test_y,
              feature_cols = feature_cols, best_K = best_K))
}

# ============================================================
# STEP 6: EVALUATION & COMPARISON
# ============================================================
evaluate_and_compare <- function(xgb_result, train_data, test_data) {
  cat("+================================================================+\n")
  cat("|  STEP 6: EVALUATION & ARIMA vs HYBRID COMPARISON             |\n")
  cat("+================================================================+\n\n")

  model <- xgb_result$model
  train_pred <- predict(model, xgb_result$dtrain)
  test_pred <- predict(model, xgb_result$dtest)

  # --- XGBoost residual prediction metrics ---
  train_rmse <- sqrt(mean((xgb_result$train_y - train_pred)^2))
  test_rmse <- sqrt(mean((xgb_result$test_y - test_pred)^2))
  test_mae <- mean(abs(xgb_result$test_y - test_pred))
  ss_res <- sum((xgb_result$test_y - test_pred)^2)
  ss_tot <- sum((xgb_result$test_y - mean(xgb_result$test_y))^2)
  test_r2 <- 1 - (ss_res / ss_tot)

  cat("XGBoost Residual Prediction (eta_hat_t):\n")
  cat("  Train RMSE:", round(train_rmse, 2), "\n")
  cat("  Test RMSE: ", round(test_rmse, 2), "\n")
  cat("  Test MAE:  ", round(test_mae, 2), "\n")
  cat("  Test R2:   ", round(test_r2, 4), "\n\n")

  # --- HYBRID vs ARIMA-only comparison ---
  # Hybrid = ARIMA_Predicted + XGBoost_Predicted_Residual
  train_hybrid <- train_data$Predicted_Premium + train_pred
  test_hybrid <- test_data$Predicted_Premium + test_pred

  # ARIMA-only errors
  train_arima_rmse <- sqrt(mean((train_data$Actual_Premium - train_data$Predicted_Premium)^2))
  test_arima_rmse <- sqrt(mean((test_data$Actual_Premium - test_data$Predicted_Premium)^2))
  test_arima_mae <- mean(abs(test_data$Actual_Premium - test_data$Predicted_Premium))
  arima_ss_res <- sum((test_data$Actual_Premium - test_data$Predicted_Premium)^2)
  arima_ss_tot <- sum((test_data$Actual_Premium - mean(test_data$Actual_Premium))^2)
  test_arima_r2 <- 1 - (arima_ss_res / arima_ss_tot)

  # Hybrid errors
  train_hybrid_rmse <- sqrt(mean((train_data$Actual_Premium - train_hybrid)^2))
  test_hybrid_rmse <- sqrt(mean((test_data$Actual_Premium - test_hybrid)^2))
  test_hybrid_mae <- mean(abs(test_data$Actual_Premium - test_hybrid))
  hybrid_ss_res <- sum((test_data$Actual_Premium - test_hybrid)^2)
  test_hybrid_r2 <- 1 - (hybrid_ss_res / arima_ss_tot)

  improvement_rmse <- ((test_arima_rmse - test_hybrid_rmse) / test_arima_rmse) * 100
  improvement_mae <- ((test_arima_mae - test_hybrid_mae) / test_arima_mae) * 100

  cat("MODEL COMPARISON (Test Set):\n")
  cat("================================================================\n")
  cat(sprintf("  %-25s | %12s | %12s | %10s\n", "Metric", "ARIMA Only", "Hybrid", "Improvement"))
  cat("  --------------------------+-------------+-------------+-----------\n")
  cat(sprintf("  %-25s | %12.2f | %12.2f | %9.2f%%\n", "RMSE", test_arima_rmse, test_hybrid_rmse, improvement_rmse))
  cat(sprintf("  %-25s | %12.2f | %12.2f | %9.2f%%\n", "MAE", test_arima_mae, test_hybrid_mae, improvement_mae))
  cat(sprintf("  %-25s | %12.4f | %12.4f |\n", "R-squared", test_arima_r2, test_hybrid_r2))
  cat("\n")

  if (improvement_rmse > 0) {
    cat("[OK] Hybrid model IMPROVES on ARIMA by", round(improvement_rmse, 2), "% RMSE reduction\n\n")
  } else {
    cat("[!] Hybrid model did not improve. Consider adding more features.\n\n")
  }

  return(list(
    train_pred = train_pred, test_pred = test_pred,
    train_hybrid = train_hybrid, test_hybrid = test_hybrid,
    test_arima_rmse = test_arima_rmse, test_hybrid_rmse = test_hybrid_rmse,
    test_arima_r2 = test_arima_r2, test_hybrid_r2 = test_hybrid_r2,
    improvement_rmse = improvement_rmse
  ))
}

# ============================================================
# STEP 7: FEATURE IMPORTANCE
# ============================================================
show_importance <- function(xgb_result) {
  cat("+================================================================+\n")
  cat("|  STEP 7: FEATURE IMPORTANCE                                   |\n")
  cat("+================================================================+\n\n")

  imp <- xgb.importance(feature_names = xgb_result$feature_cols, model = xgb_result$model)
  for (i in 1:nrow(imp)) {
    cat(sprintf("  %-15s Gain: %.4f  Cover: %.4f  Freq: %.4f\n",
                imp$Feature[i], imp$Gain[i], imp$Cover[i], imp$Frequency[i]))
  }
  cat("\n")

  if (.Platform$OS.type == "windows") { windows(width = 10, height = 6) } else { dev.new(width = 10, height = 6) }
  barplot(imp$Gain, names.arg = imp$Feature, main = "XGBoost Feature Importance (Gain)",
          ylab = "Gain", col = "steelblue", border = "navy", las = 2, cex.names = 0.8)
  cat("[OK] Feature importance plot generated!\n\n")
  return(imp)
}

# ============================================================
# STEP 8: DIAGNOSTIC PLOTS
# ============================================================
plot_diagnostics <- function(xgb_result, eval_result, train_data, test_data) {
  cat("+================================================================+\n")
  cat("|  STEP 8: DIAGNOSTIC & COMPARISON PLOTS                       |\n")
  cat("+================================================================+\n\n")

  if (.Platform$OS.type == "windows") { windows(width = 14, height = 10) } else { dev.new(width = 14, height = 10) }
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  # 1. Actual vs Predicted residuals (scatter)
  plot(xgb_result$test_y, eval_result$test_pred, pch = 16, col = "steelblue",
       main = "Actual vs Predicted Residuals (Test)", xlab = "Actual eta_t", ylab = "Predicted eta_hat_t")
  abline(0, 1, col = "red", lwd = 2, lty = 2); grid()

  # 2. XGBoost prediction errors
  errs <- xgb_result$test_y - eval_result$test_pred
  plot(errs, type = "l", col = "darkblue", main = "XGBoost Prediction Errors",
       xlab = "Test Observation", ylab = "Error"); abline(h = 0, col = "red", lty = 2); grid()

  # 3. Full timeline: Actual vs ARIMA vs Hybrid
  n_train <- nrow(train_data); n_test <- nrow(test_data); n_all <- n_train + n_test
  all_actual <- c(train_data$Actual_Premium, test_data$Actual_Premium)
  all_arima <- c(train_data$Predicted_Premium, test_data$Predicted_Premium)
  all_hybrid <- c(eval_result$train_hybrid, eval_result$test_hybrid)

  plot(1:n_all, all_actual, type = "l", lwd = 2, col = "darkblue",
       main = "Actual vs ARIMA vs Hybrid", xlab = "Time", ylab = "Premium")
  lines(1:n_all, all_arima, col = "orange", lwd = 1.5, lty = 2)
  lines(1:n_all, all_hybrid, col = "red", lwd = 2)
  abline(v = n_train + 0.5, col = "gray40", lty = 3, lwd = 2)
  legend("topright", c("Actual", "ARIMA", "Hybrid", "Split"),
         col = c("darkblue", "orange", "red", "gray40"), lty = c(1,2,1,3), lwd = 2, cex = 0.7)
  grid()

  # 4. Error distribution
  hist(errs, main = "Distribution of XGBoost Errors", xlab = "Error",
       col = "skyblue", border = "navy", breaks = 20)

  mtext("XGBoost Residual Model Diagnostics", outer = TRUE, cex = 1.1)
  par(mfrow = c(1, 1))
  cat("[OK] Diagnostic plots generated!\n\n")
}

# ============================================================
# STEP 9: SAVE RESULTS
# ============================================================
save_results <- function(eval_result, train_data, test_data) {
  cat("+================================================================+\n")
  cat("|  STEP 9: SAVE RESULTS                                        |\n")
  cat("+================================================================+\n\n")

  out <- rbind(
    data.frame(Year = train_data$Year, Month = train_data$Month,
               Actual_Premium = train_data$Actual_Premium,
               ARIMA_Predicted = train_data$Predicted_Premium,
               ARIMA_Residual = train_data$Residual,
               XGBoost_Predicted_Residual = eval_result$train_pred,
               Hybrid_Predicted = eval_result$train_hybrid, Set = "Train"),
    data.frame(Year = test_data$Year, Month = test_data$Month,
               Actual_Premium = test_data$Actual_Premium,
               ARIMA_Predicted = test_data$Predicted_Premium,
               ARIMA_Residual = test_data$Residual,
               XGBoost_Predicted_Residual = eval_result$test_pred,
               Hybrid_Predicted = eval_result$test_hybrid, Set = "Test")
  )

  write.csv(out, "objective2_hybrid_results.csv", row.names = FALSE)
  cat("[OK] Results saved to 'objective2_hybrid_results.csv'\n\n")
  return(out)
}

# ============================================================
# MAIN EXECUTION
# ============================================================
main <- function() {
  tryCatch({
    cat("\n================================================================\n")
    cat("  COMPONENT 2: NON-LINEAR RESIDUAL MODELING USING XGBOOST\n")
    cat("  OBJECTIVE 2: Train XGBoost on ARIMA Residuals with GDP & CPI\n")
    cat("================================================================\n\n")

    # Step 1: Load ARIMA residuals
    residuals_df <- load_residuals()

    # Step 2: Load macroeconomic data
    gdp_monthly <- load_and_prepare_gdp()
    cpi_data <- load_and_prepare_cpi()

    # Step 3: Merge & feature engineering
    result <- merge_and_engineer_features(residuals_df, gdp_monthly, cpi_data)
    merged_data <- result$data
    feature_cols <- result$feature_cols

    # Step 4: Train/test split
    splits <- split_data(merged_data, train_ratio = 0.8)

    # Step 5: Train XGBoost
    xgb_result <- train_xgboost(splits$train, splits$test, feature_cols)

    # Step 6: Evaluate & compare ARIMA vs Hybrid
    eval_result <- evaluate_and_compare(xgb_result, splits$train, splits$test)

    # Step 7: Feature importance
    show_importance(xgb_result)

    # Step 8: Plots
    plot_diagnostics(xgb_result, eval_result, splits$train, splits$test)

    # Step 9: Save
    save_results(eval_result, splits$train, splits$test)

    # Final summary
    cat("================================================================\n")
    cat("  FINAL: Hybrid reduced RMSE by", round(eval_result$improvement_rmse, 2), "%\n")
    cat("  ARIMA R2:", round(eval_result$test_arima_r2, 4),
        " | Hybrid R2:", round(eval_result$test_hybrid_r2, 4), "\n")
    cat("[OK] OBJECTIVE 2 COMPLETE\n")
    cat("================================================================\n\n")

  }, error = function(e) {
    cat("\n[ERROR]:", e$message, "\n\n")
  })
}

main()
