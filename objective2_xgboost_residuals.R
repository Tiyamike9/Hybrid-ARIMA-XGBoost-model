library(xgboost)
library(readxl)
library(dplyr)

# ============================================================
# COMPONENT 2: NON-LINEAR RESIDUAL MODELING USING XGBOOST
# ============================================================
# The residuals (eta_t) from the fitted ARIMA model contain
# the unexplained, non-linear information that the linear
# model failed to capture. XGBoost is employed to predict
# these residuals (eta_hat_t) using macroeconomic variables
# (xi) as input features.
#
# Prediction equation:
#   eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)
#
# Where K = total number of trees, fk(xi) = score from k-th tree
#
# Macroeconomic Variables (xi): GDP, CPI (Inflation)
# Note: Interest rate data not yet available
# ============================================================

# ============================================================
# HELPER: Load a file via file.choose()
# ============================================================

load_file <- function(prompt_msg) {
  cat("\n", prompt_msg, "\n")
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

# ============================================================
# STEP 1: LOAD ARIMA RESIDUALS FROM OBJECTIVE 1
# ============================================================

load_residuals <- function() {
  cat("\n+================================================================+\n")
  cat("|  STEP 1: LOAD ARIMA RESIDUALS FROM OBJECTIVE 1               |\n")
  cat("+================================================================+\n\n")

  cat("The residuals (eta_t) from the ARIMA model represent the\n")
  cat("unexplained non-linear patterns that XGBoost will model.\n\n")

  # Try to load from default path first, otherwise ask user
  if (file.exists("objective1_residuals.csv")) {
    residuals_df <- read.csv("objective1_residuals.csv")
    cat("[OK] Loaded residuals from 'objective1_residuals.csv'\n")
  } else {
    residuals_df <- load_file("Select the ARIMA residuals CSV file (objective1_residuals.csv)...")
    cat("[OK] Loaded residuals file\n")
  }

  cat("  Observations:", nrow(residuals_df), "\n")
  cat("  Columns:", paste(names(residuals_df), collapse = ", "), "\n")
  cat("  Year range:", min(residuals_df$Year), "to", max(residuals_df$Year), "\n")
  cat("  Residual mean:", round(mean(residuals_df$Residual), 4), "\n")
  cat("  Residual std:", round(sd(residuals_df$Residual), 2), "\n\n")

  return(residuals_df)
}

# ============================================================
# STEP 2: LOAD & PREPARE MACROECONOMIC DATA
# Handles: yearly GDP -> monthly interpolation,
#           CPI date mismatch, missing interest rates
# ============================================================

load_gdp_data <- function() {
  cat("\n+================================================================+\n")
  cat("|  STEP 2A: LOAD GDP DATA (Yearly -> Monthly Interpolation)    |\n")
  cat("+================================================================+\n\n")

  cat("GDP data is yearly. It will be interpolated to monthly\n")
  cat("frequency using cubic spline interpolation.\n\n")

  gdp_raw <- load_file("Select your GDP data file (yearly)...")

  cat("Raw GDP data preview:\n")
  print(head(gdp_raw))
  cat("\n")

  return(gdp_raw)
}

interpolate_gdp_to_monthly <- function(gdp_raw) {
  # Identify the year and GDP columns
  # Try common column name patterns
  col_names <- names(gdp_raw)
  col_lower <- tolower(col_names)

  year_col <- col_names[grep("year", col_lower)[1]]
  gdp_col <- col_names[grep("gdp", col_lower)[1]]

  if (is.na(year_col) || is.na(gdp_col)) {
    cat("[!] Could not auto-detect columns. Using first two columns.\n")
    cat("    Column 1 (Year):", col_names[1], "\n")
    cat("    Column 2 (GDP):", col_names[2], "\n\n")
    year_col <- col_names[1]
    gdp_col <- col_names[2]
  } else {
    cat("  Detected Year column:", year_col, "\n")
    cat("  Detected GDP column:", gdp_col, "\n\n")
  }

  gdp_clean <- data.frame(
    Year = as.numeric(gdp_raw[[year_col]]),
    GDP = as.numeric(gdp_raw[[gdp_col]])
  )
  gdp_clean <- gdp_clean[complete.cases(gdp_clean), ]
  gdp_clean <- gdp_clean[order(gdp_clean$Year), ]

  cat("  GDP yearly range:", min(gdp_clean$Year), "to", max(gdp_clean$Year), "\n")
  cat("  GDP values:", nrow(gdp_clean), "yearly observations\n\n")

  # Cubic spline interpolation: yearly -> monthly
  # Place yearly GDP at mid-year (month 6) for interpolation
  yearly_points <- gdp_clean$Year + 0.5  # mid-year

  # Create monthly time points for all years
  min_year <- min(gdp_clean$Year)
  max_year <- max(gdp_clean$Year)
  monthly_times <- seq(min_year + (1 - 0.5) / 12,
                       max_year + (12 - 0.5) / 12,
                       by = 1 / 12)

  # Spline interpolation
  spline_fit <- spline(yearly_points, gdp_clean$GDP,
                       xout = monthly_times, method = "natural")

  # Build monthly GDP dataframe
  gdp_monthly <- data.frame(
    Year = rep(min_year:max_year, each = 12),
    Month = rep(1:12, times = max_year - min_year + 1)
  )
  gdp_monthly <- gdp_monthly[1:length(spline_fit$y), ]
  gdp_monthly$GDP <- spline_fit$y

  cat("[OK] GDP interpolated to monthly frequency\n")
  cat("  Monthly GDP observations:", nrow(gdp_monthly), "\n\n")

  return(gdp_monthly)
}

load_cpi_data <- function() {
  cat("\n+================================================================+\n")
  cat("|  STEP 2B: LOAD CPI / INFLATION DATA (Monthly)                |\n")
  cat("+================================================================+\n\n")

  cat("CPI data is monthly. Note: CPI data may end earlier than\n")
  cat("premium data. Data will be trimmed to the common period.\n\n")

  cpi_raw <- load_file("Select your CPI / Inflation data file (monthly)...")

  cat("Raw CPI data preview:\n")
  print(head(cpi_raw))
  cat("\n")

  return(cpi_raw)
}

prepare_cpi_data <- function(cpi_raw) {
  col_names <- names(cpi_raw)
  col_lower <- tolower(col_names)

  year_col <- col_names[grep("year", col_lower)[1]]
  month_col <- col_names[grep("month", col_lower)[1]]
  cpi_col <- col_names[grep("cpi|inflation|index", col_lower)[1]]

  if (is.na(year_col) || is.na(month_col) || is.na(cpi_col)) {
    cat("[!] Could not auto-detect all columns.\n")
    cat("    Available columns:", paste(col_names, collapse = ", "), "\n")
    cat("    Using first three columns: Year, Month, CPI\n\n")
    year_col <- col_names[1]
    month_col <- col_names[2]
    cpi_col <- col_names[3]
  } else {
    cat("  Detected Year column:", year_col, "\n")
    cat("  Detected Month column:", month_col, "\n")
    cat("  Detected CPI column:", cpi_col, "\n\n")
  }

  cpi_clean <- data.frame(
    Year = as.numeric(cpi_raw[[year_col]]),
    Month = as.numeric(cpi_raw[[month_col]]),
    CPI = as.numeric(cpi_raw[[cpi_col]])
  )
  cpi_clean <- cpi_clean[complete.cases(cpi_clean), ]
  cpi_clean <- cpi_clean[order(cpi_clean$Year, cpi_clean$Month), ]

  cat("[OK] CPI data prepared\n")
  cat("  Observations:", nrow(cpi_clean), "\n")
  cat("  Range:", min(cpi_clean$Year), "M", min(cpi_clean$Month[cpi_clean$Year == min(cpi_clean$Year)]),
      " to ", max(cpi_clean$Year), "M", max(cpi_clean$Month[cpi_clean$Year == max(cpi_clean$Year)]), "\n\n")

  return(cpi_clean)
}

# ============================================================
# STEP 3: MERGE AND ALIGN DATA TO COMMON DATE RANGE
# ============================================================

merge_all_data <- function(residuals_df, gdp_monthly, cpi_data) {
  cat("\n+================================================================+\n")
  cat("|  STEP 3: MERGE & ALIGN DATA TO COMMON DATE RANGE             |\n")
  cat("+================================================================+\n\n")

  cat("Merging datasets by Year and Month...\n\n")

  # Merge residuals with GDP
  merged <- merge(residuals_df, gdp_monthly, by = c("Year", "Month"), all.x = FALSE)
  cat("  After merging with GDP:", nrow(merged), "rows\n")

  # Merge with CPI
  merged <- merge(merged, cpi_data, by = c("Year", "Month"), all.x = FALSE)
  cat("  After merging with CPI:", nrow(merged), "rows\n")

  # Sort by time
  merged <- merged[order(merged$Year, merged$Month), ]

  # Remove any rows with NA
  n_before <- nrow(merged)
  merged <- merged[complete.cases(merged), ]
  n_after <- nrow(merged)

  if (n_before != n_after) {
    cat("  Removed", n_before - n_after, "rows with missing values\n")
  }

  cat("\nALIGNED DATASET:\n")
  cat("----------------------------\n")
  cat("  Total observations:", nrow(merged), "\n")
  cat("  Date range:", min(merged$Year), "M", min(merged$Month[merged$Year == min(merged$Year)]),
      " to ", max(merged$Year), "M", max(merged$Month[merged$Year == max(merged$Year)]), "\n")
  cat("  Features: GDP, CPI\n")
  cat("  Target: ARIMA Residual (eta_t)\n\n")

  cat("DATA SUMMARY:\n")
  cat("----------------------------\n")
  cat("  Residual - Mean:", round(mean(merged$Residual), 2),
      " SD:", round(sd(merged$Residual), 2), "\n")
  cat("  GDP      - Mean:", round(mean(merged$GDP), 2),
      " SD:", round(sd(merged$GDP), 2), "\n")
  cat("  CPI      - Mean:", round(mean(merged$CPI), 2),
      " SD:", round(sd(merged$CPI), 2), "\n\n")

  cat("[!] NOTE: Interest rate data not yet available.\n")
  cat("  When available, add it as an additional feature (xi)\n")
  cat("  to improve the non-linear residual prediction.\n\n")

  return(merged)
}

# ============================================================
# STEP 4: TRAIN/TEST SPLIT
# ============================================================

split_data <- function(merged_data, train_ratio = 0.8) {
  cat("+================================================================+\n")
  cat("|  STEP 4: TRAIN/TEST SPLIT                                    |\n")
  cat("+================================================================+\n\n")

  n <- nrow(merged_data)
  train_size <- floor(train_ratio * n)
  test_size <- n - train_size

  train_data <- merged_data[1:train_size, ]
  test_data <- merged_data[(train_size + 1):n, ]

  cat("Split ratio:", train_ratio * 100, "% train /", (1 - train_ratio) * 100, "% test\n")
  cat("  Training set:", train_size, "observations\n")
  cat("  Test set:    ", test_size, "observations\n\n")

  cat("  Train period:", min(train_data$Year), "M", min(train_data$Month[train_data$Year == min(train_data$Year)]),
      " to ", max(train_data$Year), "M", max(train_data$Month[train_data$Year == max(train_data$Year)]), "\n")
  cat("  Test period: ", min(test_data$Year), "M", min(test_data$Month[test_data$Year == min(test_data$Year)]),
      " to ", max(test_data$Year), "M", max(test_data$Month[test_data$Year == max(test_data$Year)]), "\n\n")

  return(list(train = train_data, test = test_data))
}

# ============================================================
# STEP 5: XGBOOST MODEL TRAINING
# Regularized Objective Function (Eq 3.5):
#   L(phi) = sum(Loss(yi, yhat_i)) + sum(Omega(fk))
# Regularization (Eq 3.6):
#   Omega(fk) = gamma * T + 0.5 * lambda * sum(wk^2)
# ============================================================

train_xgboost <- function(train_data, test_data) {
  cat("\n+================================================================+\n")
  cat("|  STEP 5: XGBOOST MODEL TRAINING                              |\n")
  cat("+================================================================+\n\n")

  cat("XGBOOST PREDICTION EQUATION (Eq 3.4):\n")
  cat("================================================================\n")
  cat("  eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)\n\n")
  cat("  Where:\n")
  cat("    K       = total number of regression trees\n")
  cat("    fk(xi)  = score contributed by the k-th tree\n")
  cat("    xi      = macroeconomic variables (GDP, CPI)\n\n")

  cat("REGULARIZED OBJECTIVE FUNCTION (Eq 3.5):\n")
  cat("================================================================\n")
  cat("  Objective = sum(Loss(yi, yhat_i)) + sum(Omega(fk))\n\n")
  cat("  Loss: Squared error between predicted (eta_hat_t)\n")
  cat("        and actual residuals (eta_t)\n\n")

  cat("REGULARIZATION TERM (Eq 3.6):\n")
  cat("================================================================\n")
  cat("  Omega(fk) = gamma * T + 0.5 * lambda * sum(wk^2)\n\n")
  cat("  Where:\n")
  cat("    T      = number of leaves in the tree\n")
  cat("    wk     = output scores of the leaves\n")
  cat("    gamma  = min loss reduction to split a node\n")
  cat("    lambda = L2 regularization (penalizes large weights)\n\n")

  # Prepare feature matrices
  feature_cols <- c("GDP", "CPI")

  train_features <- as.matrix(train_data[, feature_cols])
  train_target <- train_data$Residual

  test_features <- as.matrix(test_data[, feature_cols])
  test_target <- test_data$Residual

  # Create DMatrix objects for xgboost
  dtrain <- xgb.DMatrix(data = train_features, label = train_target)
  dtest <- xgb.DMatrix(data = test_features, label = test_target)

  # XGBoost hyperparameters with regularization
  params <- list(
    objective = "reg:squarederror",   # Loss function: squared error
    eval_metric = "rmse",
    eta = 0.1,                        # Learning rate
    max_depth = 4,                    # Max tree depth
    min_child_weight = 3,             # Min sum of instance weight in child
    subsample = 0.8,                  # Row subsampling
    colsample_bytree = 1.0,          # Column subsampling (use all with 2 features)
    gamma = 0.1,                      # Eq 3.6: min loss reduction for split
    lambda = 1.0,                     # Eq 3.6: L2 regularization on weights
    alpha = 0.0                       # L1 regularization
  )

  cat("HYPERPARAMETERS:\n")
  cat("----------------------------\n")
  cat("  Learning rate (eta):    ", params$eta, "\n")
  cat("  Max tree depth:         ", params$max_depth, "\n")
  cat("  Min child weight:       ", params$min_child_weight, "\n")
  cat("  Subsample ratio:        ", params$subsample, "\n")
  cat("  Column sample by tree:  ", params$colsample_bytree, "\n")
  cat("  Gamma (min split loss): ", params$gamma, " (Eq 3.6)\n")
  cat("  Lambda (L2 penalty):    ", params$lambda, " (Eq 3.6)\n")
  cat("  Alpha (L1 penalty):     ", params$alpha, "\n")
  cat("  Objective:               reg:squarederror\n\n")

  cat("Training XGBoost with cross-validation...\n\n")

  # Cross-validation to find optimal number of trees (K)
  set.seed(42)
  cv_result <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 500,
    nfold = 5,
    early_stopping_rounds = 30,
    verbose = 0,
    print_every_n = 50
  )

  best_nrounds <- cv_result$best_iteration
  best_rmse <- cv_result$evaluation_log$test_rmse_mean[best_nrounds]

  cat("CROSS-VALIDATION RESULTS:\n")
  cat("----------------------------\n")
  cat("  Optimal K (number of trees):", best_nrounds, "\n")
  cat("  Best CV RMSE:", round(best_rmse, 4), "\n\n")

  # Train final model with optimal K
  cat("Training final model with K =", best_nrounds, "trees...\n\n")

  watchlist <- list(train = dtrain, test = dtest)

  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    watchlist = watchlist,
    verbose = 0
  )

  cat("[OK] XGBoost model trained successfully!\n")
  cat("  Total trees (K):", best_nrounds, "\n")
  cat("  Features used:", paste(feature_cols, collapse = ", "), "\n\n")

  return(list(
    model = xgb_model,
    dtrain = dtrain,
    dtest = dtest,
    train_target = train_target,
    test_target = test_target,
    feature_cols = feature_cols,
    best_nrounds = best_nrounds
  ))
}

# ============================================================
# STEP 6: FEATURE IMPORTANCE
# ============================================================

show_feature_importance <- function(xgb_result) {
  cat("\n+================================================================+\n")
  cat("|  STEP 6: FEATURE IMPORTANCE                                   |\n")
  cat("+================================================================+\n\n")

  cat("Feature importance shows how much each macroeconomic\n")
  cat("variable (xi) contributes to predicting the residuals.\n\n")

  importance <- xgb.importance(
    feature_names = xgb_result$feature_cols,
    model = xgb_result$model
  )

  cat("FEATURE IMPORTANCE (by Gain):\n")
  cat("================================================================\n")
  cat(sprintf("  %-15s | %-10s | %-10s | %-10s\n",
              "Feature", "Gain", "Cover", "Frequency"))
  cat("  ----------------+------------+------------+-----------\n")

  for (i in 1:nrow(importance)) {
    cat(sprintf("  %-15s | %10.4f | %10.4f | %10.4f\n",
                importance$Feature[i],
                importance$Gain[i],
                importance$Cover[i],
                importance$Frequency[i]))
  }
  cat("\n")

  cat("INTERPRETATION:\n")
  cat("  Gain:      Average improvement in loss when feature is used\n")
  cat("  Cover:     Average number of samples affected\n")
  cat("  Frequency: How often feature appears in trees\n\n")

  # Plot
  if (.Platform$OS.type == "windows") {
    windows(width = 10, height = 6)
  } else {
    dev.new(width = 10, height = 6)
  }

  barplot(importance$Gain,
          names.arg = importance$Feature,
          main = "XGBoost Feature Importance (Gain)",
          ylab = "Gain", col = "steelblue",
          border = "navy", cex.main = 1.2)

  cat("[OK] Feature importance plot generated!\n\n")

  return(importance)
}

# ============================================================
# STEP 7: PREDICTIONS & EVALUATION
# Eq 3.7: eta_hat_t = sum(eta_hat_t(k), k=1..K)
# ============================================================

evaluate_model <- function(xgb_result, train_data, test_data) {
  cat("\n+================================================================+\n")
  cat("|  STEP 7: PREDICTIONS & MODEL EVALUATION                      |\n")
  cat("+================================================================+\n\n")

  cat("PREDICTION EQUATION (Eq 3.7):\n")
  cat("  eta_hat_t = sum(eta_hat_t(k), k=1..K)\n")
  cat("  Each tree k contributes eta_hat_t(k) to the final prediction.\n\n")

  model <- xgb_result$model

  # Predictions
  train_pred <- predict(model, xgb_result$dtrain)
  test_pred <- predict(model, xgb_result$dtest)

  train_actual <- xgb_result$train_target
  test_actual <- xgb_result$test_target

  # --- Training Metrics ---
  train_rmse <- sqrt(mean((train_actual - train_pred)^2))
  train_mae <- mean(abs(train_actual - train_pred))
  train_ss_res <- sum((train_actual - train_pred)^2)
  train_ss_tot <- sum((train_actual - mean(train_actual))^2)
  train_r2 <- 1 - (train_ss_res / train_ss_tot)

  # --- Test Metrics ---
  test_rmse <- sqrt(mean((test_actual - test_pred)^2))
  test_mae <- mean(abs(test_actual - test_pred))
  test_ss_res <- sum((test_actual - test_pred)^2)
  test_ss_tot <- sum((test_actual - mean(test_actual))^2)
  test_r2 <- 1 - (test_ss_res / test_ss_tot)

  cat("TRAINING SET PERFORMANCE:\n")
  cat("----------------------------\n")
  cat("  RMSE:", round(train_rmse, 4), "\n")
  cat("  MAE: ", round(train_mae, 4), "\n")
  cat("  R-squared:", round(train_r2, 4), "\n\n")

  cat("TEST SET PERFORMANCE:\n")
  cat("----------------------------\n")
  cat("  RMSE:", round(test_rmse, 4), "\n")
  cat("  MAE: ", round(test_mae, 4), "\n")
  cat("  R-squared:", round(test_r2, 4), "\n\n")

  # Check for overfitting
  rmse_ratio <- train_rmse / test_rmse
  cat("OVERFITTING CHECK:\n")
  cat("----------------------------\n")
  cat("  Train/Test RMSE ratio:", round(rmse_ratio, 4), "\n")
  if (rmse_ratio > 0.5 && rmse_ratio < 1.5) {
    cat("  [OK] Model appears well-generalized\n\n")
  } else if (rmse_ratio < 0.5) {
    cat("  [!] Possible overfitting: train RMSE much lower than test\n\n")
  } else {
    cat("  [!] Unusual: test performs better than train\n\n")
  }

  # Preview predictions
  cat("FIRST 15 TEST PREDICTIONS:\n")
  cat("================================================================\n")
  cat(" Idx | Year | Month | Actual Resid | Predicted Resid | Error\n")
  cat("-----+------+-------+--------------+-----------------+----------\n")

  n_show <- min(15, nrow(test_data))
  for (i in 1:n_show) {
    cat(sprintf(" %3d | %4d |   %2d  | %12.2f | %15.2f | %8.2f\n",
                i, test_data$Year[i], test_data$Month[i],
                test_actual[i], test_pred[i],
                test_actual[i] - test_pred[i]))
  }
  cat("\n")

  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    train_rmse = train_rmse,
    test_rmse = test_rmse,
    train_r2 = train_r2,
    test_r2 = test_r2,
    train_mae = train_mae,
    test_mae = test_mae
  ))
}

# ============================================================
# STEP 8: DIAGNOSTIC PLOTS
# ============================================================

plot_diagnostics <- function(xgb_result, eval_result, train_data, test_data) {
  cat("\n+================================================================+\n")
  cat("|  STEP 8: DIAGNOSTIC PLOTS                                    |\n")
  cat("+================================================================+\n\n")

  if (.Platform$OS.type == "windows") {
    windows(width = 14, height = 10)
  } else {
    dev.new(width = 14, height = 10)
  }

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  # 1. Actual vs Predicted Residuals (Test Set)
  test_actual <- xgb_result$test_target
  test_pred <- eval_result$test_pred

  plot(test_actual, test_pred,
       pch = 16, col = "steelblue",
       main = "Actual vs Predicted Residuals (Test)",
       xlab = "Actual ARIMA Residual (eta_t)",
       ylab = "XGBoost Predicted (eta_hat_t)",
       cex.main = 1.0)
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  grid()

  # 2. Residuals over time (test set)
  xgb_errors <- test_actual - test_pred
  plot(1:length(xgb_errors), xgb_errors,
       type = "l", col = "darkblue",
       main = "XGBoost Prediction Errors Over Time",
       xlab = "Test Observation", ylab = "Error (Actual - Predicted)",
       cex.main = 1.0)
  abline(h = 0, col = "red", lty = 2, lwd = 2)
  grid()

  # 3. Full timeline: actual vs predicted residuals
  train_actual <- xgb_result$train_target
  train_pred <- eval_result$train_pred

  all_actual <- c(train_actual, test_actual)
  all_pred <- c(train_pred, test_pred)
  n_train <- length(train_actual)
  n_total <- length(all_actual)

  plot(1:n_total, all_actual,
       type = "l", col = "darkblue", lwd = 1.5,
       main = "ARIMA Residuals: Actual vs XGBoost Predicted",
       xlab = "Time Period", ylab = "Residual Value",
       cex.main = 1.0)
  lines(1:n_total, all_pred, col = "red", lwd = 1.5, lty = 2)
  abline(v = n_train + 0.5, col = "gray40", lty = 3, lwd = 2)
  legend("topright",
         legend = c("Actual Residual", "XGBoost Predicted", "Train/Test Split"),
         col = c("darkblue", "red", "gray40"),
         lty = c(1, 2, 3), lwd = c(1.5, 1.5, 2), cex = 0.8)
  grid()

  # 4. Distribution of XGBoost errors
  hist(xgb_errors,
       main = "Distribution of XGBoost Errors",
       xlab = "Prediction Error", col = "skyblue", border = "navy",
       breaks = 20, cex.main = 1.0)

  mtext("XGBoost Residual Model Diagnostics", outer = TRUE, cex = 1.1)
  par(mfrow = c(1, 1))

  cat("[OK] Diagnostic plots generated!\n\n")
}

# ============================================================
# STEP 9: SAVE COMBINED RESULTS FOR HYBRID MODEL
# ============================================================

save_results <- function(merged_data, eval_result, train_data, test_data) {
  cat("\n+================================================================+\n")
  cat("|  STEP 9: SAVE COMBINED RESULTS                               |\n")
  cat("+================================================================+\n\n")

  # Build output dataframe with all data
  train_out <- data.frame(
    Year = train_data$Year,
    Month = train_data$Month,
    Actual_Premium = train_data$Actual_Premium,
    ARIMA_Predicted = train_data$Predicted_Premium,
    ARIMA_Residual = train_data$Residual,
    XGBoost_Predicted_Residual = eval_result$train_pred,
    Hybrid_Predicted = train_data$Predicted_Premium + eval_result$train_pred,
    Set = "Train"
  )

  test_out <- data.frame(
    Year = test_data$Year,
    Month = test_data$Month,
    Actual_Premium = test_data$Actual_Premium,
    ARIMA_Predicted = test_data$Predicted_Premium,
    ARIMA_Residual = test_data$Residual,
    XGBoost_Predicted_Residual = eval_result$test_pred,
    Hybrid_Predicted = test_data$Predicted_Premium + eval_result$test_pred,
    Set = "Test"
  )

  combined <- rbind(train_out, test_out)

  # Hybrid model metrics
  cat("HYBRID MODEL (ARIMA + XGBoost) PERFORMANCE:\n")
  cat("================================================================\n\n")

  # Train hybrid metrics
  train_hybrid_error <- train_out$Actual_Premium - train_out$Hybrid_Predicted
  train_hybrid_rmse <- sqrt(mean(train_hybrid_error^2))
  train_hybrid_mae <- mean(abs(train_hybrid_error))
  train_ss_res <- sum(train_hybrid_error^2)
  train_ss_tot <- sum((train_out$Actual_Premium - mean(train_out$Actual_Premium))^2)
  train_hybrid_r2 <- 1 - (train_ss_res / train_ss_tot)

  # Compare with ARIMA-only on train
  train_arima_error <- train_out$Actual_Premium - train_out$ARIMA_Predicted
  train_arima_rmse <- sqrt(mean(train_arima_error^2))

  cat("Training Set:\n")
  cat("  ARIMA-only RMSE:         ", round(train_arima_rmse, 2), "\n")
  cat("  Hybrid (ARIMA+XGB) RMSE: ", round(train_hybrid_rmse, 2), "\n")
  cat("  Hybrid R-squared:        ", round(train_hybrid_r2, 4), "\n\n")

  # Test hybrid metrics
  test_hybrid_error <- test_out$Actual_Premium - test_out$Hybrid_Predicted
  test_hybrid_rmse <- sqrt(mean(test_hybrid_error^2))
  test_hybrid_mae <- mean(abs(test_hybrid_error))
  test_ss_res <- sum(test_hybrid_error^2)
  test_ss_tot <- sum((test_out$Actual_Premium - mean(test_out$Actual_Premium))^2)
  test_hybrid_r2 <- 1 - (test_ss_res / test_ss_tot)

  # Compare with ARIMA-only on test
  test_arima_error <- test_out$Actual_Premium - test_out$ARIMA_Predicted
  test_arima_rmse <- sqrt(mean(test_arima_error^2))

  cat("Test Set:\n")
  cat("  ARIMA-only RMSE:         ", round(test_arima_rmse, 2), "\n")
  cat("  Hybrid (ARIMA+XGB) RMSE: ", round(test_hybrid_rmse, 2), "\n")
  cat("  Hybrid R-squared:        ", round(test_hybrid_r2, 4), "\n\n")

  improvement <- ((test_arima_rmse - test_hybrid_rmse) / test_arima_rmse) * 100
  cat("IMPROVEMENT: XGBoost reduced test RMSE by", round(improvement, 2), "%\n\n")

  # Save
  write.csv(combined, "objective2_hybrid_results.csv", row.names = FALSE)
  cat("[OK] Results saved to 'objective2_hybrid_results.csv'\n\n")

  return(combined)
}

# ============================================================
# STEP 10: HYBRID MODEL COMPARISON PLOT
# ============================================================

plot_hybrid_comparison <- function(combined) {
  cat("\n+================================================================+\n")
  cat("|  STEP 10: HYBRID MODEL COMPARISON PLOT                        |\n")
  cat("+================================================================+\n\n")

  if (.Platform$OS.type == "windows") {
    windows(width = 14, height = 6)
  } else {
    dev.new(width = 14, height = 6)
  }

  par(mar = c(4, 4, 3, 1))

  n <- nrow(combined)
  n_train <- sum(combined$Set == "Train")

  plot(1:n, combined$Actual_Premium,
       type = "l", lwd = 2, col = "darkblue",
       main = "Premium Revenue: Actual vs ARIMA vs Hybrid (ARIMA + XGBoost)",
       xlab = "Time Period", ylab = "Premium Amount",
       cex.main = 1.1)

  lines(1:n, combined$ARIMA_Predicted,
        col = "orange", lwd = 1.5, lty = 2)

  lines(1:n, combined$Hybrid_Predicted,
        col = "red", lwd = 2, lty = 1)

  abline(v = n_train + 0.5, col = "gray40", lty = 3, lwd = 2)

  legend("topright",
         legend = c("Actual Premium", "ARIMA Only", "Hybrid (ARIMA+XGBoost)", "Train/Test Split"),
         col = c("darkblue", "orange", "red", "gray40"),
         lty = c(1, 2, 1, 3),
         lwd = c(2, 1.5, 2, 2),
         cex = 0.9)
  grid()

  cat("[OK] Hybrid comparison plot generated!\n\n")
}

# ============================================================
# MAIN EXECUTION
# ============================================================

main <- function() {
  tryCatch({
    cat("\n================================================================\n")
    cat("  COMPONENT 2: NON-LINEAR RESIDUAL MODELING USING XGBOOST\n")
    cat("  OBJECTIVE 2: Train XGBoost on ARIMA Residuals\n")
    cat("================================================================\n\n")

    cat("THEORY:\n")
    cat("  The ARIMA residuals (eta_t) contain unexplained non-linear\n")
    cat("  patterns. XGBoost uses macroeconomic variables (xi) to\n")
    cat("  predict these residuals via an ensemble of K trees:\n\n")
    cat("    eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)    (Eq 3.4)\n\n")
    cat("  Available features: GDP (interpolated), CPI (monthly)\n")
    cat("  Note: Interest rate data not yet available\n\n")

    # Step 1: Load ARIMA residuals
    residuals_df <- load_residuals()

    # Step 2A: Load and interpolate GDP
    gdp_raw <- load_gdp_data()
    gdp_monthly <- interpolate_gdp_to_monthly(gdp_raw)

    # Step 2B: Load CPI data
    cpi_raw <- load_cpi_data()
    cpi_data <- prepare_cpi_data(cpi_raw)

    # Step 3: Merge and align to common date range
    merged_data <- merge_all_data(residuals_df, gdp_monthly, cpi_data)

    # Step 4: Train/test split (80/20, time-based)
    splits <- split_data(merged_data, train_ratio = 0.8)
    train_data <- splits$train
    test_data <- splits$test

    # Step 5: Train XGBoost
    xgb_result <- train_xgboost(train_data, test_data)

    # Step 6: Feature importance
    importance <- show_feature_importance(xgb_result)

    # Step 7: Evaluate
    eval_result <- evaluate_model(xgb_result, train_data, test_data)

    # Step 8: Diagnostic plots
    plot_diagnostics(xgb_result, eval_result, train_data, test_data)

    # Step 9: Save combined results
    combined <- save_results(merged_data, eval_result, train_data, test_data)

    # Step 10: Hybrid comparison plot
    plot_hybrid_comparison(combined)

    # FINAL SUMMARY
    cat("\n================================================================\n")
    cat("                    FINAL RESULTS\n")
    cat("================================================================\n\n")

    cat("XGBOOST MODEL:\n")
    cat("  Trees (K):    ", xgb_result$best_nrounds, "\n")
    cat("  Features:      GDP, CPI\n")
    cat("  Train RMSE:   ", round(eval_result$train_rmse, 4), "\n")
    cat("  Test RMSE:    ", round(eval_result$test_rmse, 4), "\n")
    cat("  Test R-squared:", round(eval_result$test_r2, 4), "\n\n")

    cat("================================================================\n")
    cat("[OK] OBJECTIVE 2 COMPLETE\n")
    cat("[OK] ARIMA residuals modeled by XGBoost with GDP & CPI\n")
    cat("[OK] Hybrid model results saved to 'objective2_hybrid_results.csv'\n")
    cat("================================================================\n\n")

    cat("NEXT STEPS:\n")
    cat("  - Add Interest Rate data when available as additional feature\n")
    cat("  - The hybrid forecast = ARIMA prediction + XGBoost residual\n")
    cat("  - Combined model captures both linear and non-linear patterns\n\n")

  }, error = function(e) {
    cat("\n[ERROR]:", e$message, "\n\n")
  })
}

main()
