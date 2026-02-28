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

load_file <- function(prompt_msg) {
  cat("\n", prompt_msg, "\n")
  fp <- file.choose()
  ext <- tolower(tools::file_ext(fp))
  if (ext %in% c("xlsx", "xls")) return(read_excel(fp))
  if (ext == "csv") return(read.csv(fp))
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
  cat("  Obs:", nrow(df), " | Years:", min(df$Year), "-", max(df$Year), "\n")
  cat("  Residual mean:", round(mean(df$Residual), 4), " | SD:", round(sd(df$Residual), 2), "\n\n")
  return(df)
}

# ============================================================
# STEP 2: LOAD MACROECONOMIC DATA
# ============================================================
load_and_prepare_gdp <- function() {
  cat("+================================================================+\n")
  cat("|  STEP 2A: LOAD GDP (Yearly -> Monthly via Spline)             |\n")
  cat("+================================================================+\n\n")
  gdp_raw <- load_file("Select your GDP data file (yearly)...")
  cn <- names(gdp_raw); cl <- tolower(cn)
  yc <- cn[grep("year", cl)[1]]; gc <- cn[grep("gdp", cl)[1]]
  if (is.na(yc) || is.na(gc)) { yc <- cn[1]; gc <- cn[2] }
  g <- data.frame(Year = as.numeric(gdp_raw[[yc]]), GDP = as.numeric(gdp_raw[[gc]]))
  g <- g[complete.cases(g), ]; g <- g[order(g$Year), ]
  mn <- min(g$Year); mx <- max(g$Year)
  mp <- seq(mn + 0.5/12, mx + 11.5/12, by = 1/12)
  sp <- spline(g$Year + 0.5, g$GDP, xout = mp, method = "natural")
  gm <- data.frame(Year = rep(mn:mx, each = 12), Month = rep(1:12, mx - mn + 1))
  gm <- gm[1:length(sp$y), ]; gm$GDP <- sp$y
  cat("[OK] GDP interpolated:", nrow(gm), "monthly obs\n\n")
  return(gm)
}

load_and_prepare_cpi <- function() {
  cat("+================================================================+\n")
  cat("|  STEP 2B: LOAD CPI (Monthly)                                 |\n")
  cat("+================================================================+\n\n")
  raw <- load_file("Select your CPI data file (monthly)...")
  cn <- names(raw); cl <- tolower(cn)
  yc <- cn[grep("year", cl)[1]]; mc <- cn[grep("month", cl)[1]]
  cc <- cn[grep("cpi|inflation|index", cl)[1]]
  if (is.na(yc)||is.na(mc)||is.na(cc)) { yc <- cn[1]; mc <- cn[2]; cc <- cn[3] }
  d <- data.frame(Year=as.numeric(raw[[yc]]), Month=as.numeric(raw[[mc]]), CPI=as.numeric(raw[[cc]]))
  d <- d[complete.cases(d), ]; d <- d[order(d$Year, d$Month), ]
  cat("[OK] CPI loaded:", nrow(d), "obs\n\n")
  return(d)
}

# ============================================================
# STEP 3: MERGE & FEATURE ENGINEERING
# The key to hybrid improvement: create features that capture
# non-linear macro influence on premium residuals
# ============================================================
merge_and_engineer <- function(residuals_df, gdp_m, cpi_d) {
  cat("+================================================================+\n")
  cat("|  STEP 3: MERGE & FEATURE ENGINEERING                         |\n")
  cat("+================================================================+\n\n")

  m <- merge(residuals_df, gdp_m, by = c("Year", "Month"))
  m <- merge(m, cpi_d, by = c("Year", "Month"))
  m <- m[order(m$Year, m$Month), ]
  rownames(m) <- NULL
  n <- nrow(m)
  cat("  Merged:", n, "rows\n\n")

  # --- Macro rate features ---
  # GDP growth (MoM %)
  m$GDP_Growth <- c(NA, diff(m$GDP) / head(m$GDP, -1) * 100)

  # CPI inflation (MoM %)
  m$CPI_MoM <- c(NA, diff(m$CPI) / head(m$CPI, -1) * 100)

  # CPI year-over-year inflation (%)
  m$CPI_YoY <- c(rep(NA, 12), (m$CPI[13:n] - m$CPI[1:(n-12)]) / m$CPI[1:(n-12)] * 100)

  # --- Lagged residuals (proper indexing) ---
  m$Lag1 <- c(NA, m$Residual[1:(n-1)])
  m$Lag2 <- c(NA, NA, m$Residual[1:(n-2)])
  m$Lag3 <- c(NA, NA, NA, m$Residual[1:(n-3)])
  m$Lag6 <- c(rep(NA, 6), m$Residual[1:(n-6)])
  m$Lag12 <- c(rep(NA, 12), m$Residual[1:(n-12)])

  # --- Interaction: premium sensitivity to macro ---
  m$Premium_GDP_Ratio <- m$Actual_Premium / m$GDP
  m$Premium_CPI_Ratio <- m$Actual_Premium / m$CPI

  # --- Lagged premium (premium momentum) ---
  m$Prem_Lag1 <- c(NA, m$Actual_Premium[1:(n-1)])
  m$Prem_Change <- c(NA, diff(m$Actual_Premium))

  # --- Rolling residual statistics ---
  m$Resid_RollMean3 <- NA
  m$Resid_RollSD3 <- NA
  for (i in 4:n) {
    m$Resid_RollMean3[i] <- mean(m$Residual[(i-3):(i-1)])
    m$Resid_RollSD3[i] <- sd(m$Residual[(i-3):(i-1)])
  }

  # Month indicator
  m$Month_Num <- m$Month

  # Drop NAs from lagging
  m <- m[complete.cases(m), ]
  rownames(m) <- NULL

  feat <- c("GDP", "CPI", "GDP_Growth", "CPI_MoM", "CPI_YoY",
            "Lag1", "Lag2", "Lag3", "Lag6", "Lag12",
            "Premium_GDP_Ratio", "Premium_CPI_Ratio",
            "Prem_Lag1", "Prem_Change",
            "Resid_RollMean3", "Resid_RollSD3", "Month_Num")

  cat("FEATURES (xi) for XGBoost [", length(feat), "total ]:\n")
  cat("  Macro:       GDP, CPI, GDP_Growth, CPI_MoM, CPI_YoY\n")
  cat("  Lagged eta:  Lag1, Lag2, Lag3, Lag6, Lag12\n")
  cat("  Interaction: Premium_GDP_Ratio, Premium_CPI_Ratio\n")
  cat("  Momentum:    Prem_Lag1, Prem_Change\n")
  cat("  Rolling:     Resid_RollMean3, Resid_RollSD3\n")
  cat("  Seasonal:    Month_Num\n")
  cat("  Usable obs: ", nrow(m), "\n\n")

  return(list(data = m, feat = feat))
}

# ============================================================
# STEP 4: TRAIN/TEST SPLIT
# ============================================================
split_data <- function(d, ratio = 0.8) {
  cat("+================================================================+\n")
  cat("|  STEP 4: TRAIN/TEST SPLIT                                    |\n")
  cat("+================================================================+\n\n")
  n <- nrow(d); k <- floor(ratio * n)
  tr <- d[1:k, ]; te <- d[(k+1):n, ]
  cat("  Train:", k, "| Test:", n-k, "\n\n")
  return(list(train = tr, test = te))
}

# ============================================================
# STEP 5: XGBOOST TRAINING WITH TUNING
# ============================================================
train_xgboost <- function(tr, te, feat) {
  cat("+================================================================+\n")
  cat("|  STEP 5: XGBOOST TRAINING                                    |\n")
  cat("+================================================================+\n\n")
  cat("  eta_hat_t = phi(xi) = sum(fk(xi), k=1..K)  (Eq 3.4)\n")
  cat("  Objective = sum(Loss) + sum(Omega(fk))       (Eq 3.5)\n")
  cat("  Omega(fk) = gamma*T + 0.5*lambda*sum(wk^2)  (Eq 3.6)\n\n")

  dtrain <- xgb.DMatrix(data = as.matrix(tr[, feat]), label = tr$Residual)
  dtest  <- xgb.DMatrix(data = as.matrix(te[, feat]), label = te$Residual)

  # Tuning: try multiple parameter sets, pick best CV
  param_grid <- list(
    list(eta=0.1,  max_depth=3, gamma=0,   lambda=1,   subsample=0.8, colsample_bytree=0.8),
    list(eta=0.05, max_depth=4, gamma=0.1, lambda=1,   subsample=0.8, colsample_bytree=0.7),
    list(eta=0.1,  max_depth=4, gamma=0,   lambda=0.5, subsample=0.9, colsample_bytree=0.8),
    list(eta=0.08, max_depth=5, gamma=0,   lambda=1,   subsample=0.7, colsample_bytree=0.8),
    list(eta=0.1,  max_depth=3, gamma=0.1, lambda=0.5, subsample=0.8, colsample_bytree=1.0)
  )

  best_rmse <- Inf; best_params <- NULL; best_nrounds <- 100

  cat("  Tuning", length(param_grid), "parameter sets...\n")
  set.seed(42)

  for (i in seq_along(param_grid)) {
    p <- param_grid[[i]]
    p$objective <- "reg:squarederror"
    p$eval_metric <- "rmse"
    p$min_child_weight <- 3
    p$alpha <- 0

    cv <- suppressWarnings(xgb.cv(
      params = p, data = dtrain, nrounds = 500,
      nfold = 5, early_stopping_rounds = 30, verbose = 0
    ))
    rmse_val <- min(cv$evaluation_log$test_rmse_mean)
    cat(sprintf("    Set %d: eta=%.2f depth=%d gamma=%.1f lambda=%.1f -> CV RMSE=%.2f (K=%d)\n",
                i, p$eta, p$max_depth, p$gamma, p$lambda, rmse_val, cv$best_iteration))
    if (rmse_val < best_rmse) {
      best_rmse <- rmse_val; best_params <- p; best_nrounds <- cv$best_iteration
    }
  }

  cat("\n  Best CV RMSE:", round(best_rmse, 4), "with K =", best_nrounds, "trees\n")
  cat("  Best params: eta=", best_params$eta, " depth=", best_params$max_depth,
      " gamma=", best_params$gamma, " lambda=", best_params$lambda, "\n\n")

  # Train final model with best params
  model <- xgb.train(params = best_params, data = dtrain, nrounds = best_nrounds,
                      watchlist = list(train = dtrain, test = dtest), verbose = 0)

  cat("[OK] XGBoost trained with K =", best_nrounds, "trees\n\n")
  return(list(model=model, dtrain=dtrain, dtest=dtest,
              train_y=tr$Residual, test_y=te$Residual, feat=feat, K=best_nrounds))
}

# ============================================================
# STEP 6: EVALUATION & ARIMA vs HYBRID COMPARISON
# ============================================================
evaluate <- function(xr, tr, te) {
  cat("+================================================================+\n")
  cat("|  STEP 6: ARIMA vs HYBRID COMPARISON                          |\n")
  cat("+================================================================+\n\n")

  tr_pred <- predict(xr$model, xr$dtrain)
  te_pred <- predict(xr$model, xr$dtest)

  # Hybrid = ARIMA + XGBoost residual correction
  tr_hyb <- tr$Predicted_Premium + tr_pred
  te_hyb <- te$Predicted_Premium + te_pred

  # ARIMA-only metrics (test)
  a_rmse <- sqrt(mean((te$Actual_Premium - te$Predicted_Premium)^2))
  a_mae  <- mean(abs(te$Actual_Premium - te$Predicted_Premium))
  a_ss   <- sum((te$Actual_Premium - te$Predicted_Premium)^2)
  a_tot  <- sum((te$Actual_Premium - mean(te$Actual_Premium))^2)
  a_r2   <- 1 - a_ss/a_tot

  # Hybrid metrics (test)
  h_rmse <- sqrt(mean((te$Actual_Premium - te_hyb)^2))
  h_mae  <- mean(abs(te$Actual_Premium - te_hyb))
  h_ss   <- sum((te$Actual_Premium - te_hyb)^2)
  h_r2   <- 1 - h_ss/a_tot

  # ARIMA-only metrics (train)
  a_rmse_tr <- sqrt(mean((tr$Actual_Premium - tr$Predicted_Premium)^2))
  h_rmse_tr <- sqrt(mean((tr$Actual_Premium - tr_hyb)^2))

  imp_rmse <- (a_rmse - h_rmse) / a_rmse * 100
  imp_mae  <- (a_mae - h_mae) / a_mae * 100

  cat("TRAINING SET:\n")
  cat("  ARIMA RMSE:  ", round(a_rmse_tr, 2), "\n")
  cat("  Hybrid RMSE: ", round(h_rmse_tr, 2), "\n\n")

  cat("TEST SET - MODEL COMPARISON:\n")
  cat("================================================================\n")
  cat(sprintf("  %-20s | %12s | %12s | %10s\n", "Metric", "ARIMA", "Hybrid", "Improve"))
  cat("  ---------------------+-------------+-------------+-----------\n")
  cat(sprintf("  %-20s | %12.2f | %12.2f | %9.2f%%\n", "RMSE", a_rmse, h_rmse, imp_rmse))
  cat(sprintf("  %-20s | %12.2f | %12.2f | %9.2f%%\n", "MAE", a_mae, h_mae, imp_mae))
  cat(sprintf("  %-20s | %12.4f | %12.4f |\n", "R-squared", a_r2, h_r2))
  cat("\n")

  if (imp_rmse > 0) {
    cat("[OK] HYBRID IMPROVES ARIMA by", round(imp_rmse, 2), "% RMSE\n\n")
  } else {
    cat("[!] No improvement detected.\n\n")
  }

  return(list(tr_pred=tr_pred, te_pred=te_pred, tr_hyb=tr_hyb, te_hyb=te_hyb,
              a_rmse=a_rmse, h_rmse=h_rmse, a_r2=a_r2, h_r2=h_r2, imp=imp_rmse))
}

# ============================================================
# STEP 7: FEATURE IMPORTANCE
# ============================================================
show_importance <- function(xr) {
  cat("+================================================================+\n")
  cat("|  STEP 7: FEATURE IMPORTANCE                                   |\n")
  cat("+================================================================+\n\n")
  imp <- xgb.importance(feature_names = xr$feat, model = xr$model)
  for (i in 1:nrow(imp))
    cat(sprintf("  %-20s Gain:%.4f Cover:%.4f Freq:%.4f\n",
                imp$Feature[i], imp$Gain[i], imp$Cover[i], imp$Frequency[i]))
  cat("\n")
  if (.Platform$OS.type == "windows") windows(width=10,height=6) else dev.new(width=10,height=6)
  barplot(imp$Gain, names.arg=imp$Feature, main="Feature Importance (Gain)",
          ylab="Gain", col="steelblue", border="navy", las=2, cex.names=0.7)
  cat("[OK] Feature importance plot generated!\n\n")
}

# ============================================================
# STEP 8: PLOTS
# ============================================================
plot_results <- function(xr, ev, tr, te) {
  cat("+================================================================+\n")
  cat("|  STEP 8: DIAGNOSTIC PLOTS                                    |\n")
  cat("+================================================================+\n\n")
  if (.Platform$OS.type == "windows") windows(width=14,height=10) else dev.new(width=14,height=10)
  par(mfrow=c(2,2), mar=c(4,4,3,1), oma=c(0,0,2,0))

  # 1. Scatter: actual vs predicted residuals
  plot(xr$test_y, ev$te_pred, pch=16, col="steelblue",
       main="Actual vs Predicted Residuals (Test)", xlab="Actual", ylab="Predicted")
  abline(0,1,col="red",lwd=2,lty=2); grid()

  # 2. Prediction errors
  er <- xr$test_y - ev$te_pred
  plot(er, type="l", col="darkblue", main="XGBoost Errors (Test)", xlab="Obs", ylab="Error")
  abline(h=0,col="red",lty=2); grid()

  # 3. Actual vs ARIMA vs Hybrid
  nt <- nrow(tr); na <- nt + nrow(te)
  aa <- c(tr$Actual_Premium, te$Actual_Premium)
  ar <- c(tr$Predicted_Premium, te$Predicted_Premium)
  hy <- c(ev$tr_hyb, ev$te_hyb)
  plot(1:na, aa, type="l", lwd=2, col="darkblue",
       main="Actual vs ARIMA vs Hybrid", xlab="Time", ylab="Premium")
  lines(1:na, ar, col="orange", lwd=1.5, lty=2)
  lines(1:na, hy, col="red", lwd=2)
  abline(v=nt+0.5, col="gray40", lty=3, lwd=2)
  legend("topright", c("Actual","ARIMA","Hybrid","Split"),
         col=c("darkblue","orange","red","gray40"), lty=c(1,2,1,3), lwd=2, cex=0.7)
  grid()

  # 4. Error histogram
  hist(er, main="Error Distribution", xlab="Error", col="skyblue", border="navy", breaks=20)
  mtext("XGBoost Residual Model Diagnostics", outer=TRUE, cex=1.1)
  par(mfrow=c(1,1))
  cat("[OK] Plots generated!\n\n")
}

# ============================================================
# STEP 9: SAVE
# ============================================================
save_results <- function(ev, tr, te) {
  cat("+================================================================+\n")
  cat("|  STEP 9: SAVE RESULTS                                        |\n")
  cat("+================================================================+\n\n")
  out <- rbind(
    data.frame(Year=tr$Year, Month=tr$Month, Actual_Premium=tr$Actual_Premium,
               ARIMA_Predicted=tr$Predicted_Premium, ARIMA_Residual=tr$Residual,
               XGBoost_Predicted_Residual=ev$tr_pred, Hybrid_Predicted=ev$tr_hyb, Set="Train"),
    data.frame(Year=te$Year, Month=te$Month, Actual_Premium=te$Actual_Premium,
               ARIMA_Predicted=te$Predicted_Premium, ARIMA_Residual=te$Residual,
               XGBoost_Predicted_Residual=ev$te_pred, Hybrid_Predicted=ev$te_hyb, Set="Test"))
  write.csv(out, "objective2_hybrid_results.csv", row.names=FALSE)
  cat("[OK] Saved to 'objective2_hybrid_results.csv'\n\n")
}

# ============================================================
# MAIN
# ============================================================
main <- function() {
  tryCatch({
    cat("\n================================================================\n")
    cat("  COMPONENT 2: NON-LINEAR RESIDUAL MODELING USING XGBOOST\n")
    cat("================================================================\n\n")

    res <- load_residuals()
    gdp <- load_and_prepare_gdp()
    cpi <- load_and_prepare_cpi()

    fe <- merge_and_engineer(res, gdp, cpi)
    sp <- split_data(fe$data)

    xr <- train_xgboost(sp$train, sp$test, fe$feat)
    ev <- evaluate(xr, sp$train, sp$test)

    show_importance(xr)
    plot_results(xr, ev, sp$train, sp$test)
    save_results(ev, sp$train, sp$test)

    cat("================================================================\n")
    cat("  Hybrid RMSE improvement:", round(ev$imp, 2), "%\n")
    cat("  ARIMA R2:", round(ev$a_r2, 4), " | Hybrid R2:", round(ev$h_r2, 4), "\n")
    cat("[OK] OBJECTIVE 2 COMPLETE\n")
    cat("================================================================\n\n")

  }, error = function(e) cat("\n[ERROR]:", e$message, "\n\n"))
}

main()
