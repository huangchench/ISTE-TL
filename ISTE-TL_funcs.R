rm(list = ls())

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

density_opt<-function(Xs,Xt,lam){
  mm<-function(theta_0){
    theta_0<-as.matrix(theta_0,ncol=1)
    ff<-mean(exp(Xs%*%theta_0))-mean(Xt%*%theta_0)+lam*norm(theta_0,type = "2")^2
    return(ff)
  }
  return(mm)
}

run_CoxTL_BIC <- function(data_t, data_s, weights = NULL, folds_num = 5,
                          lam_set = c(0.001, 0.005, 0.01, 0.05, 0.1, seq(1, 10, by = 1))) {
  
  Ns <- nrow(data_s)
  Nt <- nrow(data_t)
  size_ratio <- Nt / Ns
  
  if (is.null(weights)) {
    weights <- rep(1, Ns)
  }
  
  cov <- c(covs_XA, paste(A, covs_X, sep = ":"))
  
  bic_matrix <- matrix(NA, nrow = length(lam_set), ncol = folds_num)
  rownames(bic_matrix) <- lam_set
  
  split_set <- sample(1:folds_num, Nt, replace = TRUE)
  
  for (lam_idx in 1:length(lam_set)) {

    lambda <- lam_set[lam_idx]  
    weighted_lambda <- size_ratio * lambda
    
    for (k in 1:folds_num) {

      data_t1 <- data_t[split_set == k, ]  
      data_t2 <- data_t[split_set != k, ]   
      data_ts <- rbind(data_s, data_t2)      
      
      functext <- paste0("cox_train <- rms::cph(",
                         "survival::Surv(Y, status) ~ ", 
                         paste(c(cov, 'strat(R)'), collapse = "+"), ",",
                         "data = data_ts, ",
                         "weights = c(", weighted_lambda, " * weights, rep(1, nrow(data_t2))), ",
                         "x = TRUE, y = TRUE, surv = TRUE)")
      eval(parse(text = functext))
      

      coef_est <- coef(cox_train)

      formula_str <- paste("~", 
                           paste(cov, collapse = " + "), 
                           "- 1")
      design_matrix <- model.matrix(as.formula(formula_str), data = data_t1)
      linear_predictor <- design_matrix %*% coef_est


      surv_obj_val <- Surv(data_t1$Y, data_t1$status)
      
      event_indices <- which(data_t1$status == 1)
      if (length(event_indices) > 0) {
        lp_events <- linear_predictor[event_indices]
        time_events <- data_t1$Y[event_indices]
        
        partial_loglik <- 0
        for (i in 1:length(event_indices)) {
          risk_set <- which(data_t1$Y >= time_events[i])
          numerator <- lp_events[i]
          denominator <- log(sum(exp(linear_predictor[risk_set])))
          partial_loglik <- partial_loglik + (numerator - denominator)
        }
        
        p <- sum(abs(coef_est) > 1e-6)
        
        n_val <- nrow(data_t1)
        
        bic_val <- -2 * partial_loglik + log(n_val) * p
        bic_matrix[lam_idx, k] <- bic_val
      }
    }
  }
  
  mean_bic <- rowMeans(bic_matrix, na.rm = TRUE)
  best_idx <- which.min(mean_bic)
  lam_best <- size_ratio * lam_set[best_idx]
  
  
  functext <- paste0("cox_final <- rms::cph(",
                     "survival::Surv(Y, status) ~ ", 
                     paste(c(cov, 'strat(R)'), collapse = "+"), ",",
                     "data = rbind(data_s, data_t), ",
                     "weights = c(", lam_best, " * weights, rep(1, Nt)), ",
                     "x = TRUE, y = TRUE, surv = TRUE)")
  eval(parse(text = functext))
  
  return(cox_final)
}

estimate_iste <-function(data,coef_vector,n){

  f <- as.formula(paste0("~ 0 + ", main_terms, " + ", int_terms))
  M0 <- model.matrix(f, transform(data, A = 0))
  M1 <- model.matrix(f, transform(data, A = 1))
  
  LP0 <- as.numeric(M0 %*% coef_vector)
  LP1 <- as.numeric(M1 %*% coef_vector)
  
  data$LP0_pred <- LP0
  data$LP1_pred <- LP1
  data$LP_pred <- ifelse(data$A == 1, data$LP1_pred, data$LP0_pred)
  

  time_evaluated <- seq(0, 150, by=0.01)
  basehaz <- gbm::basehaz.gbm(t=data$Y, delta=data$status, 
                              f.x=data$LP_pred, 
                              t.eval=time_evaluated)

  n_time <- length(time_evaluated)
  surv_trt_1 <- matrix(NA, n, n_time)
  surv_trt_0 <- matrix(NA, n, n_time)
  

  for (m in 1:n_time){
    surv_trt_1[,m] <- exp(-basehaz[m] * exp(data$LP1_pred))
    surv_trt_0[,m] <- exp(-basehaz[m] * exp(data$LP0_pred))
  }
  
  pred_mst_trt_1 <- sapply(1:n, function(i){
    time_evaluated[which.min(abs(surv_trt_1[i,] - 0.5))]
  })
  
  pred_mst_trt_0 <- sapply(1:n, function(i){
    time_evaluated[which.min(abs(surv_trt_0[i,] - 0.5))]
  })
  data$pred_mst_trt_1 <- pred_mst_trt_1
  data$pred_mst_trt_0 <- pred_mst_trt_0
  
  return_estimate_ate <- list(data, pred_mst_trt_1,  pred_mst_trt_0)
  
  
  return(return_estimate_ate)
}


run_transfer_analysis <- function(method = c("ISTE-TL", "TransCox", "Target-only", "Source-only", "Stratified-Cox"), 
                                  data_t, data_s, covs_X) {
  
  require(dplyr)
  
  suffix <- switch(method, 
                   "ISTE-TL" = "TL", 
                   "TransCox" = "transcox", 
                   "Target-only" = "t", 
                   "Source-only" = "s", 
                   "Stratified-Cox" = "str")
  
  Nt <<- nrow(data_t)
  Ns <<- nrow(data_s)
  p  <<- length(covs_X)
  
  A <<- "A" 
  cov_all <<- c(covs_X, A, paste0(covs_X, A))
  covs_XA <<- c(covs_X,A)
  cov_name <<- c(covs_X, A, paste(covs_X, A, sep = ":"))
  
  main_terms <<- paste(covs_XA, collapse = " + ")
  int_terms  <<- paste(paste(A, covs_X, sep = ":"), collapse = " + ")
  
  
  # --- 1. ISTE-TL  ---
  if (method == "ISTE-TL") {
    X_s <- as.matrix(data_s[, covs_X])
    X_t <- as.matrix(data_t[, covs_X])
    Xs_new <- cbind(1, X_s)
    Xt_new <- cbind(1, X_t)
    
    w_opt <- nlm(density_opt(Xs_new, Xt_new, 1), rep(0, p + 1), iterlim = 500)
    wei_ts <- exp(Xs_new %*% matrix(w_opt$estimate, ncol = 1))
    
    cox_tl <- run_CoxTL_BIC(data_t, data_s, weights = wei_ts)
    coef_vector <- cox_tl$coefficients
    names(coef_vector) <- cov_name
    res <- estimate_iste(data_t, coef_vector, Nt)
    data_result$pred_TL1 <- res[[2]]; data_result$pred_TL0 <- res[[3]]
  }
  
  # --- 2. TransCox  ---
  else if (method == "TransCox") {
    
    prepare_trans <- function(df) {
      df$time <- df$Y
      df$status <- ifelse(df$status == 1, 2, 1)
      for(x in covs_X) df[[paste0(x, "A")]] <- df[[x]] * df[[A]]
      return(df[, c(cov_all, "time", "status")])
    }
    pData <- prepare_trans(data_t)
    aData <- prepare_trans(data_s)
    
    Cout <- GetAuxSurv(aData, cov = cov_all)
    Pout <- GetPrimaryParam(pData, q = Cout$q, estR = Cout$estR)
    
    LRres <- SelLR_By_BIC(pData, aData, cov_all, "status", lambda1=0.1, lambda2=0.1,
                          learning_rate_vec = c(0.001, 0.002, 0.003, 0.004, 0.005), nsteps_vec = c(100, 200))
    
    BICres <- SelParam_By_BIC(pData, aData, cov_all, "status", 
                              lambda1_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                              lambda2_vec = c(0.1, 0.5, seq(1, 10, by = 0.5)),
                              learning_rate = LRres$best_lr, nsteps = LRres$best_nsteps)
    
    Tres <- runTransCox_one(Pout, l1 = BICres$best_la1, l2 = BICres$best_la2, 
                            learning_rate = LRres$best_lr, nsteps = LRres$best_nsteps, cov = cov_all)
    
    res <- estimate_iste(data_t, Tres$new_beta, Nt)
    data_result$pred_transcox1 <- res[[2]]; data_result$pred_transcox0 <- res[[3]]
  }
  
  # --- 3. Target-only / 4. Source-only / 5. Stratified-Cox ---
  else {
    f <- as.formula(paste0("Surv(Y, status) ~ ", main_terms, " + ", int_terms, 
                           ifelse(method == "Stratified-Cox", " + strata(R)", "")))

    fit_data <- switch(method,
                       "Target-only"    = data_t,
                       "Source-only"    = data_s,
                       "Stratified-Cox" = rbind(data_s, data_t))
    
    cox_fit <- survival::coxph(f, data = fit_data)
    res <- estimate_iste(data_t, coef(cox_fit), Nt)
    
    data_result[[paste0("pred_", suffix, "1")]] <- res[[2]]
    data_result[[paste0("pred_", suffix, "0")]] <- res[[3]]
  }
  
  col1 <- paste0("pred_", suffix, "1")
  col0 <- paste0("pred_", suffix, "0")
  ite_name <- paste0("ite_", suffix)
  
  data_result <- data_result %>%
    mutate(!!ite_name := !!sym(col1) - !!sym(col0))
  
  
  return(data_result)
  
}

