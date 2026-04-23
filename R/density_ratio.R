# =============================================================================
# R/density_ratio.R -- Density ratio estimation utilities.
#
# Two estimators are provided: linear-logistic ("LL") and kernel-logistic
# regression ("KLR"). The DRT family in R/hypothesis_tests.R consumes
# `estimate_r` for joint and marginal ratios; the debiased test uses the
# split helpers `estimate_marginal_ratio` and `estimate_joint_ratio`.
# =============================================================================

#' Estimate marginal and conditional density ratios on a held-out split
#'
#' Returns LL or KLR estimates of \eqn{g(x) = p_1(x)/p_2(x)} (marginal) and
#' \eqn{v(y \mid x) = p_1(y \mid x)/p_2(y \mid x)} (conditional), evaluated
#' on both groups' hold-out covariates.
#'
#' @param x11,x12,x21,x22 Train (\code{*1}) and hold-out (\code{*2}) covariate
#'   matrices for groups 1 (\code{x1}) and 2 (\code{x2}).
#' @param y11,y12,y21,y22 Corresponding response splits.
#' @param est.method Either \code{"LL"} (linear logistic) or \code{"KLR"}
#'   (kernel logistic regression).
#' @param seed Optional integer seed.
#' @return A named list with marginal ratios (\code{g12.est}, \code{g22.est})
#'   and conditional ratios (\code{v12.est}, \code{v22.est}).
estimate_r <- function(x11, x12, x21, x22, y11, y12, y21, y22, est.method="LL", seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n11 <- length(y11); n12 <- length(y12)
  n21 <- length(y21); n22 <- length(y22)
  label.fit <- factor(c(rep(0,n11), rep(1,n21)))

  if (est.method == "LL"){
    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    fit.joint <- glm(label.fit~., data=as.data.frame(xy.fit), family=binomial())
    x.fit <- rbind(x11,x21)
    fit.marginal <- glm(label.fit~., data=as.data.frame(x.fit), family=binomial())

    new_data <- rbind(x12,x22)
    prob.marginal <- predict(fit.marginal, newdata=as.data.frame(new_data), type="response")

    g12.est <- prob.marginal[1:n12]/(1-prob.marginal[1:n12])*n11/n21
    g22.est <- prob.marginal[(n12+1):(n12+n22)]/(1-prob.marginal[(n12+1):(n12+n22)])*n11/n21

    new_data <- cbind(new_data, c(y12,y22))
    prob.joint <- predict(fit.joint, newdata=as.data.frame(new_data), type="response")

    v12.est <- (1-prob.joint[1:n12])/prob.joint[1:n12]*g12.est
    v22.est <- (1-prob.joint[(n12+1):(n12+n22)])/prob.joint[(n12+1):(n12+n22)]*g22.est

  } else if (est.method == "QL"){
    xy.fit <- poly(cbind(rbind(x11,x21), c(y11,y21)), degree = 2, raw = TRUE)
    fit.joint <- glm(label.fit~., data=as.data.frame(xy.fit), family=binomial())

    x.fit <- poly(rbind(x11,x21), degree = 2, raw = TRUE)
    fit.marginal <- glm(label.fit~., data=as.data.frame(x.fit), family=binomial())

    new_data <- poly(rbind(x12,x22), degree = 2, raw = TRUE)
    prob.marginal <- predict(fit.marginal, newdata=as.data.frame(new_data), type="response")

    g12.est <- prob.marginal[1:n12] / (1 - prob.marginal[1:n12]) * n11 / n21
    g22.est <- prob.marginal[(n12 + 1):(n12 + n22)] / (1 - prob.marginal[(n12 + 1):(n12 + n22)]) * n11 / n21

    new_data <- poly(cbind(rbind(x12, x22), c(y12, y22)), degree = 2, raw = TRUE)
    prob.joint <- predict(fit.joint, newdata = as.data.frame(new_data), type = "response")

    v12.est <- (1 - prob.joint[1:n12]) / prob.joint[1:n12] * g12.est
    v22.est <- (1 - prob.joint[(n12 + 1):(n12 + n22)]) / prob.joint[(n12 + 1):(n12 + n22)] * g22.est

  } else if (est.method == "KLR"){
    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    data.fit <- constructData(xy.fit, label.fit)
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel='rbfdot', sigma=0.005, lambda=0.0005, tol=1e-6, maxiter=500)
    fit.joint <- klrlearner$learn(data.fit, params)

    x.fit <- rbind(x11, x21)
    data.fit <- constructData(x.fit, label.fit)
    fit.marginal <- klrlearner$learn(data.fit, params)

    newdata <- rbind(x12, x22)
    K <- kernelMult(fit.marginal$kernel, newdata, fit.marginal$data, fit.marginal$alpha)
    pi <- 1 / (1 + exp(-as.vector(K)))

    g12.est <- pi[1:n12]/(1-pi[1:n12])*n11/n21
    g22.est <- pi[(n12+1):(n12+n22)]/(1-pi[(n12+1):(n12+n22)])*n11/n21

    newdata <- cbind(rbind(x12, x22), c(y12, y22))
    K <- kernelMult(fit.joint$kernel, newdata, fit.joint$data, fit.joint$alpha)
    pi <- 1 / (1 + exp(-as.vector(K)))

    v12.est <- (1-pi[1:n12])/pi[1:n12]*g12.est
    v22.est <- (1-pi[(n12+1):(n12+n22)])/pi[(n12+1):(n12+n22)]*g22.est

  } else if (est.method == "NN"){
    hidden.layers <- c(10,10)
    learn.rates <- 0.001
    n.epochs <- 500
    x.fit <- rbind(x11, x21)
    newdata1 <- x12
    newdata2 <- x22

    temp <- NNfun(x.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
                  n.epochs = n.epochs, learn.rates = learn.rates)

    g12.est <- temp$prob1.fit/(1-temp$prob1.fit)*n11/n21
    g22.est <- temp$prob2.fit/(1-temp$prob2.fit)*n11/n21

    xy.fit <- cbind(rbind(x11, x21), c(y11, y21))
    newdata1 <- cbind(x12, y12)
    newdata2 <- cbind(x22, y22)

    temp <- NNfun(xy.fit, label.fit, newdata1, newdata2, nnrep = 5, hidden.layers = hidden.layers,
                  n.epochs = n.epochs, learn.rates = learn.rates)

    v12.est <- (1-temp$prob1.fit)/temp$prob1.fit*g12.est
    v22.est <- (1-temp$prob2.fit)/temp$prob2.fit*g22.est

  } else if (est.method == "KLIEP"){
    xy1 <- cbind(x11, y11)
    xy2 <- cbind(x21, y21)
    fit.joint <- densratio(xy2, xy1, method="KLIEP", verbose=FALSE)
    fit.marginal <- densratio(x21, x11, method="KLIEP", verbose=FALSE)

    g12.est <- fit.marginal$compute_density_ratio(x12)
    g22.est <- fit.marginal$compute_density_ratio(x22)
    v12.est <- fit.joint$compute_density_ratio(cbind(x12, y12))*g12.est
    v22.est <- fit.joint$compute_density_ratio(cbind(x22, y22))*g22.est
  }

  list(g12.est = g12.est, g22.est = g22.est, v12.est=v12.est, v22.est=v22.est)
}


# From Chen and Lei (2024)

# Estimates marginal density ratio given probability
# @param eta estimated probability of being in class 1
# @param n0 sample size of class 0
# @param n1 sample size of class 1
marg <- function(eta,n0,n1){
  return((n0/n1)*(eta/(1-eta)))
}

# Function to compute joint density ratio given probability
joint <- function(eta) {
  return((1 - eta) / eta)
}

# Returns a function that computes estimated marginal density ratio at a point
# @param data data to estimate density ratio
# @param n0 number of points from class 0
# @param n1 number of points from class 1
# @param type "ld" for low dimensional "hd" for high dimensional
estimate_marginal_ratio <- function(data, n0, n1, type) {
  if(type == "LL") {
    # Linear Logistic Regression
    model <- glm(class ~. -y, data = data, family = binomial())
    marg_ratio <- function(x) {
      eta <- predict(model, newdata = x, type = "response")
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "QL") {
    # Quadratic Logistic Regression
    data_numeric <- data.frame(lapply(data[, !names(data) %in% "class"], as.numeric))
    data_poly <- data.frame(poly(data_numeric, degree = 2, raw = TRUE))
    data_poly$class <- data$class

    # data_poly <- data.frame(poly(as.matrix(data[, !names(data) %in% "class"]), degree = 2, raw = TRUE))
    # data_poly$class <- data$class
    model <- glm(class ~ ., data = data_poly, family = "binomial")

    marg_ratio <- function(x) {
      x_poly <- data.frame(poly(as.matrix(x[, !names(x) %in% "class"]), degree = 2, raw = TRUE))
      eta <- predict(model, newdata = x_poly, type = "response")
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "KLR") {
    # Kernel Logistic Regression using kernelMult and kernellearner
    xy.fit <- cbind(as.matrix(data[, !names(data) %in% "class"]))
    label.fit <- as.factor(data$class)

    # Construct the KLR learner
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel = 'rbfdot', sigma=0.005, lambda=0.0005,
                   tol = 1e-6, maxiter = 500)

    # Train the model
    fit.marginal <- klrlearner$learn(constructData(xy.fit, label.fit), params)

    marg_ratio <- function(x) {
      # Predict using the learned model
      new_x <- as.matrix(x[, !names(x) %in% "class"])
      K <- kernelMult(fit.marginal$kernel, new_x, fit.marginal$data, fit.marginal$alpha)
      pi <- 1 / (1 + exp(-as.vector(K)))  # predicted probabilities
      return(sapply(pi, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "superlearner") {
    # SuperLearner approach
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100, 200, 500, 1000),
                                              max_depth = c(2, 4, 6),
                                              shrinkage = 0.3,
                                              minobspernode = 1),
                                  detailed_names = FALSE, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y = data$class,
                          X = data[, !names(data) %in% c("class", "y")],
                          SL.library = c("SL.ranger", "SL.lm", "SL.ksvm", xgb_grid$names),
                          family = binomial())

    marg_ratio <- function(x) {
      eta <- predict(model, x[, !names(x) %in% c("class", "y")], onlySL = TRUE)$pred
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  } else if(type == "nn") {
    # Neural Network approach
    rec <- recipe(class ~., data = data[,-5]) %>%
      step_normalize(all_numeric_predictors())
    model <- brulee_mlp(rec, data = data, epochs = 15000, hidden_units = c(10, 10), learn_rate = 0.001, verbose = TRUE)

    marg_ratio <- function(x) {
      eta <- predict(model, x[, 1:4], type = "prob")$.pred_1
      return(sapply(eta, function(x) { return(marg(x, n0, n1)) }))
    }
  }
  return(marg_ratio)
}
# Function to estimate joint density ratio
estimate_joint_ratio <- function(data, type) {
  if(type == "LL") {
    # Linear Logistic Regression
    model <- glm(class ~ ., data = data, family = binomial())
    joint_ratio <- function(point) {
      eta <- predict(model, newdata = point, type = "response")
      return(sapply(eta, joint))
    }
  } else if(type == "QL") {
    # Quadratic Logistic Regression
    data_poly <- data.frame(poly(as.matrix(data[, !names(data) %in% "class"]), degree = 2, raw = TRUE))
    data_poly$class <- data$class
    model <- glm(class ~ ., data = data_poly, family = "binomial")

    joint_ratio <- function(point) {
      point_poly <- data.frame(poly(as.matrix(point[, !names(point) %in% "class"]), degree = 2, raw = TRUE))
      eta <- predict(model, newdata = point_poly, type = "response")
      return(sapply(eta, joint))
    }
  } else if(type == "KLR") {
    # Kernel Logistic Regression using kernelMult and kernellearner
    xy.fit <- cbind(as.matrix(data[, !names(data) %in% "class"]))
    label.fit <- as.factor(data$class)

    # Construct the KLR learner
    klrlearner <- constructKlogRegLearner()
    params <- list(kernel = 'rbfdot', sigma = 0.005, lambda = 0.05 / nrow(data),
                   tol = 1e-6, maxiter = 500)

    # Train the joint model
    fit.joint <- klrlearner$learn(constructData(xy.fit, label.fit), params)

    joint_ratio <- function(point) {
      # Predict using the learned joint model
      new_point <- as.matrix(point[, !names(point) %in% "class"])
      K <- kernelMult(fit.joint$kernel, new_point, fit.joint$data, fit.joint$alpha)
      pi <- 1 / (1 + exp(-as.vector(K)))  # predicted probabilities
      return(sapply(pi, joint))
    }
  } else if(type == "superlearner") {
    # SuperLearner approach
    xgb_grid <- create.SL.xgboost(tune = list(ntrees = c(100, 200, 500),
                                              max_depth = c(2, 6),
                                              shrinkage = 0.3,
                                              minobspernode = 1),
                                  detailed_names = FALSE, env = .GlobalEnv,
                                  name_prefix = "SL.xgb")
    model <- SuperLearner(Y = data$class,
                          X = data[, !names(data) %in% "class"],
                          SL.library = c("SL.ranger", "SL.lm", xgb_grid$names),
                          family = binomial())

    joint_ratio <- function(point) {
      eta <- predict(model, point[, !names(point) %in% "class"], onlySL = TRUE)$pred
      return(sapply(eta, joint))
    }
  }
  return(joint_ratio)
}
