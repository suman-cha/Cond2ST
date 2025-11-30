# Install the required packages
required_pkgs <- c("glmnet",
                   "caret",
                   "ranger",
                   "mgcv",
                   "nnet",
                   "xgboost",
                   "brulee",
                   "recipes",
                   "kernlab",
                   "SuperLearner",
                   "CVST",
                   "densratio")

install_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(new_pkgs)){
    install.packages(new_pkgs, repos="http://cran.us.r-project.org")
  }
}

install_pkgs(required_pkgs)
source("./experiments/CP_FunFiles.R")

# Gaussian kernel 
gaussian.kernel <- function(x, y = NULL, h=1) {
  if (is.null(y)) {
    res <- (exp(-0.5 * (x/h)^2) / (h * sqrt(2 * pi)))
  } else {
    dist <- sum((x - y)^2)
    res <- exp(-dist / (2 * h^2))
  }
  return(res)
}

median.bandwidth <- function(x, y){
    dists <- as.vector(dist(rbind(x, y)))
    bw <- median(dists[dists > 0])
    return(bw)
}

# Linear MMD
MMDl <- function(x12, x22, y12, y22, h_x=1, h_y=1, r_X, seed=NULL){
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # x12, x22, y12, y22 are splitted data 2n -> n
  stopifnot(length(y12) == length(y22))
  n <- length(y12)
  m <- floor(n/2)
  S_hat_values <- numeric(m)

  for(i in 1:m){
    k_zz <- gaussian.kernel(x12[i, ], x12[i+m, ], h_x) * gaussian.kernel(y12[i], y12[i+m], h_y)
    k_ww <- gaussian.kernel(x22[i, ], x22[i+m, ], h_x) * gaussian.kernel(y22[i], y22[i+m], h_y)
    k_wz <- gaussian.kernel(x22[i, ], x12[i+m, ], h_x) * gaussian.kernel(y22[i], y12[i+m], h_y)
    k_zw <- gaussian.kernel(x12[i, ], x22[i+m, ], h_x) * gaussian.kernel(y12[i], y22[i+m], h_y)
  
    S_hat_values[i] <- k_zz + r_X[i]*r_X[i+m]*k_ww - r_X[i]*k_wz - r_X[i+m]*k_zw
  }
  S_bar <- mean(S_hat_values)
  sigma_hat <- sum((S_hat_values - S_bar)^2) / (m - 1)
  MMDl_hat2 <- sqrt(m) * S_bar / sqrt(sigma_hat)
  return(MMDl_hat2)
}

# Block MMD
MMDb <- function(x12, x22, y12, y22, B_size, h_x=1, h_y=1, r_X, seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    n <- length(y12)
    S <- floor(n / B_size)
    
    if (S < 2) stop("Block size too large for given sample size")
    
    MMD_values <- numeric(S)

    for (b in 1:S){
        idx <- ((b-1)*B_size+1):(b*B_size)
        xb1 <- x12[idx,,drop=FALSE]
        xb2 <- x22[idx,,drop=FALSE]
        yb1 <- y12[idx]
        yb2 <- y22[idx]
        rb <- r_X[idx]
        
        sum_k <- 0
        # (i, j) where i < j
        for (i in 1:(B_size-1)) {
            for (j in (i+1):B_size) {
                k_zz <- gaussian.kernel(xb1[i, ], xb1[j, ], h_x) * gaussian.kernel(yb1[i], yb1[j], h_y)
                k_ww <- gaussian.kernel(xb2[i, ], xb2[j, ], h_x) * gaussian.kernel(yb2[i], yb2[j], h_y)
                k_wz <- gaussian.kernel(xb2[i, ], xb1[j, ], h_x) * gaussian.kernel(yb2[i], yb1[j], h_y)
                k_zw <- gaussian.kernel(xb1[i, ], xb2[j, ], h_x) * gaussian.kernel(yb1[i], yb2[j], h_y)
                
                sum_k <- sum_k + k_zz + rb[i]*rb[j]*k_ww - rb[i]*k_wz - rb[j]*k_zw
                
            }
        }
        n_pairs <- B_size*(B_size - 1) / 2
        # n_pairs <- (B_size - 1) / 2
        MMD_values[b] <- sum_k / n_pairs
    }
    
    S_bar <- mean(MMD_values)
    sigma_hat <- sqrt(var(MMD_values))
    stat <- ifelse(sigma_hat >0, sqrt(S)*S_bar/sigma_hat, 0)
    return(stat)
}

bootstrap_MMD <- function(x12, x22, y12, y22, h_x=1, h_y=1, r_X, B, seed=NULL){
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    n <- length(y12)
    stopifnot(length(y12) == length(y22))
    
    H_hat <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n){
        for (j in 1:n) {
            if (i != j) {
                k_zz <- gaussian.kernel(x12[i,,drop=FALSE], x12[j,,drop=FALSE], 
                                        h_x) * gaussian.kernel(y12[i], y12[j], h_y)
                k_ww <- gaussian.kernel(x22[i,,drop=FALSE], x22[j,,drop=FALSE], 
                                        h_x) * gaussian.kernel(y22[i], y22[j], h_y)
                k_wz <- gaussian.kernel(x22[i,,drop=FALSE], x12[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y22[i], y12[j], h_y)
                k_zw <- gaussian.kernel(x12[i,,drop=FALSE], x22[j,,drop=FALSE],
                                        h_x) * gaussian.kernel(y12[i], y22[j], h_y)
                H_hat[i, j] <- k_zz + r_X[i]*r_X[j]*k_ww - r_X[i]*k_wz - r_X[j]*k_zw
            }
        }
    }
    obs_stat <- sum(H_hat) / (n*(n-1))
    bootstrap_stats <- numeric(B)
    
    for (b in 1:B){
        # Generate n i.i.d. Gaussian random variables
        W <- rnorm(n)
        
        # Calculate the wild boostrap statistic
        wild_bootstrap_sum <- sum(outer(W,W) * H_hat)
        bootstrap_stats[b] <- wild_bootstrap_sum / (n*(n-1))
    }
    return(list(obs_stat=obs_stat, bootstrap_stats=bootstrap_stats))
}

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

# Algorithm 1
apply_alg1 <- function(x1, x2, y1, y2, seed=NULL, epsilon=NULL){
  n1 <- length(y1)
  n2 <- length(y2)
  n <- n1 + n2

  if (is.null(epsilon)){
    epsilon <- 1/log(n)
  }

  k <- 1 - (3 * log(epsilon)) / (2 * n1) - sqrt((1 - (3 * log(epsilon)) / (2 * n1))^2 - 1)
  tilde_n <- floor(k * n)
  # print(tilde_n)
  if (!is.null(seed)){
    set.seed(seed)
    # cat("seed is set to ", seed, "\n")
  }
  tilde_n1 <- rbinom(1, size = tilde_n, prob = n1 / (n1 + n2))
  tilde_n2 <- tilde_n - tilde_n1

  # cat("tilde_n1: ", tilde_n1, "tilde_n2: ", tilde_n2, "\n")
  return(list(tilde_n1=tilde_n1, tilde_n2=tilde_n2))
}



