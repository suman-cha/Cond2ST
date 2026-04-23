# =============================================================================
# R/regression.R -- Regression method factories for CIT tests.
#
# Each factory has signature `f(X, y, ..., seed = NULL)` and returns a
# closure `predict(X_new) -> numeric vector`. The CIT runners pass these
# factories as the `regr.method` / `binary.regr.method` arguments. Naming
# convention:
#
#   <backend>_reg_method         continuous response (regression)
#   <backend>_reg_method_binary  binary response (classification probability)
# =============================================================================

#' Linear-model regression closure (continuous y)
#'
#' Fits \code{lm(y ~ X)} and returns a function that predicts on new
#' covariates.
#'
#' @param X,y Training covariate matrix and continuous response vector.
#' @param seed Ignored (deterministic). Present for interface symmetry.
#' @param ... Ignored.
#' @return A function with signature \code{function(X_new)} returning a
#'   numeric vector of predictions.
lm_reg_method <- function(X, y, seed=NULL, ...) {
    n <- length(y)
    X <- as.matrix(X, nrow = n)
    m <- lm(y ~ X)
    return(
        function(X_new) {
            X_new <- as.matrix(X_new)
            as.numeric(predict(m, list(X = X_new)))
        }
    )
}
#' Random-forest regression closure via \pkg{ranger} (continuous y)
#'
#' Fits a 500-tree ranger forest and returns a prediction closure.
#'
#' @inheritParams lm_reg_method
#' @param mtry,max.depth,min.node.size Standard ranger hyperparameters;
#'   \code{mtry} defaults to \code{sqrt(ncol(X))}.
#' @param seed Optional integer seed; controls both fitting and prediction
#'   so paired (CIT vs C2ST) calls are reproducible.
#' @return Same shape as \code{lm_reg_method}.
ranger_reg_method <- function(X, y, mtry = NULL, max.depth = NULL, min.node.size = 5, seed=NULL, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    n <- length(y)
    X <- as.matrix(X, nrow = n)
    d <- dim(X)[2]
    if (is.null(mtry)) {
        mtry <- sqrt(d)
    }
    W <- cbind(y, X)
    colnames(W) <- c("y", 1:d)
    m <- ranger::ranger(
        data = W, dependent.variable.name = "y",
        num.trees = 500, mtry = mtry, max.depth = max.depth,
        min.node.size = min.node.size,
        seed=seed, num.threads=1,
    )
    pred_func <- function(X_new) {
        X_new <- as.matrix(X_new)
        d <- dim(X_new)[2]
        colnames(X_new) <- as.character(1:d)
        as.numeric(predict(m, X_new, seed=seed)$predictions)
    }
    return(pred_func)
}

#' XGBoost regression closure (continuous y)
#'
#' Fits an \code{xgb.train} model with squared-error loss and returns a
#' prediction closure.
#'
#' @inheritParams lm_reg_method
#' @param max.depth,eta,nrounds Standard XGBoost hyperparameters.
#' @param seed Optional integer seed for reproducible boosting.
#' @return Same shape as \code{lm_reg_method}.
xgboost_reg_method <- function(X, y, max.depth = 6, eta = 0.3, nrounds = 100, seed=NULL, ...) {
    if(!is.null(seed)){
        set.seed(seed)
    }
    X <- as.matrix(X)
    colnames(X) <- paste0("V", 1:ncol(X))
    dtrain <- xgb.DMatrix(data = X, label = y, nthread = 1)

    param <- list(max_depth = max.depth, eta = eta,
                  objective = "reg:squarederror", nthread = 1)

    m <- xgb.train(params = param, data = dtrain, nrounds = nrounds, nthread = 1)

    pred_func <- function(X_new) {
        X_new <- as.matrix(X_new)
        colnames(X_new) <- colnames(X)
        as.numeric(predict(m, X_new, seed=seed))
    }

    attr(pred_func, "model") <- m
    return(pred_func)
}

#' Logistic-regression closure (binary y)
#'
#' Fits \code{glm(y ~ X, family = binomial)} and returns a closure that
#' outputs probability of class 1 on new covariates.
#'
#' @inheritParams lm_reg_method
#' @param y 0/1 binary response vector.
#' @return Closure returning class-1 probabilities.
lm_reg_method_binary <- function(X, y, seed=NULL, ...) {
    n <- length(y)
    X <- as.matrix(X, nrow = n)

    m <- glm(y ~ X, family = binomial, ...)

    return(
        function(X_new) {
            X_new <- as.matrix(X_new)
            as.numeric(predict(m, newdata=data.frame(X=X_new), type="response"))
        }
    )
}

#' XGBoost classification closure (binary y)
#'
#' Fits an \code{xgb.train} model with logistic loss and returns a
#' closure that outputs probability of class 1.
#'
#' @inheritParams xgboost_reg_method
#' @param y 0/1 binary response vector.
#' @return Closure returning class-1 probabilities.
xgboost_reg_method_binary <- function(X, y, max.depth = 6, eta = 0.3, nrounds = 100, seed=NULL, ...) {
    if(!is.null(seed)){
        set.seed(seed)
    }

    X <- as.matrix(X)
    colnames(X) <- paste0("V", 1:ncol(X))

    dtrain <- xgb.DMatrix(data = X, label = y, nthread = 1)
    param <- list(max_depth = max.depth, eta = eta,
                  objective = "binary:logistic", nthread = 1)
    m <- xgb.train(params = param, data = dtrain, nrounds = nrounds, nthread = 1)

    pred_func <- function(X_new) {
        X_new <- as.matrix(X_new)
        colnames(X_new) <- colnames(X)
        as.numeric(predict(m, X_new, seed=seed))
    }

    attr(pred_func, "model") <- m
    return(pred_func)
}
#' Random-forest classification closure via \pkg{ranger} (binary y)
#'
#' Fits a probability-forest (\code{probability = TRUE}) and returns a
#' closure that outputs probability of class 1.
#'
#' @inheritParams ranger_reg_method
#' @param y 0/1 binary response vector.
#' @param num.threads ranger thread count; default 1 to avoid contention
#'   with cdtst's own SOCK cluster.
#' @return Closure returning class-1 probabilities.
ranger_reg_method_binary <- function(X, y, mtry = NULL, max.depth = NULL,
                                     min.node.size = 10, seed = NULL,
                                     num.threads = 1, ...) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n <- length(y)
    y <- factor(y)
    X <- as.matrix(X, nrow = n)
    d <- dim(X)[2]

    if (is.null(mtry)) {
        mtry <- sqrt(d)
    }

    colnames(X) <- paste0("V", 1:d)
    W <- data.frame(y = y, X)

    m <- ranger::ranger(
        data = W,
        dependent.variable.name = "y",
        num.trees = 500,
        mtry = mtry,
        max.depth = max.depth,
        min.node.size = min.node.size,
        probability = TRUE,
        classification = TRUE,
        seed = seed,
        num.threads = num.threads
    )

    pred_func <- function(X_new) {
        X_new <- as.matrix(X_new)
        colnames(X_new) <- paste0("V", 1:d)
        if (!is.null(seed)) set.seed(seed)
        preds <- predict(m, data = as.data.frame(X_new),
                         seed = seed, num.threads = num.threads)$predictions
        as.numeric(preds[, 2])
    }

    return(pred_func)
}
