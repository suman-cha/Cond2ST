# R/classifier.R
# Learned weighted-classifier based conditional two-sample tests.
# Replaces the plugin-rule CLF_test/CV_CLF_test (preserved as
# CLF_plugin_test/CV_CLF_plugin_test in tests.R).
#
# Classifier family is exposed via the `clf.method` hyperparameter:
#   "glm", "glmnet", "ranger", "xgboost"  (default: "ranger")
#
# Semantics: train a weighted binary classifier on (x,y) pairs labelled
# {0 for group 1, 1 for group 2}. Group-2 training samples are weighted by
# r2 = P_1(x,y)/P_2(x,y) so the weighted class-2 joint aligns with class-1
# under H_0. On a held-out fold the test stat is a normalised accuracy
#   T = sqrt(n_te) * (mean(A1_hat) + mean(A2_hat) - 1) / sqrt(var(A1) + var(A2))
# which is asymptotically N(0,1) under H_0 (conditional equality P_1 = P_2).
# Larger T rejects. See compare_learned_clf_families_scenarios_1_3.R for the
# reference implementation this file was ported from.

# ── helpers ───────────────────────────────────────────────────────────────────

# Variance with NaN/NA/negative guard. Called on 0/1 vectors so the only real
# edge cases are length<=1 (returns NA) or all-identical (returns 0); the
# `<= 0` guard in `.clf_fold_stat_generic` also protects downstream.
.safe_var <- function(v) {
    s <- stats::var(v)
    if (!is.finite(s) || s < 0) 0 else s
}

# Marginal density ratio P_1(x)/P_2(x) evaluated at B-side points.
# Input convention mirrors `estimate_r`: A-side is the ratio fit, B-side is
# where the ratio is evaluated.
#
# Why MARGINAL, not joint: the classifier-accuracy test statistic is
#   T = sqrt(n) * (bar_A1 + bar_A2 - 1) / sqrt(var(A1) + var(A2))
# with A1 = 1{p <= 0.5} on group 1 and A2 = r2 * 1{p > 0.5} on group 2. If r2
# were the JOINT ratio P_1(x,y)/P_2(x,y), importance sampling would give
# bar_A2 = E_{P_1}[1{p > 0.5}] and bar_A1 + bar_A2 collapses to 1 identically
# (both terms are probabilities under the same P_1 distribution) — the test
# would have no power. The correct weighting for CONDITIONAL two-sample is
# the marginal ratio r2 = P_1(x)/P_2(x), which re-balances the x-marginal to
# P_1 while keeping the y|x conditional at P_2. Under H_0: P_1(y|x) = P_2(y|x)
# the weighted A2 is unbiased for A1's complement. Under H_A the two terms
# separate, yielding power. This also matches the (old) plugin-rule CLF_test
# in R/hypothesis_tests.R which used `ratio <- 1/Acc_ratios$g22.est`.
.estimate_ratio_targets <- function(x1_A, x1_B, x2_A, x2_B,
                                    y1_A, y1_B, y2_A, y2_B,
                                    est.method = "LL", seed = NULL) {
    er <- estimate_r(x1_A, x1_B, x2_A, x2_B,
                     y1_A, y1_B, y2_A, y2_B,
                     est.method = est.method, seed = seed)
    list(r1 = 1 / er$g12.est,
         r2 = 1 / er$g22.est)
}

# ── weighted classifier fit/predict ───────────────────────────────────────────
# Ported verbatim from compare_learned_clf_families_scenarios_1_3.R:34-122.
# `x_mat` is computed once at function entry and reused across families to
# avoid the repeated as.matrix coercion the reference code did inside every
# branch. `nthread = 1` is set for xgboost for determinism under parallel
# cdtst_sapply workers (ranger already pins threads).
.fit_generic_weighted_clf <- function(x1, x2, y1, y2, w2,
                                      clf.method = c("glm", "glmnet",
                                                     "ranger", "xgboost"),
                                      seed = NULL) {
    clf.method <- match.arg(clf.method)
    if (!is.null(seed)) set.seed(seed)

    data_1 <- as.data.frame(x1); data_1$Y <- y1
    data_2 <- as.data.frame(x2); data_2$Y <- y2

    train_x <- rbind(data_1, data_2)
    train_y <- c(rep(0, nrow(data_1)), rep(1, nrow(data_2)))
    train_w <- c(rep(1, nrow(data_1)), w2)

    if (clf.method == "glm") {
        train_df <- train_x
        train_df$Group <- train_y
        model <- glm(Group ~ ., data = train_df,
                     family = binomial(), weights = train_w)
        return(list(method = clf.method, model = model,
                    feature_names = names(train_x)))
    }

    x_mat <- as.matrix(train_x)

    if (clf.method == "glmnet") {
        model <- glmnet::cv.glmnet(
            x = x_mat, y = train_y,
            family = "binomial",
            weights = train_w,
            alpha = 1,
            nfolds = 5,
            type.measure = "deviance"
        )
        return(list(method = clf.method, model = model,
                    feature_names = colnames(x_mat)))
    }

    if (clf.method == "ranger") {
        train_df <- data.frame(Group = factor(train_y), train_x)
        model <- ranger::ranger(
            data = train_df,
            dependent.variable.name = "Group",
            probability = TRUE,
            classification = TRUE,
            case.weights = train_w,
            num.trees = 300,
            max.depth = 6,
            num.threads = 1,
            seed = seed
        )
        return(list(method = clf.method, model = model,
                    feature_names = names(train_x)))
    }

    # xgboost
    dtrain <- xgboost::xgb.DMatrix(data = x_mat, label = train_y,
                                   weight = train_w)
    model <- xgboost::xgb.train(
        params = list(
            objective = "binary:logistic",
            eval_metric = "logloss",
            max_depth = 4,
            eta = 0.1,
            subsample = 0.9,
            colsample_bytree = 0.9,
            nthread = 1
        ),
        data = dtrain,
        nrounds = 120,
        verbose = 0
    )
    list(method = clf.method, model = model,
         feature_names = colnames(x_mat))
}

# Class-2 posterior at (x, y). Only constructs the representation the chosen
# family consumes (data.frame for glm/ranger, matrix for glmnet/xgboost) to
# avoid the double coercion the reference code did.
.predict_generic_weighted_clf <- function(fit_obj, x, y) {
    if (fit_obj$method == "glm") {
        new_df <- as.data.frame(x); new_df$Y <- y
        return(as.numeric(predict(fit_obj$model,
                                  newdata = new_df, type = "response")))
    }

    if (fit_obj$method == "glmnet") {
        new_df <- as.data.frame(x); new_df$Y <- y
        new_mat <- as.matrix(new_df)
        return(as.numeric(predict(fit_obj$model,
                                  newx = new_mat, s = "lambda.min",
                                  type = "response")))
    }

    if (fit_obj$method == "ranger") {
        new_df <- as.data.frame(x); new_df$Y <- y
        preds <- predict(fit_obj$model, data = new_df,
                         num.threads = 1)$predictions
        return(as.numeric(preds[, 2]))
    }

    # xgboost
    new_df <- as.data.frame(x); new_df$Y <- y
    new_mat <- as.matrix(new_df)
    as.numeric(predict(fit_obj$model, newdata = new_mat))
}

# Per-fold normalised test statistic.
.clf_fold_stat_generic <- function(classifier,
                                   x1_test, x2_test,
                                   y1_test, y2_test,
                                   r2_test) {
    prob_1 <- .predict_generic_weighted_clf(classifier, x1_test, y1_test)
    prob_2 <- .predict_generic_weighted_clf(classifier, x2_test, y2_test)

    A1_hat <- as.integer(prob_1 <= 0.5)
    A2_hat <- r2_test * as.integer(prob_2 > 0.5)

    bar_A1 <- mean(A1_hat)
    bar_A2 <- mean(A2_hat)
    sig <- .safe_var(A1_hat) + .safe_var(A2_hat)

    if (sig <= 0) return(0)

    sqrt(length(y1_test)) * (bar_A1 + bar_A2 - 1) / sqrt(sig)
}

# ── public tests ──────────────────────────────────────────────────────────────

#' DRT: learned weighted-classifier two-sample test (single split).
#'
#' Trains a weighted binary classifier (\code{clf.method}) on group-1 vs
#' group-2 samples, with group-2 weights given by an estimated density
#' ratio. Rejects \eqn{H_0} when the held-out classification accuracy
#' significantly exceeds 0.5.
#'
#' @inheritParams debiased_test
#' @param clf.method One of \code{"glm"}, \code{"glmnet"}, \code{"ranger"},
#'   \code{"xgboost"}.
#' @param clf_split Fraction of the test fold used for classifier fitting.
#' @param clf_ratio_prop,est.prop Sample-split fractions for the inner
#'   density-ratio estimation.
#' @return Integer 0 (do not reject) or 1 (reject) at level \code{alpha}.
CLF_test <- function(x1, x2, y1, y2,
                     clf.method = "ranger",
                     clf_split = 0.5, clf_ratio_prop = 0.5,
                     est.prop = 0.8, alpha = 0.05,
                     est.method = "LL", seed = NULL) {
    clf.method <- match.arg(clf.method,
                            c("glm", "glmnet", "ranger", "xgboost"))
    if (!is.null(seed)) set.seed(seed)

    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)

    n_infer <- floor(total_n * clf_split)
    if (n_infer < 2 || n_infer >= total_n) {
        stop("clf_split must leave at least two observations per group in both D_infer and D_clf")
    }

    infer_idx <- seq_len(n_infer)
    clf_idx   <- (n_infer + 1):total_n

    x1_infer <- x1[infer_idx, , drop = FALSE]; y1_infer <- y1[infer_idx]
    x2_infer <- x2[infer_idx, , drop = FALSE]; y2_infer <- y2[infer_idx]
    x1_clf   <- x1[clf_idx, , drop = FALSE];   y1_clf   <- y1[clf_idx]
    x2_clf   <- x2[clf_idx, , drop = FALSE];   y2_clf   <- y2[clf_idx]

    n_clf <- length(y1_clf)
    n_clf_ratio <- ceiling(n_clf * clf_ratio_prop)
    n_clf_ratio <- min(max(n_clf_ratio, 1), n_clf - 1)

    x1_clf_ratio <- x1_clf[seq_len(n_clf_ratio), , drop = FALSE]
    y1_clf_ratio <- y1_clf[seq_len(n_clf_ratio)]
    x2_clf_ratio <- x2_clf[seq_len(n_clf_ratio), , drop = FALSE]
    y2_clf_ratio <- y2_clf[seq_len(n_clf_ratio)]

    x1_clf_train <- x1_clf[-seq_len(n_clf_ratio), , drop = FALSE]
    y1_clf_train <- y1_clf[-seq_len(n_clf_ratio)]
    x2_clf_train <- x2_clf[-seq_len(n_clf_ratio), , drop = FALSE]
    y2_clf_train <- y2_clf[-seq_len(n_clf_ratio)]

    clf_ratio <- .estimate_ratio_targets(
        x1_clf_ratio, x1_clf_train, x2_clf_ratio, x2_clf_train,
        y1_clf_ratio, y1_clf_train, y2_clf_ratio, y2_clf_train,
        est.method = est.method,
        seed = if (is.null(seed)) NULL else seed + 1L
    )

    classifier <- .fit_generic_weighted_clf(
        x1 = x1_clf_train, x2 = x2_clf_train,
        y1 = y1_clf_train, y2 = y2_clf_train,
        w2 = clf_ratio$r2,
        clf.method = clf.method,
        seed = if (is.null(seed)) NULL else seed + 2L
    )

    n_ratio <- ceiling(n_infer * est.prop)
    n_ratio <- min(max(n_ratio, 1), n_infer - 1)

    x1_ratio <- x1_infer[seq_len(n_ratio), , drop = FALSE]
    y1_ratio <- y1_infer[seq_len(n_ratio)]
    x2_ratio <- x2_infer[seq_len(n_ratio), , drop = FALSE]
    y2_ratio <- y2_infer[seq_len(n_ratio)]

    x1_test <- x1_infer[-seq_len(n_ratio), , drop = FALSE]
    y1_test <- y1_infer[-seq_len(n_ratio)]
    x2_test <- x2_infer[-seq_len(n_ratio), , drop = FALSE]
    y2_test <- y2_infer[-seq_len(n_ratio)]

    eval_ratio <- .estimate_ratio_targets(
        x1_ratio, x1_test, x2_ratio, x2_test,
        y1_ratio, y1_test, y2_ratio, y2_test,
        est.method = est.method,
        seed = if (is.null(seed)) NULL else seed + 3L
    )

    acc_hat <- .clf_fold_stat_generic(
        classifier = classifier,
        x1_test = x1_test, x2_test = x2_test,
        y1_test = y1_test, y2_test = y2_test,
        r2_test = eval_ratio$r2
    )

    as.integer((1 - pnorm(acc_hat)) < alpha)
}

#' DRT: K-fold cross-validated learned classifier two-sample test.
#'
#' K-fold cross-validated counterpart to \code{CLF_test}: aggregates a
#' weighted-classifier accuracy statistic across \code{K} folds, then
#' calibrates rejection at level \code{alpha}.
#'
#' @inheritParams CLF_test
#' @param K Number of cross-validation folds (default 2).
#' @return Integer 0 (do not reject) or 1 (reject) at level \code{alpha}.
CV_CLF_test <- function(x1, x2, y1, y2,
                        clf.method = "ranger",
                        clf_split = 0.5, clf_ratio_prop = 0.5,
                        alpha = 0.05, K = 2,
                        est.method = "LL", seed = NULL) {
    clf.method <- match.arg(clf.method,
                            c("glm", "glmnet", "ranger", "xgboost"))
    if (!is.null(seed)) set.seed(seed)

    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)

    n_infer <- floor(total_n * clf_split)
    if (n_infer < K * 2 || n_infer >= total_n) {
        stop("clf_split and K must leave at least two observations per group in each D_infer fold and nonempty D_clf")
    }

    infer_idx <- seq_len(n_infer)
    clf_idx   <- (n_infer + 1):total_n

    x1_infer <- x1[infer_idx, , drop = FALSE]; y1_infer <- y1[infer_idx]
    x2_infer <- x2[infer_idx, , drop = FALSE]; y2_infer <- y2[infer_idx]
    x1_clf   <- x1[clf_idx, , drop = FALSE];   y1_clf   <- y1[clf_idx]
    x2_clf   <- x2[clf_idx, , drop = FALSE];   y2_clf   <- y2[clf_idx]

    n_clf <- length(y1_clf)
    n_clf_ratio <- ceiling(n_clf * clf_ratio_prop)
    n_clf_ratio <- min(max(n_clf_ratio, 1), n_clf - 1)

    x1_clf_ratio <- x1_clf[seq_len(n_clf_ratio), , drop = FALSE]
    y1_clf_ratio <- y1_clf[seq_len(n_clf_ratio)]
    x2_clf_ratio <- x2_clf[seq_len(n_clf_ratio), , drop = FALSE]
    y2_clf_ratio <- y2_clf[seq_len(n_clf_ratio)]

    x1_clf_train <- x1_clf[-seq_len(n_clf_ratio), , drop = FALSE]
    y1_clf_train <- y1_clf[-seq_len(n_clf_ratio)]
    x2_clf_train <- x2_clf[-seq_len(n_clf_ratio), , drop = FALSE]
    y2_clf_train <- y2_clf[-seq_len(n_clf_ratio)]

    clf_ratio <- .estimate_ratio_targets(
        x1_clf_ratio, x1_clf_train, x2_clf_ratio, x2_clf_train,
        y1_clf_ratio, y1_clf_train, y2_clf_ratio, y2_clf_train,
        est.method = est.method,
        seed = if (is.null(seed)) NULL else seed + 11L
    )

    classifier <- .fit_generic_weighted_clf(
        x1 = x1_clf_train, x2 = x2_clf_train,
        y1 = y1_clf_train, y2 = y2_clf_train,
        w2 = clf_ratio$r2,
        clf.method = clf.method,
        seed = if (is.null(seed)) NULL else seed + 12L
    )

    folds <- sample(rep(1:K, length.out = n_infer))
    acc_values <- numeric(K)

    for (j in seq_len(K)) {
        test_idx  <- which(folds == j)
        ratio_idx <- which(folds != j)

        eval_ratio <- .estimate_ratio_targets(
            x1_infer[ratio_idx, , drop = FALSE],
            x1_infer[test_idx,  , drop = FALSE],
            x2_infer[ratio_idx, , drop = FALSE],
            x2_infer[test_idx,  , drop = FALSE],
            y1_infer[ratio_idx], y1_infer[test_idx],
            y2_infer[ratio_idx], y2_infer[test_idx],
            est.method = est.method,
            seed = if (is.null(seed)) NULL else seed + 100L + j
        )

        acc_values[j] <- .clf_fold_stat_generic(
            classifier = classifier,
            x1_test = x1_infer[test_idx, , drop = FALSE],
            x2_test = x2_infer[test_idx, , drop = FALSE],
            y1_test = y1_infer[test_idx],
            y2_test = y2_infer[test_idx],
            r2_test = eval_ratio$r2
        )
    }

    acc_final <- sum(acc_values) / sqrt(K)
    as.integer((1 - pnorm(acc_final)) < alpha)
}
