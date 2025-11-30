# =====================================================================
# Title: DRT Tests — Statistic Distributions & Type-I/Power by Split Ratio
# Author: (your project)
# Description:
#   For tests {LinearMMD_test, CLF_test, CP_test, debiased_test,
#              bootstrap_MMD_test, BlockMMD_test}, run simulations under
#   Null/Alternative; collect (statistic, p-value, rejection),
#   visualize statistic distributions, and summarize Type-I error / Power.
# Outputs:
#   results/drt_stat_records.csv
#   results/drt_type1_power_summary.csv
#   results/fig_stats_density.png
#   results/fig_stats_box.png
#   results/fig_stats_qq.png
# =====================================================================

rm(list = ls())
suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(pbapply)
    library(tmvtnorm)
})

# ---- Project sources (adjust paths if needed) ---------------------------------
# utils.R: MMDl, MMDb, estimate_r, estimate_marginal_ratio, estimate_joint_ratio, etc.
# CP_FunFiles.R: getFinalStat (CP test)
# all_tests.R: reference tests (LinearMMD_test, CLF_test, CP_test, debiased_test, ...)
source("./experiments/utils.R")
source("./experiments/CP_FunFiles.R")
source("./experiments/all_tests.R")

# ---- Data generators -----------------------------------------------------------
generate_data <- function(n, p, group) {
    mu <- c(1, 1, -1, -1, rep(0, p - 4))
    sigma <- diag(1, p)
    lb <- rep(-.5, p); ub <- rep(.5, p)
    if (group == 1) {
        tmvtnorm::rtmvnorm(n, mean = mu, sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
    } else {
        tmvtnorm::rtmvnorm(n, mean = rep(0, p), sigma = sigma, lower = lb, upper = ub, algorithm = "gibbs")
    }
}

generate_y <- function(x, is_null = TRUE, sigma = 2) {
    n <- nrow(x)
    epsilon <- stats::rt(n, df = sigma)
    f0 <- x %*% c(1, -1, 1, 1, rep(0, ncol(x) - 4))
    mean_shift <- if (is_null) 0 else 0.5
    as.numeric(f0 + epsilon + mean_shift)
}

# ---- Helper: safe outer over closures -----------------------------------------
.make_a_outer <- function(marg_ratio, joint_ratio) {
    a_calc <- function(r0, r1) if (r0 < r1) 1 else 0
    function(point0, point1) {
        cr0 <- joint_ratio(point0) * marg_ratio(point0)
        cr1 <- joint_ratio(point1) * marg_ratio(point1)
        return(outer(cr0, cr1, Vectorize(a_calc)))
    }
}

# ---- Test wrappers returning (statistic, pvalue, rejection) --------------------
LinearMMD_test_stat <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, bandwidth=1,
                                est.method="LL", seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    stopifnot(length(y1) == length(y2))
    total_n <- length(y1)
    n <- ceiling(total_n * prop)
    x11 <- x1[1:n,,drop=FALSE]; x12 <- x1[-(1:n),,drop=FALSE]
    y11 <- y1[1:n];              y12 <- y1[-(1:n)]
    x21 <- x2[1:n,,drop=FALSE];  x22 <- x2[-(1:n),,drop=FALSE]
    y21 <- y2[1:n];              y22 <- y2[-(1:n)]
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    stat <- MMDl(x12, x22, y12, y22, h_x=bandwidth, h_y=bandwidth, r_X=r_X2, seed=seed)  # ~N(0,1)
    pval <- 1 - pnorm(stat)
    rej <- as.integer(pval < alpha)
    list(statistic = stat, pvalue = pval, rejection = rej)
}

CLF_test_stat <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, est.method="LL", seed=NULL) {
    # 분류기 기반 테스트의 Acc_hat 정규근사에 기반한 통계량 복원
    if (!is.null(seed)) set.seed(seed)
    stopifnot(length(y1) == length(y2))
    n1 <- length(y1); n2 <- length(y2)
    n12 <- ceiling(n1 * prop); n22 <- ceiling(n2 * prop)
    n11 <- n1 - n12; n21 <- n2 - n22
    x11 <- x1[1:n11,,drop=FALSE]; x12 <- x1[-(1:n11),,drop=FALSE]
    y11 <- y1[1:n11];              y12 <- y1[-(1:n11)]
    x21 <- x2[1:n21,,drop=FALSE];  x22 <- x2[-(1:n21),,drop=FALSE]
    y21 <- y2[1:n21];              y22 <- y2[-(1:n21)]
    
    # 테스트 분할을 통한 분류기 결정과 비율 보정: estimate_r를 테스트 구간에 적용
    te_ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    g12.est <- 1/te_ratios$g12.est
    g22.est <- 1/te_ratios$g22.est
    J12.est <- te_ratios$v12.est * g12.est
    J22.est <- te_ratios$v22.est * g22.est
    
    pred_V1 <- as.integer((g12.est/(J12.est+g12.est)) <= 0.5)  # class 0 로 판정
    pred_V2 <- as.integer((g22.est/(J22.est+g22.est)) >  0.5)  # class 1 로 판정
    
    acc_ratios <- estimate_r(x12, x12, x22, x22, y12, y12, y22, y22, est.method, seed)
    ratio <- 1/acc_ratios$g22.est
    
    A1_hat <- as.integer(pred_V1 == 0)
    A2_hat <- ratio * as.integer(pred_V2 == 1)
    
    bar_A1 <- mean(A1_hat)
    bar_A2 <- mean(A2_hat)
    sigma1_sq <- var(A1_hat)
    sigma2_sq <- var(A2_hat)
    sig <- sigma1_sq + sigma2_sq
    
    n1_te <- length(y12)
    stat <- if (sig > 0) sqrt(n1_te) * (bar_A1 + bar_A2 - 1) / sqrt(sig) else 0.0
    pval <- 1 - pnorm(stat)
    rej <- as.integer(pval < alpha)
    list(statistic = stat, pvalue = pval, rejection = rej)
}

CP_test_stat <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, est.method="LL", seed=NULL) {
    if (!is.null(seed)) set.seed(seed)
    n1 <- length(y1); n2 <- length(y2)
    n12 <- ceiling(n1 * prop); n22 <- ceiling(n2 * prop)
    n11 <- n1 - n12; n21 <- n2 - n22
    x11 <- x1[1:n11,,drop=FALSE]; x12 <- x1[-(1:n11),,drop=FALSE]
    y11 <- y1[1:n11];              y12 <- y1[-(1:n11)]
    x21 <- x2[1:n21,,drop=FALSE];  x22 <- x2[-(1:n21),,drop=FALSE]
    y21 <- y2[1:n21];              y22 <- y2[-(1:n21)]
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    tmp <- getFinalStat(ratios$g12.est, ratios$g22.est, ratios$v12.est, ratios$v22.est)
    stat <- stats::median(tmp$z.hm)               # 대표 통계량
    pval <- pnorm(stat)
    rej <- as.integer(pval < alpha)
    list(statistic = stat, pvalue = pval, rejection = rej)
}

debiased_test_stat <- function(x1, x2, y1, y2, alpha=0.05, est.method="LL", S=2, seed=NULL) {
    # Chen-Lei (2024) DCP 기반 디바이어스 테스트의 교차적합 구현 (통계량/ p-value/ 기각)
    if (!is.null(seed)) set.seed(seed)
    n1 <- length(y1); n2 <- length(y2)
    d <- if (is.null(dim(x1))) 1 else ncol(x1)
    
    data0 <- as.data.frame(x1); data0$y <- y1; data0$class <- 0
    data1 <- as.data.frame(x2); data1$y <- y2; data1$class <- 1
    colnames(data1) <- colnames(data0)
    
    # 1/2 split for ratio training vs test (as in all_tests.R style)
    split0_ind <- sample(seq_len(nrow(data0)), size = floor(0.5 * nrow(data0)))
    split1_ind <- sample(seq_len(nrow(data1)), size = floor(0.5 * nrow(data1)))
    split0 <- data0[split0_ind, ]; split1 <- data1[split1_ind, ]
    test0  <- data0[-split0_ind, ]; test1  <- data1[-split1_ind, ]
    
    # ratio models on split
    split <- rbind(split0, split1)
    marg_ratio  <- estimate_marginal_ratio(split, nrow(split0), nrow(split1), est.method)
    joint_ratio <- estimate_joint_ratio(split, est.method)
    a_outer <- .make_a_outer(marg_ratio, joint_ratio)
    
    # cross-fit folds on test portion
    ind0 <- sample(seq_len(nrow(test0)))
    ind1 <- sample(seq_len(nrow(test1)))
    fold0 <- split(ind0, ceiling(seq_along(ind0) / (nrow(test0) / S)))
    fold1 <- split(ind1, ceiling(seq_along(ind1) / (nrow(test1) / S)))
    
    store_values <- matrix(NA_real_, nrow = nrow(test0), ncol = nrow(test1))
    
    for (j in seq_len(S)) {
        for (k in seq_len(S)) {
            est_data0 <- test0[fold0[[j]], ]
            est_data1 <- test1[fold1[[k]], ]
            nuisance0 <- test0[-fold0[[j]], ]
            nuisance1 <- test1[-fold1[[k]], ]
            nuisance  <- rbind(nuisance0, nuisance1)
            
            # alpha model using nuisance folds
            alpha_data <- as.data.frame(nuisance0[, 1:d])
            as_mat <- a_outer(nuisance0, nuisance1)
            alpha_data$successes <- rowSums(as_mat)
            alpha_data$failures  <- nrow(nuisance1) - rowSums(as_mat)
            colnames(alpha_data) <- c(paste0("V", 1:d), "successes", "failures")
            alpha_model <- glm(cbind(successes, failures) ~ ., data = alpha_data, family = binomial())
            
            gamma <- estimate_marginal_ratio(nuisance, nrow(nuisance0), nrow(nuisance1), est.method)
            
            # estimate on est folds
            as_est <- a_outer(est_data0, est_data1)
            gamma0 <- gamma(est_data0)
            colnames(est_data0)[1:d] <- paste0("V", 1:d)
            colnames(est_data1)[1:d] <- paste0("V", 1:d)
            alpha0 <- predict(alpha_model, est_data0, type = "response")
            alpha1 <- predict(alpha_model, est_data1, type = "response")
            
            for (u in seq_along(fold0[[j]])) {
                for (v in seq_along(fold1[[k]])) {
                    store_values[fold0[[j]][u], fold1[[k]][v]] <-
                        gamma0[u] * as_est[u, v] + alpha1[v] - gamma0[u] * alpha0[u]
                }
            }
        }
    }
    
    value <- sum(store_values) / (nrow(test0) * nrow(test1))
    
    # variance estimate
    variance0_vec <- rep(0, nrow(test0))
    variance1_vec <- rep(0, nrow(test1))
    for (u in seq_len(nrow(test0)))  variance0_vec[u] <- mean(store_values[u, ]) - 0.5
    for (u in seq_len(nrow(test1)))  variance1_vec[u] <- mean(store_values[, u]) - 0.5
    variance0 <- sum(variance0_vec^2) / (nrow(test0) - 1)
    variance1 <- sum(variance1_vec^2) / (nrow(test1) - 1)
    variance  <- (2 * variance0 + 2 * variance1) / (nrow(test0) + nrow(test1))
    
    stat <- (0.5 - value) / sqrt(variance)
    pval <- 1 - pnorm(stat)
    rej  <- as.integer(stat >= qnorm(1 - alpha))
    list(statistic = stat, pvalue = pval, rejection = rej)
}

BlockMMD_test_stat <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, bandwidth=1,
                               block.power=0.5, est.method="LL", seed=NULL) {
    if (!exists("MMDb", mode = "function")) return(NULL)
    if (!is.null(seed)) set.seed(seed)
    total_n <- length(y1); n <- ceiling(total_n * prop)
    x11 <- x1[1:n,,drop=FALSE]; x12 <- x1[-(1:n),,drop=FALSE]
    y11 <- y1[1:n];              y12 <- y1[-(1:n)]
    x21 <- x2[1:n,,drop=FALSE];  x22 <- x2[-(1:n),,drop=FALSE]
    y21 <- y2[1:n];              y22 <- y2[-(1:n)]
    ratios <- estimate_r(x11, x12, x21, x22, y11, y12, y21, y22, est.method, seed)
    r_X2 <- 1/ratios$g22.est
    n_test <- length(y12)
    B <- max(2, floor(n_test ^ block.power)); if (B >= n_test/2) B <- max(2, floor(n_test/3))
    stat <- MMDb(x12, x22, y12, y22, B, h_x=bandwidth, h_y=bandwidth, r_X=r_X2, seed=seed)  # ~N(0,1)
    pval <- 1 - pnorm(stat)
    rej  <- as.integer(pval < alpha)
    list(statistic = stat, pvalue = pval, rejection = rej, block_size = B)
}

bootstrap_MMD_test_stat <- function(x1, x2, y1, y2, prop=0.5, alpha=0.05, est.method="LL", seed=NULL) {
    if (!exists("bootstrap_MMD_test", mode = "function")) return(NULL)
    if (!is.null(seed)) set.seed(seed)
    # 인터페이스가 구현마다 달라 통계량 복원은 일반화 어렵습니다: p-value/기각만 수집
    rej <- bootstrap_MMD_test(x1, x2, y1, y2, prop=prop, alpha=alpha, est.method=est.method, seed=seed)
    list(statistic = NA_real_, pvalue = NA_real_, rejection = as.integer(rej))
}

# ---- Driver: simulate & collect -----------------------------------------------
simulate_and_collect <- function(n = 2000, d = 10, n_sim = 500, seed0 = 1203,
                                 prop_vals = c(0.1, 0.3, 0.5, 0.7, 0.9),
                                 tests = c("LinearMMD_test","CLF_test","CP_test",
                                           "debiased_test","bootstrap_MMD_test","BlockMMD_test"),
                                 est.method = "LL", alpha = 0.05, block.power = 0.5) {
    dispatch <- list(
        LinearMMD_test     = LinearMMD_test_stat,
        CLF_test           = CLF_test_stat,
        CP_test            = CP_test_stat,
        debiased_test      = function(...) debiased_test_stat(..., alpha = alpha),
        bootstrap_MMD_test = bootstrap_MMD_test_stat,
        BlockMMD_test      = function(...) BlockMMD_test_stat(..., alpha = alpha, block.power = block.power)
    )
    records <- list()
    for (hyp in c("null","alternative")) {
        is_null <- (hyp == "null")
        for (prop in prop_vals) {
            cat(sprintf("==> Hyp=%s | prop=%.1f\n", hyp, prop))
            rec_prop <- pblapply(seq_len(n_sim), function(i) {
                seed <- seed0 + i
                set.seed(seed)
                # data
                x1 <- generate_data(n, d, group = 1)
                y1 <- generate_y(x1, is_null = TRUE)
                set.seed(seed + n_sim)
                x2 <- generate_data(n, d, group = 2)
                y2 <- generate_y(x2, is_null = !is_null)  # Alternative일 때만 shift
                
                row_list <- list()
                for (tn in tests) {
                    FUN <- dispatch[[tn]]
                    if (is.null(FUN)) next
                    out <- try({
                        if (tn %in% c("LinearMMD_test","CLF_test","CP_test"))
                            FUN(x1,x2,y1,y2,prop=prop, est.method=est.method, seed=seed)
                        else if (tn == "BlockMMD_test")
                            FUN(x1,x2,y1,y2,prop=prop, est.method=est.method, seed=seed)
                        else if (tn == "bootstrap_MMD_test")
                            FUN(x1,x2,y1,y2,prop=prop, est.method=est.method, seed=seed)
                        else if (tn == "debiased_test")
                            FUN(x1,x2,y1,y2, est.method=est.method, seed=seed)
                        else NULL
                    }, silent = TRUE)
                    if (inherits(out, "try-error") || is.null(out)) {
                        row_list[[length(row_list)+1]] <- data.table(
                            hypothesis = hyp, prop = prop, test = tn,
                            statistic = NA_real_, pvalue = NA_real_, rejection = NA_integer_,
                            block_size = NA_integer_, sim = i
                        )
                    } else {
                        block_size <- if (!is.null(out$block_size)) out$block_size else NA_integer_
                        row_list[[length(row_list)+1]] <- data.table(
                            hypothesis = hyp, prop = prop, test = tn,
                            statistic = out$statistic, pvalue = out$pvalue, rejection = out$rejection,
                            block_size = block_size, sim = i
                        )
                    }
                }
                data.table::rbindlist(row_list)
            })
            records[[length(records)+1]] <- data.table::rbindlist(rec_prop)
        }
    }
    data.table::rbindlist(records, fill = TRUE)
}

# ---- Run ----------------------------------------------------------------------
dir.create("results", showWarnings = FALSE, recursive = TRUE)

res <- simulate_and_collect(
    n = 2000, d = 10, n_sim = 500, seed0 = 1203,
    prop_vals = c(0.1, 0.3, 0.5, 0.7, 0.9),
    tests = c("LinearMMD_test","CLF_test","CP_test","debiased_test","bootstrap_MMD_test","BlockMMD_test"),
    est.method = "LL", alpha = 0.05, block.power = 0.5
)

fwrite(res, "results/drt_stat_records.csv")
cat("Saved: results/drt_stat_records.csv\n")

# ---- Summaries: Type-I error (Null) & Power (Alt) -----------------------------
summ_type1 <- res[hypothesis=="null", .(
    type1_error = mean(rejection==1, na.rm=TRUE),
    valid = sum(!is.na(rejection)),
    total = .N
), by = .(test, prop)]

summ_power <- res[hypothesis=="alternative", .(
    power = mean(rejection==1, na.rm=TRUE),
    valid = sum(!is.na(rejection)),
    total = .N
), by = .(test, prop)]

summ <- merge(summ_type1, summ_power, by = c("test","prop"), all = TRUE)
fwrite(summ, "results/drt_type1_power_summary.csv")
cat("Saved: results/drt_type1_power_summary.csv\n")

# ---- Visualization: statistic distributions (per-test files) ------------------
# 입력: res (simulate_and_collect의 반환값)
# 출력: results/plots_by_test/<TEST>_{density,box,qq}.png

library(ggplot2)
library(data.table)

make_per_test_plots <- function(res, outdir = "results/plots_by_test") {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    # 통계량이 존재하는 것만 대상으로
    viz_all <- res[is.finite(statistic)]
    if (nrow(viz_all) == 0L) {
        message("No finite statistics to plot. Did your wrappers return 'statistic'?")
        return(invisible(NULL))
    }
    
    tests <- sort(unique(viz_all$test))
    for (tn in tests) {
        sub <- viz_all[test == tn]
        if (nrow(sub) == 0L) {
            message(sprintf("Skip %s (no finite statistics)", tn))
            next
        }
        tn_san <- gsub("[^A-Za-z0-9_]+", "_", tn)  # 파일명 안전 처리
        
        # 1) Density (prop으로 facet, Null vs Alt 색/linetype)
        p_den <- ggplot(sub, aes(x = statistic, color = hypothesis, linetype = hypothesis)) +
            geom_density() +
            facet_grid(. ~ factor(prop), labeller = labeller(prop = label_both)) +
            labs(title = sprintf("[%s] Statistic Distribution by Splitting Ratio", tn),
                 x = "Test statistic", y = "Density",
                 color = "Hypothesis", linetype = "Hypothesis") +
            theme_minimal()
        ggsave(file.path(outdir, sprintf("%s_density.png", tn_san)), p_den, width = 10, height = 4.5, dpi = 200)
        
        # 2) Boxplot (y-축은 테스트마다 다를 수 있어 free_y 비추천; 비교를 위해 고정)
        p_box <- ggplot(sub, aes(x = factor(prop), y = statistic, fill = hypothesis)) +
            geom_boxplot(outlier.size = 0.6) +
            labs(title = sprintf("[%s] Statistic Distribution (Boxplot)", tn),
                 x = "Splitting ratio (prop)", y = "Test statistic", fill = "Hypothesis") +
            theme_minimal()
        ggsave(file.path(outdir, sprintf("%s_box.png", tn_san)), p_box, width = 7.5, height = 4.5, dpi = 200)
        
        # 3) Q-Q plot (정규근사 확인; prop으로 facet)
        p_qq <- ggplot(sub, aes(sample = statistic, color = hypothesis)) +
            stat_qq() + stat_qq_line() +
            facet_grid(. ~ factor(prop), labeller = labeller(prop = label_both)) +
            labs(title = sprintf("[%s] Q-Q plots vs Normal(0,1) by Splitting Ratio", tn),
                 x = "Theoretical Quantiles", y = "Sample Quantiles", color = "Hypothesis") +
            theme_minimal()
        ggsave(file.path(outdir, sprintf("%s_qq.png", tn_san)), p_qq, width = 10, height = 4.5, dpi = 200)
        
        message(sprintf("Saved plots for %s → %s", tn, outdir))
    }
    
    invisible(NULL)
}

# ---- Run per-test plotting on your simulation output --------------------------
# res: simulate_and_collect(...) 결과 테이블이어야 합니다.
make_per_test_plots(res)

# ---- (예시) 첫 번째 테스트만 따로 샘플 그림 저장 --------------------------------
#   - 'tests' 중 첫 번째 테스트만 골라 소용량 샘플로 미리보기 그림을 저장합니다.
save_example_plots <- function(res, outdir = "results/plots_by_test_example", max_per_group = 200) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    viz_all <- res[is.finite(statistic)]
    if (nrow(viz_all) == 0L) {
        message("No finite statistics to plot for example.")
        return(invisible(NULL))
    }
    first_test <- sort(unique(viz_all$test))[1]
    sub <- viz_all[test == first_test]
    # 예시가 너무 크면 축약
    set.seed(1L)
    sub <- sub[, .SD[sample(.N, min(.N, max_per_group))], by = .(hypothesis, prop)]
    tn_san <- gsub("[^A-Za-z0-9_]+", "_", first_test)
    
    p_den <- ggplot(sub, aes(x = statistic, color = hypothesis, linetype = hypothesis)) +
        geom_density() +
        facet_grid(. ~ factor(prop), labeller = labeller(prop = label_both)) +
        labs(title = sprintf("[Example | %s] Density by Splitting Ratio", first_test),
             x = "Test statistic", y = "Density",
             color = "Hypothesis", linetype = "Hypothesis") +
        theme_minimal()
    ggsave(file.path(outdir, sprintf("%s_example_density.png", tn_san)), p_den, width = 10, height = 4.5, dpi = 200)
    
    p_box <- ggplot(sub, aes(x = factor(prop), y = statistic, fill = hypothesis)) +
        geom_boxplot(outlier.size = 0.6) +
        labs(title = sprintf("[Example | %s] Boxplot", first_test),
             x = "Splitting ratio (prop)", y = "Test statistic", fill = "Hypothesis") +
        theme_minimal()
    ggsave(file.path(outdir, sprintf("%s_example_box.png", tn_san)), p_box, width = 7.5, height = 4.5, dpi = 200)
    
    p_qq <- ggplot(sub, aes(sample = statistic, color = hypothesis)) +
        stat_qq() + stat_qq_line() +
        facet_grid(. ~ factor(prop), labeller = labeller(prop = label_both)) +
        labs(title = sprintf("[Example | %s] Q-Q vs Normal(0,1)", first_test),
             x = "Theoretical Quantiles", y = "Sample Quantiles", color = "Hypothesis") +
        theme_minimal()
    ggsave(file.path(outdir, sprintf("%s_example_qq.png", tn_san)), p_qq, width = 10, height = 4.5, dpi = 200)
    
    message(sprintf("Saved example plots for '%s' → %s", first_test, outdir))
}