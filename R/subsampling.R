# R/subsampling.R
# Algorithm 1 for sample size adjustment.

#' Algorithm 1: subsampling-based reduction to conditional independence testing
#'
#' Given two independent samples of sizes \eqn{n_1} and \eqn{n_2} drawn from
#' the marginal-confounder distributions, returns subsample sizes
#' \eqn{(\tilde n_1, \tilde n_2)} such that the merged sample is approximately
#' i.i.d. from the conditional distribution required by a CIT, up to error
#' \eqn{\epsilon}. \eqn{k} is computed in closed form from the Chernoff bound.
#'
#' @param x1,x2 Conditioning vectors (rows = observations) for groups 1 and 2.
#'   Only their lengths are used; values are not consumed by the subsampler.
#' @param y1,y2 Response vectors for groups 1 and 2.
#' @param seed Optional integer seed for reproducible Binomial draws.
#' @param epsilon Subsampling slack \eqn{\epsilon \in (0, 1)}; defaults to
#'   \code{1/log(n1 + n2)}. Smaller \eqn{\epsilon} gives a tighter
#'   approximation at the cost of a smaller subsample.
#' @return A named list with integer fields \code{tilde_n1} and \code{tilde_n2}
#'   summing to \eqn{\lfloor k(n_1+n_2) \rfloor}.
apply_alg1 <- function(x1, x2, y1, y2, seed=NULL, epsilon=NULL){
  n1 <- length(y1)
  n2 <- length(y2)
  n <- n1 + n2

  if (is.null(epsilon)){
    epsilon <- 1/log(n)
  }
  n_min <- min(n1, n2)

  k <- 1 - (3 * log(epsilon)) / (2 * n_min) - sqrt((1 - (3 * log(epsilon)) / (2 * n_min))^2 - 1)
  tilde_n <- floor(k * n)
  if (!is.null(seed)){
    set.seed(seed)
  }
  tilde_n1 <- rbinom(1, size = tilde_n, prob = n1 / n)
  tilde_n2 <- tilde_n - tilde_n1

  return(list(tilde_n1=tilde_n1, tilde_n2=tilde_n2))
}

#' Algorithm 1 with refined k* (exact Binomial tail bound)
#'
#' Same interface as \code{apply_alg1} but searches for the largest
#' \eqn{k \in (0, 1)} satisfying the exact Binomial-tail constraint
#' \eqn{P(\mathrm{Bin}(\lceil kn \rceil, n_j/n) > n_j) \le \epsilon} for both
#' \eqn{j = 1, 2}. Yields strictly larger subsamples than the closed-form
#' Chernoff bound used by \code{apply_alg1} at the cost of a binary search.
#'
#' @inheritParams apply_alg1
#' @return Same shape as \code{apply_alg1}.
apply_alg1_refined <- function(x1, x2, y1, y2, seed=NULL, epsilon=NULL){
  n1 <- length(y1)
  n2 <- length(y2)
  n  <- n1 + n2

  if (is.null(epsilon)){
    epsilon <- 1/log(n)
  }

  if (!is.null(seed)){
    set.seed(seed)
  }

  # Binary search for a single group
  find_k_single <- function(nj) {
    k_low  <- 0.501
    k_high <- 0.999
    tol    <- 1e-6
    while (k_high - k_low > tol) {
      k_mid   <- (k_low + k_high) / 2
      tilde_n <- ceiling(k_mid * n)
      tail_p  <- pbinom(nj, size = tilde_n, prob = nj / n, lower.tail = FALSE)
      if (tail_p <= epsilon) k_low <- k_mid
      else                   k_high <- k_mid
    }
    k_low
  }

  # Must satisfy both groups: take the minimum
  k <- min(find_k_single(n1), find_k_single(n2))
  tilde_n <- floor(k * n)

  tilde_n1 <- rbinom(1, size = tilde_n, prob = n1 / n)
  tilde_n2 <- tilde_n - tilde_n1

  return(list(tilde_n1=tilde_n1, tilde_n2=tilde_n2))
}
