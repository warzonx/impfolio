#' Interior-Point Portfolio Optimizer (Logarithmic Barrier)
#'
#' Solves a long-only portfolio optimization problem using a primal
#' logarithmic-barrier interior-point method with Newton steps.
#'
#' @param mu Numeric vector of expected returns.
#' @param cov Optional covariance matrix. If provided, the problem becomes
#'   mean-variance (maximize return - lambda * variance). If NULL, pure
#'   maximum-return problem is solved.
#' @param lambda Risk-aversion parameter (only used when cov is not NULL).
#' @param x0 Starting point (must be strictly positive and sum to 1).
#'   Defaults to the uniform portfolio.
#' @param barrier_start Starting value of the barrier parameter µ.
#' @param barrier_final Final value of µ (smaller → closer to the boundary).
#' @param reduction_factor Factor by which µ is reduced (0 < factor < 1).
#' @param max_iter Maximum number of outer iterations.
#' @param tol Convergence tolerance on the Newton step.
#' @param verbose Logical. Print iteration log?
#'
#' @return A list with:
#'   \item{weights}{Optimal portfolio weights}
#'   \item{expected_return}{Portfolio expected return}
#'   \item{variance}{Portfolio variance (if cov was provided)}
#'   \item{iterations}{Number of iterations performed}
#'
#' @examples
#' mu <- c(0.10, 0.05, 0.07)
#' res <- ipm_portfolio(mu)
#' res$weights
#'
#' @export
ipm_portfolio <- function(mu,
                          cov = NULL,
                          lambda = 1,
                          x0 = NULL,
                          barrier_start = 0.1,
                          barrier_final = 1e-8,
                          reduction_factor = 0.2,
                          max_iter = 50,
                          tol = 1e-10,
                          verbose = TRUE) {

  n <- length(mu)
  if (is.null(x0)) {
    x <- rep(1/n, n)
  } else {
    x <- x0
  }

  # Safety checks
  stopifnot(all(x > 0), abs(sum(x) - 1) < 1e-10)
  if (!is.null(cov)) {
    stopifnot(is.matrix(cov), nrow(cov) == n, ncol(cov) == n)
  }

  mu_barrier <- barrier_start
  history <- list()

  if (verbose) {
    cat("Iter | Barrier µ   | Return    | Step norm\n")
    cat("-------------------------------------------\n")
  }

  for (iter in 1:max_iter) {

    # Objective gradient and Hessian of the barrier problem
    if (is.null(cov)) {
      # Pure maximum return: f = -mu'x - µ sum(log x)
      grad <- -mu - mu_barrier / x
      Hess <- diag(mu_barrier / x^2)
    } else {
      # Mean-variance: f = -mu'x + lambda * x'Σx - µ sum(log x)
      grad <- -mu + 2 * lambda * as.numeric(cov %*% x) - mu_barrier / x
      Hess <- 2 * lambda * cov + diag(mu_barrier / x^2)
    }

    # Newton step on the equality-constrained problem (sum x = 1)
    # We solve the KKT system
    ones <- rep(1, n)
    KKT <- rbind(
      cbind(Hess, ones),
      c(ones, 0)
    )
    rhs <- c(-grad, 0)
    sol <- solve(KKT, rhs)
    delta <- sol[1:n]

    # Backtracking line search to stay strictly feasible
    alpha <- 1
    while (any(x + alpha * delta <= 1e-12)) {
      alpha <- alpha * 0.5
    }
    alpha <- min(alpha * 0.99, 1)   # safety factor

    # Update
    x <- x + alpha * delta
    x <- x / sum(x)                # re-normalize for numerical stability

    step_norm <- sqrt(sum(delta^2))

    if (verbose) {
      port_ret <- sum(mu * x)
      cat(sprintf("%4d | %.2e | %.6f | %.2e\n",
                  iter, mu_barrier, port_ret, step_norm))
    }

    # Reduce barrier parameter when Newton step is small
    if (step_norm < tol) {
      mu_barrier <- mu_barrier * reduction_factor
    }

    if (mu_barrier < barrier_final) break
  }

  # Final results
  result <- list(
    weights = x,
    expected_return = sum(mu * x),
    iterations = iter
  )
  if (!is.null(cov)) {
    result$variance <- as.numeric(t(x) %*% cov %*% x)
  }

  return(result)
}
