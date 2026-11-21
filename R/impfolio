impfolio <- function(mu) {

  # Constantes fixées
  n <- 3

  Sigma <- matrix(c(
    0.04, 0.01, 0.02,
    0.01, 0.09, 0.03,
    0.02, 0.03, 0.16
  ), nrow = n, byrow = TRUE)

  # Initialisation
  x <- rep(1 / n, n)

  tol <- 1e-8
  max_iter <- 100

  for (k in 1:max_iter) {

    # Gradient de la fonction objectif avec barrière
    grad <- as.vector(Sigma %*% x) - mu / x

    # Hessienne
    H <- Sigma + diag(mu / x^2)

    # Contrainte : somme des poids = 1
    A <- matrix(1, nrow = 1, ncol = n)

    # Système de Newton
    KKT <- rbind(
      cbind(H, t(A)),
      cbind(A, 0)
    )

    rhs <- c(-grad, 0)

    sol <- solve(KKT, rhs)

    dx <- sol[1:n]

    # Critère d'arrêt
    if (max(abs(dx)) < tol) {
      break
    }

    # Recherche du pas assurant x > 0
    alpha <- 1

    while (any(x + alpha * dx <= 0)) {
      alpha <- alpha / 2
    }

    # Mise à jour
    x <- x + alpha * dx
  }

  return(x)
}

# Exemple
impfolio(mu = 0.01)
