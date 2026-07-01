#' Utility: convert to vector format
#'
#' @param x value(s)
#' @param n number of replicates
#' @export

to_v <- function(x, n) {

  if (length(x) == 1L) {
    return(rep(x, n))
  }

  if (length(x) != n) {
    stop("length(x) must be 1 or equal to n")
  }

  x
}

#' Utility: get maximum trophic position
#'
#' @param w Square binary adjacency matrix where \code{w[i, j] = 1}
#'   indicates species i consumes species j.
#' @param occupancy Numeric vector of species occupancies at equilibrium.
#' @param weight Logical. If \code{TRUE}, trophic positions are weighted by
#'   relative occupancies.
#'
#' @author Akira Terui \email{hanabi0111@gmail.com}
#' @export

maxtp <- function(w,
                  occupancy,
                  weight = TRUE) {

  # check input -------------------------------------------------------------

  absfwb <- abs(w)

  if (!all(absfwb == t(absfwb)))
    stop("Input w is invalid (abs(w) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("Input w is invalid (abs(w) must be binary)")

  if (any(dim(w) != length(occupancy)))
    stop("Input w or occupancy has invalid dimension.")

  # trophic position --------------------------------------------------------

  if (any(occupancy > 0)) {
    ## with at least one persisting species

    ## declare objects
    ## - v_id_p, index of persisting species
    ## - m_sub_fwb, subset the w by persisting species,
    ## -- then empty lower triangle
    ## - v_tp, initialize trophic position vector
    v_id_p <- which(occupancy > 0)
    m_sub_fwb <- as.matrix(absfwb[v_id_p, v_id_p])
    m_sub_fwb[lower.tri(m_sub_fwb)] <- 0
    v_tp <- rep(-1, ncol(m_sub_fwb))

    ## - v_o, occupancy of persisting species
    ## - v_id_b, index of basal species
    ## - n_b, number of basal species
    ## - v_n_prey, number of prey
    ## - v_sum_o, summed occupancy of prey
    v_o <- occupancy[v_id_p]
    v_n_prey <- colSums(m_sub_fwb)
    v_id_b <- which(v_n_prey == 0)
    n_b <- sum(v_n_prey == 0)
    v_sum_o <- drop(v_o %*% m_sub_fwb)

    if (any(v_n_prey > 0)) {
      ## v_tp = 1 for basal
      v_tp[v_id_b] <- 1

      ## update v_tp recursively for consumers
      if (weight) {

        ## - weight = T; calculate v_tp based on relative occupancies
        for (i in (n_b + 1):length(v_tp)) {
          v_tp_prime <- drop(v_tp %*% (m_sub_fwb * v_o))
          v_tp[i] <- v_tp_prime[i] / v_sum_o[i] + 1
        }

      } else {

        ## - weight = F; calculate v_tp based on presence absence
        for (i in (n_b + 1):length(v_tp)) {
          v_tp_prime <- drop(v_tp %*% m_sub_fwb)
          v_tp[i] <-  v_tp_prime[i] / v_n_prey[i] + 1
        }

      }

    } else {

      v_tp[v_id_b] <- 1

    }

    ## - FCL as the maximum trophic position of persisting species
    fcl <- max(v_tp)
    attr(fcl, "tp") <- v_tp

  } else {

    ## all extinct
    fcl <- 0

  }

  return(fcl)
}

#' Utility: Poisson distribution conditioned on even outcomes
#'
#' @param lambda Numeric. Branching rate.
#' @param size Numeric. Total network length.
#' @param min_z Integer. Minimum value of realization.
#'
#' @return Named probability vector for even z = 0, 2, 4, ...
#' @export

cpois <- function(lambda, size, min_z = 0) {

  ## effective Poisson mean
  mu <- lambda * size

  ## handle degenerate case
  if (mu <= 0) {
    pz <- numeric(1)
    names(pz) <- "0"
    return(pz)
  }

  ## truncation for numerical support
  pois_max <- stats::qpois(1 - 1e-10, lambda = mu)

  v_z <- min_z:pois_max
  pr_z <- stats::dpois(v_z, lambda = mu)

  ## keep only even states (conditioning step)
  keep <- (v_z %% 2 == 0)

  pz <- pr_z[keep]
  pz <- pz / sum(pz)
  z <- v_z[keep]

  return(cbind(z, pz))
}

#' Utility: Laplace transform of a Rayleigh distribution
#'
#' Computes \code{E[exp(-delta * d)]} where d follows a Rayleigh distribution
#' with mean distance `mu`.
#'
#' Uses the numerically stable representation involving the scaled
#' complementary error function (`erfcx`).
#'
#' @param delta Non-negative rate parameter.
#' @param mu Mean of the Rayleigh distribution.
#'
#' @return A numeric vector equal to \code{E[exp(-delta * d)]}.
#'
#' @export

laplace_rayleigh <- function(delta, mu) {
  1 - delta * mu * pracma::erfcx(delta * mu / sqrt(pi))
}

#' Utility: Laplace transform of downstream distance
#'
#' Computes the Laplace transform of the downstream distance from a randomly
#' located disturbance to an affected habitat, assuming root-to-leaf distances
#' follow a Rayleigh distribution.
#'
#' @param nu Positive rate parameter of the Laplace transform.
#' @param mu Mean of the Rayleigh distribution for root-to-leaf distance.
#' @param exact Logical.
#'   If \code{TRUE}, computes the transform by numerical integration.
#'   If \code{FALSE}, uses the asymptotic approximation for large
#'   \code{nu * mu}.
#'
#' @return A numeric vector giving \eqn{E[\exp(-\nu d)]}, where
#'   \eqn{d} is the downstream distance from the disturbance to a randomly
#'   selected affected habitat.
#'
#' @author Akira Terui
#'
#' @export

laplace_rt <- function(nu, mu, exact = TRUE) {

  stopifnot(nu >= 0, mu > 0)

  if (nu == 0)
    return(1)

  if (!exact) {

    y <- pi / (nu * mu) -
      pi * log(nu * mu) / (nu * mu)^2

    if (y > 1)
      stop("Asymptotic approximation invalid: consider `exact = TRUE`")

    if (nu * mu < 30)
      warning("Asymptotic approximation unreliable: consider `exact = TRUE`")

    return(y)
  }

  sigma <- mu * sqrt(2 / pi)

  f <- function(l) {
    (-expm1(-nu * l) / l) *
      exp(-l^2 / (2 * sigma^2))
  }

  psi <- integrate(
    f,
    lower = 0,
    upper = Inf,
    rel.tol = 1e-10,
    subdivisions = 1000
  )$value / sigma^2

  (pi * nu / mu - 2 * psi) / nu^2
}
