#' Calculate the probability of drawing a link of magnitude m
#'
#' @param m Integer vector of link magnitudes.
#' @param M Integer network magnitude.
#' @param exact Logical. If FALSE, use the asymptotic approximation.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

p_mag <- function(m, M, exact = TRUE) {

  # check inputs ------------------------------------------------------------

  if (any(m > M))
    stop("M must be >= m")

  if (any(m <= 0) || M <= 0)
    stop("m and M must be positive")

  if (any(m %% 1 != 0) || M %% 1 != 0)
    stop("m and M must be integers")

  if (length(M) > 1)
    stop("M must be a scalar")

  # probability calculation -------------------------------------------------

  if (exact) {
    # exact probability
    pr_m <- exp(
      -log(2 * m - 1) +
        lchoose(2 * m, m) +
        lchoose(2 * (M - m), M - m) -
        lchoose(2 * M, M)
    )
  } else {
    ## use Stirling's formula for the central binomial coefficient
    pr_m <- exp(
      -(1 / (8 * m)) +
        (1 / (192 * m^3))
    ) / ((2 * m - 1) * sqrt(pi * m))
  }

  return(pr_m)
}

#' Calculate the expected value of upstream river length
#'
#' @param lambda Numeric. Branching rate.
#' @param size Numeric. Total network length.
#' @param exact Logical. If FALSE, use the asymptotic approximation.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

u_length <- function(lambda, size, exact = TRUE) {

  if (any(c(lambda, size) <= 0))
    stop("lambda and size must be > 0")

  if (exact) {
    ## z: number of links minus 1
    ## pz: probability of b - 1 (= z) links
    m_pz <- cpois(lambda = lambda, size = size)
    z <- m_pz[ ,"z"]
    pz <- m_pz[, "pz"]

    n_z <- exp((z + 2) * log(2) - lchoose(z + 2, 0.5 * (z + 2))) - 2
    l_z <- size / (z + 1)
    u_z <- l_z * (n_z + 0.5)

    u_hat <- sum(u_z * pz)
  } else {
    ## order of approximation is ^-1/2
    u_hat <- sqrt(pi / 2) * sqrt(size / lambda) - (3 / 2) * (1 / lambda)
  }

  return(u_hat)
}

#' Calculate the expected value of network diameter
#'
#' @param lambda Numeric. Branching rate.
#' @param size Numeric. Total network length.
#' @param exact Logical. If FALSE, use the asymptotic approximation.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

diameter <- function(lambda, size, exact = TRUE) {

  if (exact) {
    m_pz <- cpois(lambda = lambda, size = size)
    z <- m_pz[ ,"z"]
    pz <- m_pz[, "pz"]

    log_v <- log(size) - log(z + 2) + (z + 2) * log(2) - lchoose(z + 2, 0.5 * (z + 2))
    d <- sum(exp(log_v) * pz)
  } else {
    d <- sqrt(pi/2) * sqrt(size / lambda)
  }

  return(d)
}

#' Calculate the expected value of pairwise network distance
#'
#' @param lambda Numeric. Branching rate.
#' @param size Numeric. Total network length.
#' @param exact Logical. If FALSE, use the asymptotic approximation.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

pdist <- function(lambda, size, exact = TRUE) {

  if (exact) {
    m_pz <- cpois(lambda = lambda, size = size)
    z <- m_pz[ ,"z"]
    pz <- m_pz[, "pz"]

    A <- log(size) + log(z + 2) - (log(z) + log(z + 1))
    B <- log(3/2) + (z + 2) * log(2) - lchoose(z + 2, 0.5 * (z + 2))
    nd <- exp(A + B) - 2 * exp(A)
  } else {
    nd <- (3/2) * sqrt(pi/2) * sqrt(size / lambda) - (2 / lambda)
  }

  return(nd)
}

#' Equilibrium occupancy for basal species
#'
#' @inheritParams u_length
#' @param h Numeric. Habitat patch density (or rate) per unit distance.
#' @param delta Numeric. Species' dispersal capability.
#' @param rsrc Numeric. Minimum establishment probability for basal species.
#' @param zeta Numeric. Effect of stream size on establishment probability for basal species.
#' @param mu Numeric. Disturbance rate.
#' @param rho Numeric. Synchrony probability of disturbance.
#' @param g Numeric. Propagule size.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

p_base <- function(lambda,
                   size,
                   h = 1,
                   delta = 1,
                   rsrc = 1,
                   zeta = 0,
                   mu = 1,
                   rho = 0.5,
                   g = 10) {

  ## check input
  l_par <- sapply(list(h, delta, rsrc, mu, rho, g),
                  function(x) length(x) > 1)

  v_par <- sapply(list(h, delta, mu, g),
                  function(x) any(x < 0))

  v_zero_one <- sapply(list(rsrc, rho),
                       function(x) any(x < 0) || any(x > 1))

  if (any(c(l_par, v_par, v_zero_one)))
    stop("invalid parameter input")

  ## define upstream river length
  u <- u_length(lambda = lambda, size = size)

  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## - s: scalar, survival probability during migration
  s <- 1 - exp(- delta * h)
  sxg <- s * g

  pgle <- ifelse(sxg < n_patch,
                 yes = sxg,
                 no = n_patch)

  ## - stream size dependency in establishment prob.
  ## - establishment probability changes with u, but truncated at zero and one
  p_r <- rsrc + zeta * u
  r0 <- max(min(c(p_r, 1)), 0)

  ## clnz: colonization rate
  clnz <- r0 * pgle

  ## extn: extinction rate
  extn <- mu * (1 + rho * u)

  ## equilibrium patch occupancy
  if (extn == 0 && clnz == 0)
    stop("colonization and extinction rates are both zero; equilibrium undefined.")

  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)

  return(p_hat)
}


#' Equilibrium occupancy for consumer species
#'
#' @inheritParams u_length
#' @inheritParams p_base
#' @param prey Numeric. Prey availability for consumer species.
#' @param max_prey Numeric. Species richness of possible prey.
#' @param mu Numeric vector of disturbance rates.
#'  \code{mu[1]} = base rate;
#'  \code{mu[2]} = prey-induced rate;
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

p_cnsm <- function(lambda,
                   size,
                   h = 1,
                   delta = 1,
                   prey,
                   max_prey,
                   mu = 1,
                   rho = 0.5,
                   g = 10) {

  ## check input
  l_par <- sapply(list(h, delta, rho, g, prey, max_prey),
                  function(x) length(x) > 1)

  v_par <- sapply(list(h, delta, g, prey, max_prey),
                  function(x) any(x < 0))

  v_zero_one <- sapply(list(rho),
                       function(x) any(x < 0) || any(x > 1))

  if (any(c(l_par, v_par, v_zero_one)))
    stop("invalid parameter input")

  ## define upstream river length
  u <- u_length(lambda = lambda, size = size)

  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## s: scalar, survival probability during migration
  s <- 1 - exp(-delta * h)
  sxg <- s * g

  pgle <- ifelse(sxg < n_patch,
                 yes = sxg,
                 no = n_patch)

  ## clnz: colonization rate
  clnz <- (prey / max_prey) * pgle

  ## extn: extinction rate
  v_mu <- to_v(mu, 2)

  extn <- v_mu[1] * (1 + rho * u) + v_mu[2] * (1 - (prey / max_prey))

  ## equilibrium patch occupancy
  if (extn == 0 && clnz == 0)
    stop("colonization and extinction rates are both zero; equilibrium undefined.")

  p_hat <- 1 - (extn / clnz)
  p_hat <- ifelse(p_hat > 0, p_hat, 0)

  return(p_hat)
}

#' Numerical solver for equilibrium occupancies
#'
#' @inheritParams u_length
#' @inheritParams p_base
#' @param w Matrix. Binary food web matrix from \code{mcbrnet::ppm()}
#' @param mu0 Numeric scalar or vector of base extinction rates.
#' @param mu_p Numeric scalar or vector of prey-induced extinction rates.
#' @param mu_c Numeric scalar or vector of consumer-induced extinction rates.
#' @param rho Numeric scalar or vector of synchrony probability of disturbance.
#' @param x0 Numeric. Initial occupancy.
#' @param n_timestep Integer. Number of time steps.
#' @param interval Numeric. Interval for numerical solver.
#' @param threshold Numeric. Threshold value for absorbing condition.
#' @param ... Additional arguments for \code{deSolve::ode()}
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

npom <- function(w,
                 size,
                 lambda,
                 h = 1,
                 delta = 1,
                 rsrc = 1,
                 zeta = 0,
                 g = 1,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 rho = 0.5,
                 x0 = 0.5,
                 n_timestep = 100,
                 interval = 0.01,
                 threshold = 1E-5,
                 ...) {

  # check input -------------------------------------------------------------

  absfwb <- abs(w)
  if (!all(absfwb == t(absfwb)))
    stop("the input w is invalid (abs(w) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input w is invalid (abs(w) must be binary)")

  l_par <- sapply(list(h, rsrc),
                  function(x) length(x) > 1)

  v_par <- sapply(list(h, delta, rsrc, g, mu0, mu_p, mu_c),
                  function(x) any(x < 0))

  v_zero_one <- sapply(list(rho, x0),
                       function(x) any(x < 0) || any(x > 1))

  if (any(c(l_par, v_par, v_zero_one)))
    stop("invalid parameter input")

  # constant setup ----------------------------------------------------------

  ## number of species
  n_species <- nrow(absfwb)

  ## prey (Mp) and consumption (Mc) matrix
  Mp <- Mc <- absfwb
  Mp[upper.tri(Mp)] <- 0
  Mc[lower.tri(Mc)] <- 0

  ## number of prey for each species
  s_prey <- rowSums(Mp)

  ## inverse of maximum prey numbers
  inv_s_prey <- ifelse(s_prey > 0,
                       yes = 1 / s_prey,
                       no = -1)

  ## number of basal species
  id_b <- which(s_prey == 0)
  n_b <- length(id_b)

  ## number of consumer species
  id_c <- which(s_prey > 0)
  n_c <- length(id_c)


  # parameter setup ---------------------------------------------------------

  ## spatial parameters
  v_rho <- to_v(rho, n = n_species)
  u <- u_length(lambda = lambda, size = size)

  ## colonization rate
  ## - propagule survival
  v_delta <- to_v(delta, n = n_species)
  v_s <- 1 - exp(-v_delta * h)

  ## - resource availability
  v_zeta <- to_v(zeta, n = n_b)
  v_rsrc <- to_v(rsrc, n = n_b)
  v_p_r <- v_rsrc + v_zeta * u
  v_r0 <- sapply(v_p_r,
                 function(y) max(min(c(y, 1)), 0)
                 )

  v_r <- c(v_r0, rep(0, n_c))

  ## - propagule
  n_patch <- h * size
  v_g <- to_v(g, n = n_species)
  v_pgle <- rep(NA, n_species)

  v_sxg <- v_s * v_g

  id_sxg <- which(v_sxg < n_patch)
  id_n_patch <- which(v_sxg >= n_patch)
  v_pgle[id_sxg] <- v_sxg[id_sxg]
  v_pgle[id_n_patch] <- n_patch

  ## extinction rate
  ## - base rate
  v_mu0 <- to_v(mu0, n = n_species)

  ## - prey availability
  v_mu_p <- to_v(mu_p, n = n_species)
  v_mu_p[id_b] <- 0

  ## - consumption (predation) influence
  if (is.matrix(mu_c)) {
    ## - if matrix
    if (any(dim(mu_c) != n_species))
      stop(paste("if mu_c is a matrix, its dimensions must be", n_species, "x", n_species))

    m_mu_c <- mu_c
  } else {
    ## - if non-matrix
    if (!(length(mu_c) == 1 || length(mu_c) == n_species))
      stop("if mu_c is non-matrix, it must be a scalar or have a length of n_species")

    m_mu_c <- to_v(mu_c, n = n_species)
  }

  # run ode -----------------------------------------------------------------

  ## derivative
  derivr <- function(t, x, parms) {
    with(parms, {

      ## - colonization
      clnz <- pgle * ((Mp %*% x) * inv_s_prey + r)

      ## - extinction
      extn <- mu0 * (1 + rho * u) +
        mu_p * (1 - (Mp %*% x) * inv_s_prey) +
        mu_c * Mc %*% x

      dx <- clnz * x * (1 - x) - extn * x

      list(dx)
    })
  }

  ## set parameters for ode()
  parms <- list(pgle = v_pgle,
                g = v_g,
                Mp = Mp,
                r = v_r,
                mu0 = v_mu0,
                mu_p = v_mu_p,
                inv_s_prey = inv_s_prey,
                mu_c = m_mu_c,
                Mc = Mc,
                rho = v_rho,
                u = u)

  x_init <- to_v(x0, n_species)
  times <- seq(0, n_timestep, by = interval)

  ## define absorbing condition
  ## - root function
  rootfun <- function(t, x, parms) {
    return(x - threshold)
  }

  ## - extinction: triggered when "x - threshold = 0"
  eventfun <- function(t, x, parms) {
    x <- ifelse(x <= threshold, 0, x)
    return(x)
  }

  # run ode solver
  cout <- deSolve::ode(y = x_init,
                       times = times,
                       func = derivr,
                       parms = parms,
                       events = list(func = eventfun,
                                     root = TRUE),
                       rootfun = rootfun,
                       ...)

  return(cout)
}

#' Equilibrium food chain length (analytical)
#'
#' @inheritParams u_length
#' @inheritParams npom
#' @param weight Logical.
#'  If \code{TRUE}, maximum trophic position is weighted by relative occupancies.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

fcl <- function(w,
                lambda,
                size,
                h = 1,
                delta = 1,
                rsrc = 1,
                zeta = 0,
                g = 10,
                mu0 = 0.1,
                mu_p = 0.1,
                rho = 0.5,
                weight = TRUE) {

  # check input -------------------------------------------------------------

  absfwb <- abs(w)
  if (!all(absfwb == t(absfwb)))
    stop("the input w is invalid (abs(w) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input w is invalid (abs(w) must be binary)")


  # transform input ---------------------------------------------------------

  ## w: matrix, consumer-resource matrix. produce with ppm()
  fwb <- absfwb
  fwb[lower.tri(fwb)] <- 0
  max_prey <- colSums(fwb)

  ## constant terms, rsrc, zeta for basal species
  ## - create vectors with n-basal elements
  n_b <- sum(max_prey == 0)
  list_b <- lapply(list(rsrc, zeta),
                   FUN = function(x) to_v(x, n = n_b))

  names(list_b) <- c("rsrc",
                     "zeta")


  ## constant terms, delta, g, mu0, mu_p, rho
  ## - create vectors with n-species elements
  n_sp <- unique(dim(w))
  list_all <- lapply(list(delta, g, mu0, mu_p, rho),
                     FUN = function(x) to_v(x, n = n_sp))

  names(list_all) <- c("delta",
                       "g",
                       "mu0",
                       "mu_p",
                       "rho")

  list_parms <- c(list_b, list_all)

  ## p_hat: vector initialized with -1, equilibrium occupancy
  ## max_prey: vector, maximum number of prey items for consumer j
  p_hat <- rep(-1, n_sp)

  # occupancies -------------------------------------------------------------

  ## sequential determination of equilibrium occupancies
  for (j in seq_len(n_sp)) {

    if (max_prey[j] == 0) {
      ## basal species
      p_hat[j] <- with(list_parms,
                       p_base(lambda = lambda,
                              size = size,
                              h = h,
                              delta = delta[j],
                              rsrc = rsrc[j],
                              zeta = zeta[j],
                              mu = mu0[j],
                              rho = rho[j],
                              g = g[j])
                       )

    } else {
      ## consumers

      ## index of prey species for consumer j
      index_prey <- which(fwb[, j] == 1)

      ## mean-field prey richness
      prey <- sum(p_hat[index_prey])

      ## possible maximum of prey richness
      n_prey <- max_prey[j]

      p_hat[j] <- with(list_parms,
                       p_cnsm(lambda = lambda,
                              size = size,
                              h = h,
                              delta = delta[j],
                              prey = prey,
                              max_prey = n_prey,
                              mu = c(mu0[j], mu_p[j]),
                              rho = rho[j],
                              g = g[j])
                       )

    } # ifelse
  } # for j


  # food chain length -------------------------------------------------------

  fcl <- maxtp(w = w,
               occupancy = p_hat,
               weight = weight)

  attr(fcl, "p_hat") <- p_hat

  return(fcl)
}


#' Equilibrium food chain length (numerical)
#'
#' @inheritParams u_length
#' @inheritParams npom
#' @param n_plus Number of additional runs to check convergence to equilibrium.
#' @param weight Logical.
#'  If \code{TRUE}, maximum trophic position is weighted by relative occupancies.
#' @param tol Numeric.
#'  Tolerance value for convergence.
#'  If the difference in the final values of
#'  the main and additional runs is less than the tolerance value,
#'  the function returns successful convergence as \code{0}
#'  in attribute "convergence," otherwise \code{1}.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

nfcl <- function(w,
                 lambda,
                 size,
                 h = 1,
                 delta = 1,
                 rsrc = 1,
                 zeta = 0,
                 g = 10,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 rho = 0.5,
                 x0 = 0.5,
                 n_timestep = 100,
                 interval = 0.01,
                 threshold = 1e-05,
                 n_plus = 10,
                 weight = TRUE,
                 tol = 1e-06,
                 ...) {

  # numerical solution ------------------------------------------------------

  ## main run
  cout <- npom(w = w,
               lambda = lambda,
               size = size,
               h = h,
               delta = delta,
               rsrc = rsrc,
               zeta = zeta,
               g = g,
               mu0 = mu0,
               mu_p = mu_p,
               mu_c = mu_c,
               rho = rho,
               x0 = x0,
               n_timestep = n_timestep,
               interval = interval,
               threshold = threshold,
               ...)

  p_hat <- cout[nrow(cout), -1]
  p_hat[p_hat < 0] <- 0
  p_hat[p_hat > 1] <- 1

  ## additional run to check equilibrium
  cout_plus <- npom(w = w,
                    lambda = lambda,
                    size = size,
                    h = h,
                    delta = delta,
                    rsrc = rsrc,
                    zeta = zeta,
                    g = g,
                    mu0 = mu0,
                    mu_p = mu_p,
                    mu_c = mu_c,
                    rho = rho,
                    x0 = c(p_hat),
                    n_timestep = n_plus,
                    interval = interval,
                    threshold = threshold)

  p_hat_plus <- cout_plus[nrow(cout_plus), -1]

  ## check difference
  z <- abs(p_hat - p_hat_plus)
  conv <- ifelse(any(z > tol), 1, 0)

  # food chain length -------------------------------------------------------

  ## in case npom() returns zeros as non-zero values (e.g, 1e-30)
  if (threshold > 0) p_hat <- floor(p_hat / threshold) * threshold

  fcl <- maxtp(w = w,
               occupancy = p_hat,
               weight = weight)

  attr(fcl, "p_hat") <- p_hat
  attr(fcl, "convergence") <- conv

  return(fcl)
}
