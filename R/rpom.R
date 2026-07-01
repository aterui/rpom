#' Probability of drawing a link with magnitude m
#'
#' Calculates the probability that a randomly selected link has magnitude
#' `m` in a branching river network of magnitude `M`.
#'
#' @param m Integer vector of link magnitudes.
#' @param M Integer. River network magnitude.
#' @param exact Logical. If `TRUE`, compute the exact probability. If `FALSE`,
#'   use the asymptotic approximation.
#'
#' @return A numeric vector of probabilities corresponding to `m`.
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

#' Expected upstream river length
#'
#' Calculates the expected total length of the upstream river network from a
#' randomly selected location in a branching river network.
#'
#' @param lambda Numeric. Branching rate of the river network.
#' @param size Numeric. Total river network length.
#' @param exact Logical. If `TRUE`, compute the exact expectation. If `FALSE`,
#'   use the asymptotic approximation.
#'
#' @return A numeric scalar giving the expected upstream river length.
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
    z <- m_pz[, "z"]
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

#' Expected network diameter
#'
#' Calculates the expected diameter of a branching river network.
#'
#' @param lambda Numeric. Branching rate of the river network.
#' @param size Numeric. Total river network length.
#' @param exact Logical. If `TRUE`, compute the exact expectation. If `FALSE`,
#'   use the asymptotic approximation.
#'
#' @return A numeric scalar giving the expected network diameter.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

diameter <- function(lambda, size, exact = TRUE) {

  if (any(c(lambda, size) <= 0))
    stop("lambda and size must be > 0")

  if (exact) {
    m_pz <- cpois(lambda = lambda, size = size)
    z <- m_pz[, "z"]
    pz <- m_pz[, "pz"]

    log_v <- log(size) - log(z + 2) + (z + 2) * log(2) - lchoose(z + 2, 0.5 * (z + 2))
    d <- sum(exp(log_v) * pz)
  } else {
    d <- sqrt(pi / 2) * sqrt(size / lambda)
  }

  return(d)
}

#' Expected pairwise network distance
#'
#' Calculates the expected distance between two randomly selected locations in
#' a branching river network.
#'
#' @param lambda Numeric. Branching rate of the river network.
#' @param size Numeric. Total river network length.
#' @param exact Logical. If `TRUE`, compute the exact expectation. If `FALSE`,
#'   use the asymptotic approximation.
#'
#' @return A numeric scalar giving the expected pairwise network distance.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

pdist <- function(lambda, size, exact = TRUE) {

  if (any(c(lambda, size) <= 0))
    stop("lambda and size must be > 0")

  if (exact) {
    m_pz <- cpois(lambda = lambda, size = size)
    z <- m_pz[, "z"]
    pz <- m_pz[, "pz"]

    ## "between" links cases for z >= 2 (>= 3 links present)
    v <- numeric(length(z))
    idx <- z > 0

    A <- log(size) + log(z[idx] + 2) - log(z[idx]) - log(z[idx] + 1)

    B <- (z[idx] + 2) * log(2) - lchoose(z[idx] + 2, 0.5 * (z[idx] + 2))

    v[idx] <- exp(A + B) - 2 * exp(A)

    pw <- z / (z + 2)

    nd_z <- pw * v + (1 - pw) * (2 * size / ((z + 2) * (z + 3)))
    nd <- sum(nd_z * pz)
  } else {
    nd <- sqrt(pi / 2) * sqrt(size / lambda) - (2 / lambda)
  }

  return(nd)
}

#' Equilibrium occupancy of a basal species
#'
#' Computes the equilibrium patch occupancy of a basal species in a branching
#' river network under colonization–extinction dynamics influenced by network
#' structure.
#'
#' @inheritParams u_length
#' @param h Numeric. Habitat patch density (patches per unit stream length).
#' @param delta Numeric. Distance-decay rate of propagule dispersal.
#' @param r0 Numeric. Baseline establishment probability (0--1).
#' @param b Numeric. Effect of upstream river length on establishment
#'   probability.
#' @param mu Numeric. Baseline extinction rate.
#' @param rho0 Numeric. Baseline synchrony probability.
#' @param nu Numeric. Distance-decay rate of disturbance synchrony.
#' @param kernel Character. Functional form of distance decay in disturbance
#'   synchrony. Either `"exp"` for exponential decay or `"linear"` for
#'   linear decay with distance.
#' @param g Numeric. Propagule production rate.
#' @param exact Logical. If `TRUE`, use exact network calculations. If `FALSE`,
#'   use asymptotic approximations.
#'
#' @return A numeric scalar giving the equilibrium occupancy.
#'
#' @author Akira Terui
#' @export

p_base <- function(lambda,
                   size,
                   h = 1,
                   delta = 0.1,
                   r0 = 1,
                   b = 0,
                   mu = 1,
                   rho0 = 1,
                   nu = 0,
                   kernel = c("exp", "linear"),
                   g = 1,
                   exact = FALSE) {

  ## check input (basic scalar/validity checks)
  l_par <- sapply(list(h, delta, r0, b, mu, nu, g),
                  function(x) length(x) > 1)

  if (any(l_par))
    stop("All parameters must be scalar.")

  if (any(c(h, delta, r0, b, mu, nu, g) < 0))
    stop("All parameters must be non-negative.")

  if (r0 < 0 || r0 > 1)
    stop("r0 must be between 0 and 1.")

  if (rho0 < 0 || rho0 > 1)
    stop("rho0 must be between 0 and 1.")

  kernel <- match.arg(kernel)

  ## define upstream river length
  u <- u_length(lambda = lambda,
                size = size,
                exact = exact)

  ## pairwise distance
  d <- pdist(lambda = lambda,
             size = size,
             exact = exact)

  ## network diameter
  diam <- diameter(lambda = lambda,
                   size = size,
                   exact = exact)

  if (any(c(u, d, diam) < 0))
    stop("Invalid approximation; increase lambda * size")

  ## number of habitat patches
  n_patch <- h * size

  ## propagule pressure (ensure non-negative)
  pgle <- (g * n_patch) * laplace_rayleigh(delta = delta, mu = d)

  if (pgle < 0 || is.nan(pgle))
    stop("pgle = ", pgle, "; invalid parameter values")

  ## establishment probability (bounded 0–1)
  r <- r0 + b * u
  if (r < 0 || r > 1)
    stop("r = ", r, "; invalid parameter values")

  ## disturbance synchrony (bounded 0–1)
  rho <- switch(
    kernel,
    linear = rho0 * (1 - nu * (diam / 3)),
    exp = rho0 * laplace_rt(nu = nu, mu = diam, exact = exact),
    stop("Unknown kernel type: ", kernel)
  )

  if (rho < 0 || rho > 1)
    stop("rho = ", rho, "; invalid parameter values")

  ## colonization and extinction
  clnz <- r * pgle
  extn <- mu * (1 + rho * u)

  if (clnz == 0 && extn == 0)
    stop("Colonization and extinction are both zero; equilibrium undefined.")

  ## equilibrium occupancy
  p_hat <- 1 - (extn / clnz)
  p_hat <- max(p_hat, 0)

  return(p_hat)
}

#' Equilibrium occupancy of a consumer species
#'
#' Computes the equilibrium patch occupancy of a consumer species in a
#' branching river network given the equilibrium occupancy of its prey.
#'
#' @inheritParams u_length
#' @inheritParams p_base
#' @param prey Numeric vector of equilibrium occupancy probabilities for prey
#'   species.
#' @param max_prey Numeric. Total number of potential prey species.
#' @param mu Numeric. Extinction-rate parameter(s). If a scalar, the same value
#'   is used for both extinction components. If a length-2 vector,
#'   `mu[1]` is the baseline extinction rate and `mu[2]` is the
#'   prey-dependent extinction rate.
#'
#' @return A numeric scalar giving the equilibrium occupancy.
#'
#' @author Akira Terui
#'
#' @export

p_cnsm <- function(lambda,
                   size,
                   h = 1,
                   delta = 0.1,
                   prey,
                   max_prey,
                   mu = 1,
                   rho0 = 1,
                   nu = 0,
                   kernel = c("exp", "linear"),
                   g = 1,
                   exact = FALSE) {

  ## check input (basic scalar/validity checks)
  l_par <- sapply(list(h, delta, max_prey, nu, g),
                  function(x) length(x) > 1)

  if (any(l_par))
    stop("All parameters but 'mu' must be scalar.")

  if (any(c(h, delta, prey, max_prey, mu, nu, g) < 0))
    stop("All parameters must be non-negative.")

  if (!(length(mu) %in% c(1, 2)))
    stop("'mu' must have length 1 or 2.")

  if (rho0 < 0 || rho0 > 1)
    stop("rho0 must be between 0 and 1.")

  kernel <- match.arg(kernel)

  ## define upstream river length
  u <- u_length(lambda = lambda,
                size = size,
                exact = exact)

  ## pairwise distance
  d <- pdist(lambda = lambda,
             size = size,
             exact = exact)

  ## network diameter
  diam <- diameter(lambda = lambda,
                   size = size,
                   exact = exact)

  if (any(c(u, d, diam) < 0))
    stop("Invalid approximation; increase lambda * size")

  ## number of habitat patches
  n_patch <- h * size

  ## propagule pressure (ensure non-negative)
  pgle <- (g * n_patch) * laplace_rayleigh(delta = delta, mu = d)

  if (pgle < 0 || is.nan(pgle))
    stop("pgle = ", pgle, "; invalid parameter values")

  ## summed prey occupancy
  s <- sum(prey)
  if (s > max_prey)
    stop("Summed prey occupancy = ", sum(prey), "; must be smaller than 'max_prey'")

  ## fraction of colonizable habitat with at least one prey
  ## log1p(-prey) = log(1 - prey)
  log_eta <- sum(log1p(-prey))
  eta <- 1 - exp(log_eta)

  ## disturbance synchrony (bounded 0–1)
  rho <- switch(
    kernel,
    linear = rho0 * (1 - nu * (diam / 3)),
    exp = rho0 * laplace_rt(nu = nu, mu = diam, exact = exact),
    stop("Unknown kernel type: ", kernel)
  )

  if (rho < 0 || rho > 1)
    stop("rho = ", rho, "; invalid parameter values")

  ## clnz: colonization rate
  clnz <- (s / max_prey) * pgle

  ## extn: extinction rate
  v_mu <- to_v(mu, 2)
  extn <- v_mu[1] * (1 + rho * u) + v_mu[2] * (1 - (s / max_prey))

  ## equilibrium patch occupancy
  if (extn == 0 && clnz == 0)
    stop("colonization and extinction rates are both zero; equilibrium undefined.")

  p_hat <- eta - (extn / clnz)
  p_hat <- max(p_hat, 0)

  return(p_hat)
}

#' Equilibrium food chain length
#'
#' Computes the equilibrium food chain length of a branching river food web
#' using the analytical colonization--extinction model.
#'
#' @inheritParams u_length
#' @inheritParams npom
#' @param weight Logical. If `TRUE`, the maximum trophic position is weighted
#'   by equilibrium species occupancies. If `FALSE`, food chain length is the
#'   maximum trophic position.
#'
#' @return A numeric scalar giving the equilibrium food chain length. The
#'   equilibrium occupancies are stored as the `"p_hat"` attribute of the
#'   returned value.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

fcl <- function(w,
                lambda,
                size,
                h = 1,
                delta = 0.1,
                r0 = 1,
                b = 0,
                rho0 = 1,
                nu = 0,
                kernel = c("exp", "linear"),
                g = 1,
                mu0 = 0.1,
                mu_p = 0.1,
                weight = TRUE,
                exact = FALSE) {

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

  ## constant terms, r0, b for basal species
  ## - create vectors with n-basal elements
  n_b <- sum(max_prey == 0)
  list_b <- lapply(list(r0 = r0,
                        b = b),
                   FUN = function(x) to_v(x, n = n_b))

  ## constant terms, delta, g, mu0, mu_p, rho
  ## - create vectors with n-species elements
  n_sp <- unique(dim(w))
  list_all <- lapply(list(delta = delta,
                          g = g,
                          mu0 = mu0,
                          mu_p = mu_p,
                          rho0 = rho0,
                          nu = nu,
                          kernel = kernel),
                     FUN = function(x) to_v(x, n = n_sp))

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
                              r0 = r0[j],
                              b = b[j],
                              mu = mu0[j],
                              rho0 = rho0[j],
                              nu = nu[j],
                              kernel = kernel[j],
                              g = g[j],
                              exact = exact)
      )

    } else {
      ## consumers

      ## index of prey species for consumer j
      index_prey <- which(fwb[, j] == 1)

      ## mean-field prey richness
      prey <- p_hat[index_prey]
      if (any(prey < 0) || any(prey > 1))
        stop("Invalid prey occupancy.")

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
                              rho0 = rho0[j],
                              nu = nu[j],
                              kernel = kernel[j],
                              g = g[j],
                              exact = exact)
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

#' Numerical solver for equilibrium occupancies
#'
#' Solves the colonization--extinction model numerically to obtain equilibrium
#' occupancies of all species in a branching river food web.
#'
#' @inheritParams u_length
#' @inheritParams p_base
#' @param w Matrix. Binary food web matrix produced by `ecotools::ppm()`.
#' @param mu0 Numeric scalar or vector. Baseline extinction rate.
#' @param mu_p Numeric scalar or vector. Prey-dependent extinction rate.
#' @param mu_c Numeric scalar, vector, or matrix. Consumer-dependent extinction
#'   rate. If a matrix, each element specifies the per-capita effect of a
#'   consumer on its prey.
#' @param x0 Numeric scalar or vector. Initial occupancy.
#' @param n_timestep Integer. Number of time steps for numerical integration.
#' @param intv Numeric. Time interval between integration steps.
#' @param threshold Numeric. Occupancy threshold below which a species is
#'   considered extinct.
#' @param ... Additional arguments passed to `deSolve::ode()`.
#'
#' @return A matrix returned by `deSolve::ode()`, with one row per time step and
#'   one column per species (plus the time column).
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

npom <- function(w,
                 size,
                 lambda,
                 h = 1,
                 delta = 1,
                 r0 = 1,
                 b = 0,
                 g = 1,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 rho0 = 1,
                 nu = 0,
                 kernel = c("exp", "linear"),
                 x0 = 0.5,
                 n_timestep = 100,
                 intv = 0.01,
                 threshold = 1E-5,
                 exact = TRUE,
                 ...) {

  # check input -------------------------------------------------------------

  absfwb <- abs(w)

  l_par <- sapply(list(h, r0, b), function(x) length(x) > 1)

  zo <- any(x0 < 0) || any(x0 > 1)

  if (!all(absfwb == t(absfwb)))
    stop("Input w is invalid (abs(w) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("Input w is invalid (abs(w) must be binary)")

  if (any(c(lambda, size, h, delta, r0, g, mu0, mu_p, mu_c, nu) < 0))
    stop("All parameters must be positive")

  if (any(l_par))
    stop("Parameters h, r0, b must be scalar input")

  if (zo)
    stop("x0 must be a fraction, i.e., x0 in [0, 1]")

  if (rho0 < 0 || rho0 > 1)
    stop("rho0 must be between 0 and 1.")

  kernel <- match.arg(kernel)

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

  ## geometric parameters
  u <- u_length(lambda = lambda, size = size, exact = exact)
  diam <- diameter(lambda = lambda, size = size, exact = exact)
  d <- pdist(lambda = lambda, size = size, exact = exact)

  ## spatial parameters
  v_nu <- to_v(nu, n = n_species)
  v_laplace_rt <- Vectorize(laplace_rt)

  v_rho <-
    switch(
      kernel,
      linear = rho0 * (1 - v_nu * (diam / 3)),
      exp = rho0 * v_laplace_rt(nu = v_nu, mu = diam, exact = exact),
      stop("Unknown kernel type: ", kernel)
    )

  if (any(c(v_rho < 0, v_rho > 1)))
    stop("rho must be a probability i.e., rho in [0, 1].")

  ## colonization rate
  ## - propagule survival
  v_delta <- to_v(delta, n = n_species)
  v_s <- laplace_rayleigh(delta = v_delta, mu = d)

  ## - resource availability
  v_b <- to_v(b, n = n_b)
  v_r0b <- to_v(r0, n = n_b)
  v_rb <- v_r0b + v_b * u

  if (any(c(v_rb < 0, v_rb > 1)))
    stop("r (= r0 + b * u) must be a probability i.e., r in [0, 1].")

  v_r <- c(v_rb, rep(0, n_c))

  ## - propagule
  n_patch <- h * size
  v_g <- to_v(g, n = n_species)
  v_pgle <- v_g * n_patch * v_s

  ## extinction rate
  ## - base rate
  v_mu0 <- to_v(mu0, n = n_species)

  ## - prey availability effect (set zero for basal species)
  v_mu_p <- to_v(mu_p, n = n_species)
  v_mu_p[id_b] <- 0

  ## - consumption (predation) influence
  if (is.matrix(mu_c)) {
    ## - if matrix
    if (any(dim(mu_c) != n_species))
      stop("if mu_c is a matrix, its dimensions must be ", n_species, " by ", n_species)

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
      extn <-
        mu0 * (1 + rho * u) +
        mu_p * (1 - (Mp %*% x) * inv_s_prey) +
        mu_c * Mc %*% x

      ## eta
      ## - subtract small value `xi` to avoid exact 1
      ## - basal species's eta = 1
      xi <- 1e-12
      k <- pmin(x, 1 - xi)

      eta <- 1 - exp(Mp %*% log1p(-k))
      eta[id_b] <- 1

      ## dx/dt
      dx <- clnz * x * (eta - x) - extn * x

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
                u = u,
                id_b = id_b)

  x_init <- to_v(x0, n_species)
  times <- seq(0, n_timestep, by = intv)

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

#' Equilibrium food chain length
#'
#' Computes the equilibrium food chain length of a branching river food web
#' using the numerical colonization--extinction model.
#'
#' @inheritParams u_length
#' @inheritParams npom
#' @param n_plus Integer. Number of additional time steps used to assess
#'   convergence to equilibrium.
#' @param weight Logical. If `TRUE`, the maximum trophic position is weighted
#'   by equilibrium species occupancies. If `FALSE`, food chain length is the
#'   maximum trophic position.
#' @param tol Numeric. Convergence tolerance. Convergence is considered
#'   successful if the absolute difference between the final occupancies from
#'   the main and additional runs is less than `tol` for all species.
#'
#' @return A numeric scalar giving the equilibrium food chain length. The
#'   equilibrium occupancies are stored as the `"p_hat"` attribute. The
#'   `"convergence"` attribute equals `0` if convergence was achieved and `1`
#'   otherwise.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

nfcl <- function(w,
                 lambda,
                 size,
                 h = 1,
                 delta = 1,
                 r0 = 1,
                 b = 0,
                 g = 10,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 rho0 = 1,
                 nu = 0,
                 kernel = c("exp", "linear"),
                 x0 = 0.5,
                 n_timestep = 100,
                 intv = 0.01,
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
               r0 = r0,
               b = b,
               g = g,
               mu0 = mu0,
               mu_p = mu_p,
               mu_c = mu_c,
               rho0 = rho0,
               nu = nu,
               kernel = kernel,
               x0 = x0,
               n_timestep = n_timestep,
               intv = intv,
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
                    r0 = r0,
                    b = b,
                    g = g,
                    mu0 = mu0,
                    mu_p = mu_p,
                    mu_c = mu_c,
                    rho0 = rho0,
                    nu = nu,
                    kernel = kernel,
                    x0 = c(p_hat),
                    n_timestep = n_plus,
                    intv = intv,
                    threshold = threshold,
                    ...)

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
