#' Calculate a probability of drawing a link with m magnitude
#'
#' @param m Integer. Magnitude of a given link.
#' @param M Integer. Magnitude of a network.
#'
#' @importFrom stats dpois qpois
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

p_mag <- function(m, M) {

  # check inputs ------------------------------------------------------------

  if (any(m > M))
    stop("M must be >= m")

  if (any(m <= 0) || M <= 0)
    stop("m and M must be positive integer")

  if (length(M) > 1)
    stop("M must be a scalar")

  # probability calculation -------------------------------------------------

  if (M <= 500) {
    ## if M <= 500, exact calculation
    x <- (2 * m - 1)^(-1) * choose(2 * m, m) * choose(2 * (M - m), M - m)
    pr_m <- x / choose(2 * M, M)
  } else {
    ## if M > 500, approximation as if M -> infinity
    pr_m <- 2 ^ (-(2 * m - 1)) * (2 * m - 1) ^ (-1) * choose(2 * m - 1, m)
  }

  return(pr_m)
}

#' Calculate an expected value of upstream distance
#'
#' @param lambda Numeric. Branching rate of a network.
#' @param size Numeric. Total length of a network.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

u_length <- function(lambda, size) {
  ## check input
  if (lambda < 0)
    stop("invalid input: lambda must be >= 0")

  if (size <= 0)
    stop("invalid input: size must be > 0")

  ## z: number of links minus 1
  ## pr_z: probability of b - 1 (= z) links
  pois_max <- stats::qpois(1 - 1e-10, lambda = lambda * size)
  v_z <- 0:pois_max
  pr_z <- stats::dpois(v_z, lambda = lambda * size)

  ## pr_z_tr: truncate probabilities for z taking odd numbers
  ## - note, when z is an odd number, n links is an even number
  even_id <- which(v_z %% 2 == 0)
  odd_id <- which(v_z %% 2 == 1)
  pr_even <- sum(pr_z[even_id])
  pr_z_tr <- pr_z / pr_even
  pr_z_tr[odd_id] <- 0

  ## expected total length of upstream links, conditional on z
  u_z <- sapply(v_z, function(z) {
    if (z %% 2 == 0) {
      ## when z is even = n links is odd

      ## b: number of links
      ## l_hat: expected length of a link/branch, conditional on b
      ## - l_hat derived from a Beta distribution Beta(1, z)
      b <- (z + 1)
      l_hat <- size / b

      ## maximum magnitude in a network
      M <- 0.5 * (b + 1)

      if (M > 514)
        stop(paste0("Stream magnitude M exceeds 514,
                    which will return Inf in choose(2 * M, M);
                    consider smaller values of lambda and/or size"))

      ## weighted values for the number of upstream links
      ## - `2m - 2` is the number of upstream links, ub
      ## - weighted by probability of drawing a link with m magnitude, w_ub
      m <- 1:M
      w_ub <- (2 * m - 2) * p_mag(m, M)

      ## expected value for the number of upstream links ub
      ub_hat <- sum(w_ub)

      ## expected value for upstream stream length
      ## - if poisson distributed,
      ## - the arrival time (distance to a given patch) will be a uniform dist
      ## - thus, expectation is (max - min) / 2 = l / 2
      u <- ub_hat * l_hat + 0.5 * l_hat
    } else {
      ## when z is odd = n links is even
      u <- -1
    }

    return(u)
  })

  ## sum over z to get an expected value of upstream link length
  ## u = -1 will be cancelled by multiplying pr_z_tr = 0
  u_hat <- sum(pr_z_tr * u_z)

  return(u_hat)
}


#' Equilibrium occupancy for basal species
#'
#' @inheritParams u_length
#' @param h Numeric. Habitat patch density (or rate) per unit distance.
#' @param delta Numeric. Species' dispersal capability.
#' @param rsrc Numeric. Establishment probability for basal species.
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

  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## - s: scalar, survival probability during migration
  s <- 1 - exp(- delta * h)
  pgle <- ifelse(s * g < n_patch,
                 yes = s * g,
                 no = n_patch)

  ## clnz: colonization rate
  clnz <- rsrc * pgle

  ## extn: extinction rate
  extn <- mu * (1 + rho * u_length(lambda = lambda, size = size))

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


  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## s: scalar, survival probability during migration
  s <- 1 - exp(-delta * h)

  pgle <- ifelse(s * g < n_patch,
                 yes = s * g,
                 no = n_patch)

  ## clnz: colonization rate
  clnz <- (prey / max_prey) * pgle

  ## extn: extinction rate
  v_mu <- to_v(mu, 2)

  extn <-
    v_mu[1] * (1 + rho * u_length(lambda = lambda, size = size)) +
    v_mu[2] * (1 - (prey / max_prey))

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
#' @param foodweb Matrix. Binary food web matrix from \code{mcbrnet::ppm()}
#' @param mu0 Numeric scalar or vector of base extinction rates.
#' @param mu_p Numeric scalar or vector of prey-induced extinction rates.
#' @param mu_c Numeric scalar or vector of consumer-induced extinction rates.
#' @param rho Numeric scalar or vector of synchrony probability of disturbance.
#' @param x0 Numeric. Initial occupancy.
#' @param n_timestep Integer. Number of time steps.
#' @param interval Numeric. Interval for numerical solver.
#' @param threshold Numeric. Threshold value for absorbing condition.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

npom <- function(foodweb,
                 size,
                 lambda,
                 h = 1,
                 delta = 1,
                 rsrc = 1,
                 g = 1,
                 mu0 = 1,
                 mu_p = 1,
                 mu_c = 1,
                 rho = 0.5,
                 x0 = 0.5,
                 n_timestep = 100,
                 interval = 0.01,
                 threshold = 1E-5) {

  # check input -------------------------------------------------------------

  absfwb <- abs(foodweb)
  if (!all(absfwb == t(absfwb)))
    stop("the input foodweb is invalid (abs(foodweb) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input foodweb is invalid (abs(foodweb) must be binary)")

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
  n_species <- nrow(foodweb)

  ## prey (Mp) and consumption (Mc) matrix
  Mp <- Mc <- abs(foodweb)
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

  ## colonization rate
  ## - propagule survival
  v_delta <- to_v(delta, n = n_species)
  v_s <- 1 - exp(-v_delta * h)

  ## - resource availability
  r0 <- to_v(rsrc, n = n_b)
  v_r <- c(r0, rep(0, n_c))

  ## - propagule
  n_patch <- h * size
  v_g <- to_v(g, n = n_species)
  v_phi <- rep(NA, n_species)

  v_sxg <- v_s * v_g

  id_sxg <- which(v_sxg < n_patch)
  id_n_patch <- which(v_sxg >= n_patch)
  v_phi[id_sxg] <- v_sxg[id_sxg]
  v_phi[id_n_patch] <- n_patch

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
      stop(paste("mu_c must be", n_species, "x", n_species, "matrix"))

    m_mu_c <- mu_c
  } else {
    ## - if non-matrix
    if (length(mu_c) != 1)
      stop("mu_c must be a scalar or matrix")

    m_mu_c <- mu_c
  }

  ## - spatial
  v_rho <- to_v(rho, n = n_species)
  u <- u_length(lambda = lambda, size = size)

  # run ode -----------------------------------------------------------------

  ## derivative
  derivr <- function(t, x, parms) {
    with(parms, {

      ## - colonization
      clnz <- phi * ((Mp %*% x) * inv_s_prey + r)

      ## - extinction
      extn <-
        mu0 * (1 + rho * u) +
        mu_p * (1 - (Mp %*% x) * inv_s_prey) +
        mu_c * Mc %*% x

      dx <- clnz * x * (1 - x) - extn * x

      list(dx)
    })
  }

  ## set parameters for ode()
  parms <- list(phi = v_phi,
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
                       rootfun = rootfun)

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

fcl <- function(foodweb,
                lambda,
                size,
                h = 1,
                delta = 1,
                rsrc = 1,
                g = 10,
                mu0 = 0.1,
                mu_p = 0.1,
                rho = 0.5,
                weight = TRUE) {

  # check input -------------------------------------------------------------

  absfwb <- abs(foodweb)
  if (!all(absfwb == t(absfwb)))
    stop("the input foodweb is invalid (abs(foodweb) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input foodweb is invalid (abs(foodweb) must be binary)")


  # transform input ---------------------------------------------------------

  ## foodweb: matrix, consumer-resource matrix. produce with ppm()
  fwb <- absfwb
  fwb[lower.tri(fwb)] <- 0

  ## constant terms, delta, rsrc, g, mu0, mu_p, rho
  ## - create vectors with n-species elements
  n_sp <- unique(dim(foodweb))
  list_parms <- lapply(list(delta, g, mu0, mu_p, rho),
                       FUN = to_v, n_sp)

  names(list_parms) <- c("delta",
                         "g",
                         "mu0",
                         "mu_p",
                         "rho")

  ## p_hat: vector initialized with -1, equilibrium occupancy
  ## max_prey: vector, maximum number of prey items for consumer j
  p_hat <- rep(-1, n_sp)
  max_prey <- colSums(fwb)

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
                              rsrc = rsrc,
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

  fcl <- maxtp(foodweb = foodweb,
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

nfcl <- function(foodweb,
                 lambda,
                 size,
                 h = 1,
                 delta = 1,
                 rsrc = 1,
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
                 tol = 1e-06) {

  # numerical solution ------------------------------------------------------

  ## main run
  cout <- npom(foodweb = foodweb,
               lambda = lambda,
               size = size,
               h = h,
               delta = delta,
               rsrc = rsrc,
               g = g,
               mu0 = mu0,
               mu_p = mu_p,
               mu_c = mu_c,
               rho = rho,
               x0 = x0,
               n_timestep = n_timestep,
               interval = interval,
               threshold = threshold)

  p_hat <- cout[nrow(cout), -1]
  p_hat[p_hat < 0] <- 0
  p_hat[p_hat > 1] <- 1

  ## additional run to check equilibrium
  cout_plus <- npom(foodweb = foodweb,
                    lambda = lambda,
                    size = size,
                    h = h,
                    delta = delta,
                    rsrc = rsrc,
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

  fcl <- maxtp(foodweb = foodweb,
               occupancy = p_hat,
               weight = weight)

  attr(fcl, "p_hat") <- p_hat
  attr(fcl, "convergence") <- conv

  return(fcl)
}
