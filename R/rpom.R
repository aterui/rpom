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
#' @param rsrc Resource availability for basal species.
#' @param mu Numeric vector of disturbance rates:
#'  \code{mu[1]} = base rate;
#'  \code{mu[2]} = spatial rate
#' @param rho Synchrony probability of disturbance.
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
                   mu = c(1, 1),
                   rho = 0.5,
                   g = 100) {

  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  s <- 1 - exp(- delta * h)
  pgle <- ifelse(s * g < n_patch,
                 yes = s * g,
                 no = n_patch)
  clnz <- rsrc * pgle

  ## extn: extinction rate
  if (length(mu) == 1) {
    mu <- rep(mu, 2)
  } else {
    if (length(mu) != 2) stop("error in mu")
  }

  extn <- mu[1] + mu[2] * rho * u_length(lambda = lambda, size = size)

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
#' @param prey Resource availability for basal species.
#' @param max_prey Resource availability for basal species.
#' @param mu Numeric vector of disturbance rates.
#'  \code{mu[1]} = base rate;
#'  \code{mu[2]} = prey rate;
#'  \code{mu[3]} = spatial rate.
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
                   mu = c(1, 1, 1),
                   rho = 0.5,
                   g = 100) {

  ## n_patch: scalar, # habitat patches
  n_patch <- h * size

  ## clnz: colonization rate
  ## - s: scalar, survival probability during migration
  s <- 1 - exp(-delta * h)

  pgle <- ifelse(s * g < n_patch,
                 yes = s * g,
                 no = n_patch)

  clnz <- prey * pgle

  ## extn: extinction rate
  if (length(mu) == 1) {
    mu <- rep(mu, 3)
  } else {
    if (length(mu) != 3) stop("error in mu")
  }

  extn <- mu[1] +
    mu[2] * (1 - (prey / max_prey)) +
    mu[3] * rho * u_length(lambda = lambda, size = size)

  ## equilibrium patch occupancy
  if (extn == 0 && clnz == 0)
    stop("colonization and extinction rates are both zero;
         equilibrium undefined.")

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
#' @param mu_s Numeric scalar or vector of spatial extinction rates.
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
                 mu_s = 1,
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

  ## check inputs
  l_par <- sapply(list(h, rsrc), function(x) length(x))
  if (any(l_par != 1))
    stop("Invalid length in h or rsrc.")

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
  v_mu_s <- to_v(mu_s, n = n_species)
  v_rho <- to_v(rho, n = n_species)
  u <- u_length(lambda = lambda, size = size)

  # run ode -----------------------------------------------------------------

  ## derivative
  derivr <- function(t, x, parms) {
    with(parms, {

      dx <- phi * (Mp %*% x + r) * x * (1 - x) -
        (mu0 +
           mu_p * (1 - (Mp %*% x) * inv_s_prey) +
           mu_c * Mc %*% x +
           mu_s * rho * u) * x

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
                mu_s = v_mu_s,
                rho = v_rho,
                u = u)

  x_init <- rep(x0, n_species)
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
                mu_s = 0.1,
                rho = 0.5) {

  # check input -------------------------------------------------------------

  absfwb <- abs(foodweb)
  if (!all(absfwb == t(absfwb)))
    stop("the input foodweb is invalid (abs(foodweb) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input foodweb is invalid (abs(foodweb) must be binary)")

  l_par <- sapply(list(h, rsrc), function(x) length(x))
  if (any(l_par != 1))
    stop("Invalid length in h or rsrc.")

  # transform input ---------------------------------------------------------

  ## foodweb: matrix, consumer-resource matrix. produce with ppm()
  fwb <- absfwb
  fwb[lower.tri(fwb)] <- 0

  ## p_hat: vector initialized with -1, equilibrium occupancy
  ## max_prey: vector, maximum number of prey items for consumer j
  p_hat <- rep(-1, ncol(fwb))
  max_prey <- colSums(fwb)

  ## constant terms, delta, rsrc, g, mu0, mu_p, mu_s, rho
  ## - create vectors with n-species elements
  n_sp <- unique(dim(foodweb))
  list_parms <- lapply(list(delta, g, mu0, mu_p, mu_s, rho),
                       FUN = to_v, n_sp)

  names(list_parms) <- c("delta",
                         "g",
                         "mu0",
                         "mu_p",
                         "mu_s",
                         "rho")


  # occupancies -------------------------------------------------------------

  ## sequential determination of equilibrium occupancies
  for (j in seq_len(ncol(fwb))) {

    if (max_prey[j] == 0) {
      ## basal species
      p_hat[j] <- with(list_parms,
                       p_base(lambda = lambda,
                              size = size,
                              h = h,
                              delta = delta[j],
                              rsrc = rsrc,
                              mu = c(mu0[j], mu_s[j]),
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
                              mu = c(mu0[j], mu_p[j], mu_s[j]),
                              rho = rho[j],
                              g = g[j])
                       )

    } # ifelse
  } # for j


  # food chain length -------------------------------------------------------

  ## report fcl as the maximum binary trophic position in the landscape
  if (any(p_hat > 0)) {
    ## at least one species persist

    ## index of persistent species
    ## subset the foodweb by persistent species
    index_p <- which(p_hat > 0)
    sub_fw <- as.matrix(fwb[index_p, index_p])
    tp <- rep(-1, ncol(sub_fw))

    ## index of basal species, number of basal, number of prey
    index_b <- which(colSums(sub_fw) == 0)
    n_b <- sum(colSums(sub_fw) == 0)
    n_p <- colSums(sub_fw)

    if (any(n_p > 0)) {
      ## tp = 1 if basal
      tp[index_b] <- 1

      for (i in (n_b + 1):length(tp)) {
        ## update tp recursively if consumers
        tp[i] <-  (drop(tp %*% sub_fw)[i]) / n_p[i] + 1
      }

    } else {

      tp[index_b] <- 1

    }

    fcl <- max(tp)
    attr(fcl, "tp") <- tp
  } else {
    ## no species persist
    fcl <- 0
  }

  attr(fcl, "p_hat") <- p_hat

  return(fcl)
}


#' Equilibrium food chain length (numerical)
#'
#' @inheritParams u_length
#' @inheritParams npom
#' @param foodweb Matrix. Binary food web matrix from \code{mcbrnet::ppm()}
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
                 mu_s = 1,
                 rho = 0.5,
                 x0 = 0.5,
                 n_timestep = 100,
                 interval = 0.01,
                 threshold = 1E-5) {

  ## numerical simulation
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
               mu_s = mu_s,
               rho = rho,
               x0 = x0,
               n_timestep = n_timestep,
               interval = interval,
               threshold = threshold)

  p_hat <- cout[nrow(cout), -1]

  ## report fcl as the maximum binary trophic position in the landscape
  fwb <- abs(foodweb)
  fwb[lower.tri(fwb)] <- 0

  if (any(p_hat > 0)) {
    ## at least one species persist

    ## index of persistent species
    ## subset the foodweb by persistent species
    index_p <- which(p_hat > 0)
    sub_fw <- as.matrix(fwb[index_p, index_p])
    tp <- rep(-1, ncol(sub_fw))

    ## index of basal species, number of basal, number of prey
    index_b <- which(colSums(sub_fw) == 0)
    n_b <- sum(colSums(sub_fw) == 0)
    n_p <- colSums(sub_fw)

    if (any(n_p > 0)) {
      ## tp = 1 if basal
      tp[index_b] <- 1

      for (i in (n_b + 1):length(tp)) {
        ## update tp recursively if consumers
        tp[i] <-  (drop(tp %*% sub_fw)[i]) / n_p[i] + 1
      }

    } else {

      tp[index_b] <- 1

    }

    fcl <- max(tp)
    attr(fcl, "tp") <- tp
  } else {
    ## no species persist
    fcl <- 0
  }

  attr(fcl, "p_hat") <- p_hat

  return(fcl)
}
