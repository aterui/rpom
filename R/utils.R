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
#' @param w Binary food web matrix.
#' @param occupancy Equilibrium occupancies of constituent species
#' @param weight Logical.
#'  If \code{TRUE}, maximum trophic position is weighted by relative occupancies.
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

maxtp <- function(w,
                  occupancy,
                  weight = TRUE) {

  # check input -------------------------------------------------------------

  absfwb <- abs(w)

  if (!all(absfwb == t(absfwb)))
    stop("the input w is invalid (abs(w) must be symmetric)")

  if (any(!(absfwb %in% c(0, 1))))
    stop("the input w is invalid (abs(w) must be binary)")

  if (any(length(occupancy) != dim(w)))
    stop("the input w or occupancy has invalid dimension.")

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
    v_id_b <- which(colSums(m_sub_fwb) == 0)
    n_b <- sum(colSums(m_sub_fwb) == 0)
    v_n_prey <- colSums(m_sub_fwb)
    v_sum_o <- drop(v_o %*% m_sub_fwb)

    if (any(v_n_prey > 0)) {
      ## v_tp = 1 for basal
      v_tp[v_id_b] <- 1

      if (weight) {

        ## update v_tp recursively for consumers
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
