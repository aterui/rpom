#' Utility: convert to vector format
#'
#' @param x value(s)
#' @param n number of replicates
#'
#' @author Akira Terui, \email{hanabi0111@gmail.com}
#'
#' @export

to_v <- function(x, n) {

  if (length(x) == 1) {

    v_x <- rep(x, n)

  } else {

    if (n != length(x))
      stop("incorrect input in one or more of the parameters")

    v_x <- x
  }

  return(v_x)
}
