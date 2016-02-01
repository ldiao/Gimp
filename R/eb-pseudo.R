#' @title empirical Bayes adjustment of sparse counts matrix
#'
#' @description
#' \code{eb.pseudo} performs empirical Bayes adjustment for sparse count matrix "counts",
#' where each sample is independently adjusted.
#'
#' @details
#' Essentially a wrapper funtion for \code{eb.pseudo.sample}, which performs empirical
#' Bayes adjustment for a feature counts vector for a single sample.
#'
#' @param counts Integer count matrix with features in rows and samples in columns
#'
#' @seealso \code{\link{eb.pseudo.sample}} for adjustment of a single sample. See
#'    \code{\link{sim.groups}} for a simple example.
#' @export
eb.pseudo <- function(counts) {
  s <- colSums(counts)
  X <- sapply(X   = 1:ncol(counts),
              FUN = function(i) {s[i] * eb.pseudo.sample(x = counts[, i])})
  colnames(X) <- colnames(counts)
  rownames(X) <- rownames(counts)
  X
}

#' @title empirical Bayes adjustment of sparse counts vector
#' @description Primarily a helper function for \code{eb.pseudo}, but can be used
#'    on its own to normalize a single sample
#' @param x Feature count vector
#' @seealso \code{\link{eb.pseudo}}
#' @return List with objects \code{z}, \code{psi}, and \code{mu}
#' @export
eb.pseudo.sample <- function(x) {
  x <- as.vector(x)
  s <- sum(x)
  fm <- fac.moment(x, return.scaled = TRUE)
  mu <- fm[1]
  sigma <- sqrt(fm[2] - mu^2)
  var.x <- mean((x / s) ^ 2) - (mean(x / s) ) ^ 2
  psi <- sigma / sqrt(var.x)
  z <- psi * x / s + (1 - psi) * mu
  z
}

#' @title calculates the first two factorial moments of a numeric vector
#' @description Helper function for \code{eb.pseudo.sample}
#' @param x Integer count matrix with features in rows and samples in columns
#' @param return.scaled Whether or not to return the scaled factorial moments
#' @seealso \code{\link{eb.pseudo.sample}}
fac.moment <- function(x, return.scaled = TRUE) {
  x <- sort(as.numeric(x))
  s <- sum(x)

  fm <- c(0, 0)
  fm[1] <- ifelse(return.scaled, mean(x) / s, mean(x) )

  x0 <- x
  x <- x * (x0 - 1)
  fm[2] <- ifelse(return.scaled, mean(x) / s ^ 2, mean(x))
  fm
}
