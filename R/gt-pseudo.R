#' @title Good-Turing adjustment of sparse counts matrix
#'
#' @description
#' \code{gt.pseudo} performs Good-Turing adjustment for sparse count matrix "counts",
#' where each sample is independently adjusted.
#'
#' @details
#' Essentially a wrapper funtion for \code{gt.pseudo.sample}, which performs
#' Good-Turing adjustment for a feature counts vector for a single sample.
#'
#' @param counts Integer count matrix with features in rows and samples in columns
#' @seealso \code{\link{gt.pseudo.sample}} for adjustment of a single sample. See
#'    \code{\link{sim.groups}} for a simple example.
#' @export
gt.pseudo <- function(counts) {
  apply(counts, 2, gt.pseudo.sample)
}

#' @title Good-Turing adjustment of sparse counts vector
#' @description Primarily a helper function for \code{gt.pseudo}, but can be used
#'    on its own to normalize a single sample
#' @param x Feature count vector
#' @seealso \code{\link{gt.pseudo}}
#' @export
gt.pseudo.sample <- function(x) {
  N <- sum(x)
  i0 <- which(x == 0)
  N0 <- length(i0)
  if (1 %in% x) {
    m0 <- pred.mass(x, k=0)
  } else{
    m0 <- pred.mass.smooth(x, k=0)
  }
  p0 <- m0 / N0
  c0 <- N * (p0 / (1 + p0))
  newx <- x / (1 + m0)
  newx[i0] <- c0
  newx
}

#' @title Good-Turing frequency estimator
#' @description Estimates the mass which should be assigned to all features with count
#'    k in sample x.
#' @param x vector of feature counts
#' @param k count value for which mass should be estimated
pred.mass <- function(x, k = 0) {
  (k + 1) * length(which(x == (k + 1))) / sum(x)
}

#' @title Simple Good-Turing frequency estimator, smoothed
#' @description Estimates the mass which should be assigned to all features with count
#'    k in sample x, smoothed by Gale's method (see "Good-Turing Smoothing Without Tears")
#' @param x vector of feature counts
#' @param k count value for which mass should be estimated
pred.mass.smooth <- function(x, k = 0) {
  f <- gt.fit.lm(x)
  pred.Nr <- exp(f$coefficients[1]) * (k + 1) ^ (f$coefficients[2])
  pred.Nr <- as.numeric(pred.Nr)
  (k + 1) * pred.Nr / sum(x)
}

#' @title Linear model for average transformed counts
#' @description Low-level function for smoothed Good-Turing frequency estimation. Fits
#'    a linear model to counts subjected to "averaging transformation" in order to obtain
#'    better frequency estimates for samples where no 1s are observed.
#' @param x vector of feature counts
#' @seealso \code{\link{avg.transform}}
gt.fit.lm <- function(x) {
  Z <- avg.transform(x, log.transform = TRUE)
  fit <- glm(formula = as.formula(Zr ~ r), data = Z[-1, ])
  fit
}

#' @title Averaging transformation
#' @description Averaging transformation for a vector of feature counts.
#' For more details, see "Good-Turing Smoothing Without Tears", by William A. Gale.
#' @param x vector of feature counts
#' @param log.transform whether or not log transform should be applied
avg.transform <- function(x, log.transform = TRUE) {
  M <- table(x)
  M <- matrix(c(as.numeric(names(M)), M), ncol = length(M), byrow = TRUE, dimnames = NULL)
  M <- t(M)

  newM <- M
  nr <- nrow(newM)
  newM[nr, 2] <- newM[nr, 2] / (0.5 * (M[nr, 1] - M[nr - 1, 1]))

  for (i in 2:(nrow(M) - 1)) {
    newM[i, 2] <- M[i, 2] / (0.5 * (M[i + 1, 1] - M[i - 1, 1]))
  }

  colnames(newM) <- c("r", "Zr")
  if (log.transform == TRUE)
    newM <- log(newM)

  data.frame(newM)
}

