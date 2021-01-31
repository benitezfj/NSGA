#Three Objectives
DTLZ1 <- function (x, nobj = 3) {
  if (is.null(dim(x))) {
    x <- matrix(x, 1)
  }
  n <- ncol(x)
  y <- matrix(x[, 1:(nobj - 1)], nrow(x))
  z <- matrix(x[, nobj:n], nrow(x))
  g <- 100 * (n - nobj + 1 + rowSums((z - 0.5)^2 - cos(20 * 
      pi * (z - 0.5))))
  tmp <- t(apply(y, 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  tmp2 <- cbind(1, t(apply(1 - y, 1, rev)))
  f <- tmp * tmp2 * 0.5 * (1 + g)
  return(f)
}

DTLZ2  <- function (x, nobj = 3) {
  if (is.null(dim(x))) {
    x <- matrix(x, 1)
  }
  n <- ncol(x)
  y <- matrix(x[, 1:(nobj - 1)], nrow(x))
  z <- matrix(x[, nobj:n], nrow(x))
  g <- rowSums((z - 0.5)^2)
  tmp <- t(apply(cos(y * pi/2), 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
  f <- tmp * tmp2 * (1 + g)
}

DTLZ3 <-function (x, nobj = 3) {
  if (is.null(dim(x))) {
    x <- matrix(x, 1)
  }
  n <- ncol(x)
  y <- matrix(x[, 1:(nobj - 1)], nrow(x))
  z <- matrix(x[, nobj:n], nrow(x))
  g <- 100 * (n - nobj + 1 + rowSums((z - 0.5)^2 - cos(20 * 
      pi * (z - 0.5))))
  tmp <- t(apply(cos(y * pi/2), 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
  f <- tmp * tmp2 * (1 + g)
}

DTLZ7 <- function (x, nobj = 3) {
  if (is.null(dim(x))) {
    x <- matrix(x, 1)
  }
  n <- ncol(x)
  y <- matrix(x[, 1:(nobj - 1)], nrow(x))
  z <- matrix(x[, nobj:n], nrow(x))
  g <- 1 + 9 * rowSums(z/(1:(n - nobj + 1)))
  tmp <- cbind(y, 1)
  tmp2 <- cbind(matrix(1, nrow(x), nobj - 1), 
    (1 + g) * (nobj - rowSums(y/(1 + g) * (1 + sin(3 * pi * y)))))
  f <- tmp * tmp2
}