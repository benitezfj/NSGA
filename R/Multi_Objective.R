#Test Functions

#Test "Real-Valued"
#Two Objectives
ZDT1 <- function (x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  n <- ncol(x)
  g <- 1 + rowSums(x[, 2:n, drop = FALSE]) * 9/(n - 1)
  return(cbind(x[, 1], g * (1 - sqrt(x[, 1]/g))))
} 

ZDT2 <- function (x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  n <- ncol(x)
  g <- 1 + rowSums(x[, 2:n, drop = FALSE]) * 9/(n - 1)
  return(cbind(x[, 1], g * (1 - (x[, 1]/g)^2)))
}

ZDT3 <- function (x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  n <- ncol(x)
  g <- 1 + rowSums(x[, 2:n, drop = FALSE]) * 9/(n - 1)
  return(cbind(x[, 1], 
               g * (1 - sqrt(x[, 1]/g) - x[, 1]/g * sin(10 * pi * x[, 1]))))
}

ZDT4 <- function (x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  n <- ncol(x)
  g <- 1 + 10 * (n - 1) + rowSums((x[, 2:n, drop = FALSE] * 
      10 - 5)^2 - 10 * cos(4 * pi * (x[, 2:n, drop = FALSE] * 
          10 - 5)))
  return(cbind(x[, 1], g * (1 - sqrt(x[, 1]/g))))
}

#Test "Binary"
ZDT5 <- function (x, m = 10, n = 5, normal = TRUE) {
  if (length(x)%%5 != 0 && length(x) < 35){
    stop("Bit vector's length must contain at least 35.")
  }
  x1 <- x[1:30]
  xm <- x[31:length(x)]
  
  g <- 0
  vu <- function(v){
    if (v < 5){
      x <- 2 + v
    }else if (v == 5){
      x <- 1
    }
    return(x)
  }
  
  for (i in seq(m)) {
    x <- sum(xm[((n*(i-1))+1):(n*i)])
    v <- vu(x)
    g <- g + v
  }
  f1 <- 1 + sum(x1 == 1)
  f2 <- g * (1 / f1)
  
  normalize <- function(x, x_min, x_max) {
    # calculate the denominator
    denom <- x_max - x_min
    # we can not divide by zero -> plus small epsilon
    denom <- denom + 1e-30
    # normalize the actual values
    N <- (x - x_min) / denom
    return(N)
  }
  
  if (normal == TRUE){
    f1 <- normalize(f1, 1, 30)
    f2 <- normalize(f2, (m-1) * 1/30, (m-1))
  }
  
  return(cbind(f1,f2))
}

ZDT6 <- function (x) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }
  n <- ncol(x)
  f1 <- 1 - exp(-4 * x[, 1]) * (sin(6 * pi * x[, 1]))^6
  g <- 1 + 9 * (1/(n - 1) * rowSums(x[, 2:n, drop = FALSE]))^(0.25)
  return(cbind(f1, g * (1 - (f1/g)^2)))
}
