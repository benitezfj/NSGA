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


#Test "Permutation"
#motsp
##City 1
kroA100 <- read.table("kroA100.txt", header = TRUE)
kroA100 <- as.matrix(kroA100[,2:3])
kroA100dist <- dist(kroA100, method = "euclidean", diag = TRUE, upper = TRUE)
kroA100dist <- as.matrix(kroA100dist)

##City 2
kroB100 <- read.table("kroB100.txt", header = TRUE)
kroB100 <- as.matrix(kroB100[,2:3])
kroB100dist <- dist(kroB100, method = "euclidean", diag = TRUE, upper = TRUE)
kroB100dist <- as.matrix(kroB100dist)

##City 3
kroC100 <- read.table("kroC100.txt", header = TRUE)
kroC100 <- as.matrix(kroC100[,2:3])
kroC100dist <- dist(kroC100, method = "euclidean", diag = TRUE, upper = TRUE)
kroC100dist <- as.matrix(kroC100dist)

distArray <- array(NA, dim=c(nrow(kroA100dist), ncol(kroA100dist), 3))
distArray[,,1] <- kroA100dist
distArray[,,2] <- kroB100dist
distArray[,,3] <- kroC100dist

#tours <- sapply("nearest_insertion", FUN = function(m) solve_TSP(kroA100dist, method = m), simplify = FALSE)
#Lust T. and Teghem J.
tourLength <- function(tour, distArray) {
  cost1 <- distArray[,,1]
  cost2 <- distArray[,,2]
  cost3 <- distArray[,,3]
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  cbind(sum(cost1[route]), sum(cost2[route]), sum(cost3[route]))
}

#Keller. and Goodchild
tourLength <- function(tour, distArray) {
  distCities <- distArray[,1:100]
  prizeCities <- as.vector(distArray[,101])
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  cbind(sum(distCities[route]), -sum(prizeCities[tour]))
}

motspFitness <- function(tour, ...) 1/tourLength(tour, ...)

kroABC100 <- nsga3(type = "permutation", 
  fitness = motspFitness, 
  lower = 1, upper = 100, 
  seed = 123, 
  maxiter = 3000, 
  popSize = 50, 
  monitor = TRUE, summary = FALSE,
  n_partitions = 9, 
  pcrossover = 0.5, pmutation = 0.2, 
  distArray = distArray, nObj = 3)

kroAB100 <- nsga2(type = "permutation", 
  fitness = motspFitness, 
  lower = 1, upper = 100, 
  seed = 123, 
  maxiter = 3000, 
  popSize = 50, 
  monitor = TRUE, summary = FALSE,
  pcrossover = 0.5, pmutation = 0.2, 
  distArray = distArray, nObj = 2)

kroAB100 <- nsga3(type = "permutation", 
  fitness = motspFitness, 
  lower = 1, upper = 100, 
  seed = 123, 
  maxiter = 3000, 
  popSize = 50, n_partitions = 50,
  monitor = TRUE, summary = FALSE,
  pcrossover = 0.5, pmutation = 0.2, 
  distArray = distArray, nObj = 2)
