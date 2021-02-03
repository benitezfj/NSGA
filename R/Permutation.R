#Test "Permutation"
#motsp
##City 1
city1 <- as.matrix(rmoo::kroA100)
kroA100dist <- dist(city1, method = "euclidean", diag = TRUE, upper = TRUE)
kroA100dist <- as.matrix(kroA100dist)

##City 2
city2 <- as.matrix(rmoo::kroB100)
kroB100dist <- dist(city2, method = "euclidean", diag = TRUE, upper = TRUE)
kroB100dist <- as.matrix(kroB100dist)

##City 3
city3 <- as.matrix(rmoo::kroC100)
kroC100dist <- dist(city3, method = "euclidean", diag = TRUE, upper = TRUE)
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
