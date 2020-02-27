#Creamos una nueva clase de la union de numeric, logical y matrix
setClassUnion("numberOrNAOrMatrix", members = c("numeric", "logical", "matrix"))

#Modificado para NSGA-II
setClass(Class = "nsga",
         representation(call = "language",
                        type = "character",
                        lower = "numberOrNAOrMatrix",
                        upper = "numberOrNAOrMatrix",
                        nBits = "numberOrNAOrMatrix",
                        names = "character",
                        popSize = "numeric",
                        front = "numberOrNAOrMatrix",
                        f = "list",
                        iter = "numeric",
                        run = "numeric",
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = "numberOrNAOrMatrix",
                        elitism = "numeric",
                        pcrossover = "numeric",
                        pmutation = "numberOrNAOrMatrix",
                        optim = "logical",
                        dumFitness = "numberOrNAOrMatrix", #NSGA-I
                        dShare = "numeric", #NSGA-I
                        deltaDummy = "numeric", #NSGA-I
                        crowdingDistance = "numberOrNAOrMatrix", #NSGA-II
                        fitness = "numberOrNAOrMatrix",
                        summary = "matrix",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix"))
#,algorithm = "character"

###Se utilizará como bandera para saber cual version del NSGA utilizar.
