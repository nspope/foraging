library(Rcpp)
library(rstan)
library(plyr)
library(doMC)

doMC::registerDoMC(cores=8) # or however many cores you have access to

## input sseed, lid


source("simulate.R")

load("stanModels.RData") # NOT SHOWN--COMPILED MODELS IN STAN
set.seed(seedList[lid]) # NOT SHOWN--LIST OF RANDOM SEEDS FOR EACH SIMULATION

pars <- landscapeParSim(1000)

traptypes <- c("grid", "transect", "cross", "random")
buffers <- c(250, 100)
dims <- c(1000,1000)
snapset <- c(4)
radius <- matrix(c(0,0), nrow=1)

library(RandomFields)
RFoptions(maxGB=3)
nesting <- RFsimulate(model = RMfbm(pars[[1]][1]), x=1:dims[1], y=1:dims[2], grid=TRUE)
foraging <- RFsimulate(model = RMfbm(pars[[1]][2]), x=1:dims[1], y=1:dims[2], grid=TRUE)

nesting@data[,1] <- scale(nesting@data[,1])
foraging@data[,1] <- scale(foraging@data[,1])

execution_times <- adply( pars[[2]], 1, function(x){
	starttime <- Sys.time()
	set.seed(x["seed"])
	if( x["toc"]==1 ){
		stopping <- function(counts) sum(counts) >= x["total"]
		tcnt <- x["total"]
	} else {
		stopping <- function(counts) mean(counts) >= x["countsPer"]
		tcnt <- x["countsPer"] * x["dens"]
	}
	
	if( x["snap"]==1 ){
		snap <- snapset
	} else {
		snap <- NULL
	}

	Traps <- Trapping(radius, x["dens"], stopping, traptypes[x["shape"]])
	repeat{
		nL <- Landscape(fpar=c(x["lambda"],x["beta"],x["theta"],x["phi"]), 
				lpar=c(1,1), 
				dim=dims, 
				buffer=buffers, 
				ncolonies=x["ncolonies"], 
				Traps=Traps, 
				snap=snap, 
				Nesting=nesting, 
				Foraging=foraging)

		#### check that simulation is reasonable: number of bees in landscape exceeds number required by stopping rule
		if (tcnt < nL@totalbees)
			break
	}
	nL <- runSim(nL)
	outputname <- paste0("L", lid, "_", x["id"], ".RData")
	runModel(nL, outputname, models, chains=3, iter=2000, settings=x) ## need some error catching here

	stoptime <- Sys.time()
	return( c(x, "timing"=as.numeric(stoptime - starttime)) )
}, .parallel = TRUE) ## add parallelization

nmout <- paste0("L", lid, "_simInfo.tab")
write.table(execution_times, file=nmout, row.names=FALSE, sep="\t")
