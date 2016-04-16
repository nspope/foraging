library(RandomFields)
library(Rcpp)
####### start

.coordToCell <- function(coord, dim)
	dim[1]*(coord[1]-1) + coord[2]

.cellToCoord <- function(cell, dim)
	c( floor( (cell-1)/dim[1] ) + 1, ((cell-1) %% dim[1]) + 1)


.distance <- function(x, y)
	sqrt(sum((x-y)^2))

.trapCoords <- function(Traps, Landscape, snap=NULL){
	dim <- Landscape@dim
	buffer <- Landscape@buffer[1]
	if(Traps@setup == "transect"){
		start = buffer + 1
		end = dim[2] - buffer
		Traps@coords[,1] = round( seq(start, end, length.out = Traps@ntraps) )
		if( any(Traps@coords[,1] %% 1 != 0) ){
			stop("Number of traps does not match dimension")
		}
		Traps@coords[,2] = floor(median(1:dim[1]))
	} else if(Traps@setup == "grid") {
		if( sqrt(Traps@ntraps) %% 1 != 0 )
			stop("Number of traps must have an integer square root")
		sqT = sqrt(Traps@ntraps)
		startx = buffer
		endx = dim[2] - buffer
		starty = buffer
		endy = dim[1] - buffer
		seqx = round( seq(startx, endx, length.out = sqT) )
		seqy = round( seq(starty, endy, length.out = sqT) )
		if( any(c(seqx, seqy) %% 1 != 0 ) )
			stop("Number of traps does not match dimensions")
		Traps@coords = as.matrix(expand.grid(seqx, seqy))
	} else if(Traps@setup == "cross") {
		lenx = ceiling(Traps@ntraps / 2)
		leny = floor(Traps@ntraps / 2)
		startx = buffer
		endx = dim[2] - buffer
		vert <- floor(median(1:dim[1]))		
		seqx = round( seq(startx, endx, length.out = lenx) )
		
		midpoint = round(median(seqx))
		starty = buffer
		endy = dim[1] - buffer
		seqy = round( seq(starty, endy, length.out = leny + 1) )
		seqy = seqy[-c(leny/2 + 1)]
		if( any(c(seqx, seqy) %% 1 != 0 ) )
			stop("Number of traps does not match dimensions")
		Traps@coords = rbind( cbind(seqx, vert), cbind(midpoint, seqy))
	} else if(Traps@setup == "random") {
		yc <- sample(buffer:(dim[1]-buffer), Traps@ntraps, replace=TRUE)
		xc <- sample(buffer:(dim[2]-buffer), Traps@ntraps, replace=TRUE)
		Traps@coords = cbind( xc, yc)
	} else {
		stop("Trap type not recognized.")
	}
	if(!is.null(snap)){
		foraging <- matrix(Landscape@Foraging@data[,1], nrow=Landscape@dim[1], ncol=Landscape@dim[2])
		for(i in 1:Traps@ntraps){
			x <- Traps@coords[i,1]
			y <- Traps@coords[i,2]
			moveto <- arrayInd(which.max(foraging[(y-snap):(y+snap), (x-snap):(x+snap)]), .dim=c(snap*2+1, snap*2+1))
			adjust <- moveto - snap - 1
			Traps@coords[i,] <- Traps@coords[i,] + adjust[2:1]
		}
	}
	Traps@target = apply(Traps@coords, 1, .coordToCell, dim=dim)
	return(Traps)
}

importCpp <- function(infile, output_dir="lib", rebuild=FALSE){
    output_dir = ifelse(is.null(output_dir), ".", output_dir)
    dir.create(output_dir, recursive=T, showWarnings=FALSE)
    outfile = file.path(output_dir, paste0(infile, ".R"))

    if (!file.exists(outfile) || file.info(infile)$mtime > file.info(outfile)$mtime || rebuild){
        Rcpp::sourceCpp(infile, rebuild=rebuild)
        context = .Call("sourceCppContext", PACKAGE = "Rcpp",
            normalizePath(infile, winslash = "/"), code=NULL, 0L, .Platform)
        scriptfile = file.path(context$buildDirectory, context$rSourceFilename)
        content = readLines(scriptfile)
	ext = .Platform$dynlib.ext        
	m = regexpr(paste0("(?<=dyn.load\\(').*", ext), content[1], perl=TRUE)
        shlibfile = file.path(output_dir, paste0(infile, ext))
        shlibfile0 = regmatches(content[1], m)
        content[1] = sub(shlibfile0, shlibfile, content[1])

        f = file(outfile, "w+")
        writeLines(content, f)
        close(f)
        file.copy(shlibfile0, shlibfile, overwrite=TRUE)
    }else{
        source(outfile)
    }
    invisible(outfile)
}

library(Rcpp)
code = '
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".cellDistance")]]
List cellDistance(double beta, double theta, vec coords, uvec dim, vec kernel){
	int i,j;
	double d;
	long double mnd = 0;
	for(i=1;i<=dim[1];i++){
		for(j=1;j<=dim[0];j++){
			d = sqrt( pow(double(j) - coords[1], 2) + pow(double(i) - coords[0], 2) );
			kernel[(i-1)*dim[0] + j - 1] = exp(d*beta + theta*kernel[(i-1)*dim[0] + j - 1]);
			mnd += d * kernel[(i-1)*dim[0] + j - 1];
		}
	}
	List out(2);
	out[0] = double(mnd/sum(kernel));
	out[1] = kernel/sum(kernel);
	return( out );
}

// [[Rcpp::export(".colLik")]]
vec colLik(double beta, double theta, uvec dim, vec floral, mat traps, uvec trapcounts){
	int i,j,k;
	int K = trapcounts.n_elem;
	double d, ll, l, L;
	vec kernel(dim[0]*dim[1]);
	for(i=1;i<=dim[1];i++){
		for(j=1;j<=dim[0];j++){
			ll = 0;
			L = 0;
			for(k=0;k<K;k++){
				d = sqrt( pow(double(j) - traps(k,1), 2) + pow(double(i) - traps(k,0), 2) );
				l = beta * sqrt( pow(double(j) - traps(k,1), 2) + pow(double(i) - traps(k,0), 2) ) + theta*floral[k];
				ll += double(trapcounts[k])*l;
				L += exp(l);
			}
			kernel[(i-1)*dim[0] + j - 1] = ll - L;			
		}
	}
	return( kernel );
}
'
#sourceCpp(code=code, verbose=TRUE)
cat(code, file="mycode.cpp")
importCpp("mycode.cpp")

.makeBlacklist <- function(dim, buffersize){
	
	blacklist <- c()
	for(i in 1:dim[2])
		if( (i <= buffersize) | (i > (dim[2]-buffersize) ) )
			blacklist <- c(blacklist, ((i-1)*dim[1]+1):(i*dim[1]) )
		else
			blacklist <- c(blacklist, ((i-1)*dim[1]+1):((i-1)*dim[1]+buffersize), (i*dim[1]-buffersize+1):(i*dim[1]) )
	blacklist
}

setClass("Colony",
	## slots
	slots = c(
		coords = "numeric",
		bees = "numeric",
		beesStart = "numeric",
		kernel = "numeric",
		expDist = "numeric",
		par = "numeric"
	)
)

Colony <- function( nesting, foraging, dim, target, blacklist, par=c("lambda" = 100, "beta" = 0, "theta" = 0, "phi" = 0), coords = NULL, beesStart = NULL ){
	## add blacklist to buffer against here ...

	## location
	if(is.null(coords)){
		coords = sample( (1:length(nesting))[-blacklist], 1, prob=nesting[-blacklist])	
		coords = .cellToCoord(coords, dim)
	}
	#coords = coords + runif(2)

	## size of colony
	if(is.null(beesStart))
		#bees = beesStart = rpois( 1, par["lambda"] )
		bees = beesStart = sample(10:par["lambda"], 1 )
	else
		bees = beesStart

	## foraging kernel
	kernel = .cellDistance(par["beta"], par["theta"], coords, dim, foraging)
	expDist = kernel[[1]]
	kernel = kernel[[2]]
	kernel = kernel[target]
	
	new("Colony", coords = coords, bees = bees, beesStart = beesStart, kernel = kernel, expDist = expDist, par = par)
}

setClass("Trapping",
	## slots
	slots = c(
		coords = "matrix",
		ntraps = "numeric",
		radius = "matrix",
		trapcounts = "numeric",
		target = "numeric",
		trappedbees = "matrix",
		stoppingRule = "function",
		setup = "character"
	)
)

Trapping <- function(radius, ntraps, stoppingRule, setup){
	coords <- matrix(NA, nrow=ntraps, ncol=2)
	trapcounts <- rep(0, ntraps)
	target <- rep(0, ntraps)
	trappedbees <- matrix(nrow=0, ncol=3)
	new("Trapping", coords=coords, ntraps=ntraps, radius=radius, trapcounts=trapcounts, stoppingRule=stoppingRule, setup=setup, target=target, trappedbees=trappedbees)
}

setClass("Landscape",
	## slots
	slots = c(
		Foraging = "RFspatialGridDataFrame",
		Nesting = "RFspatialGridDataFrame",
		Traps = "Trapping",
		Colonies = "list",
		dim = "numeric",
		ncolonies = "numeric",
		totalbees = "numeric",
		fpar = "numeric",
		lpar = "numeric",
		blacklist = "numeric",
		buffer = "numeric"
	)
)

Landscape <- function(fpar, lpar, dim, buffer, ncolonies, Traps, snap = NULL, Nesting = NULL, Foraging = NULL){
	# landscapes
	library(RandomFields)
	if(is.null(Nesting)){
		Nesting = RFsimulate(model = RMfbm(lpar[1]), x=1:dim[1], y=1:dim[2], grid=TRUE)
		Nesting@data[,1] <- scale(Nesting@data[,1])
	}
	if(is.null(Foraging)){
		Foraging = RFsimulate(model = RMfbm(lpar[2]), x=1:dim[1], y=1:dim[2], grid=TRUE)
		Foraging@data[,1] <- scale(Foraging@data[,1])
	}
	Landscape <- new("Landscape", Foraging=Foraging, Nesting=Nesting, Traps=Traps, Colonies=list(), dim=dim, totalbees=0, ncolonies=ncolonies, fpar=fpar, lpar=lpar, buffer=buffer, blacklist=c(-1))

	# setup trapping grid
	Traps = .trapCoords(Traps, Landscape, snap = snap)

	# blacklist
	blacklist <- .makeBlacklist(dim, buffer[2])

	# colonies
	Colonies = list()
	totalbees = 0
	for(i in 1:ncolonies){
		Colonies[[i]] = Colony( exp(Nesting@data[,1]*fpar["phi"]), Foraging@data[,1], dim, Traps@target, blacklist, fpar)
		totalbees = totalbees + Colonies[[1]]@beesStart
	}

	Landscape@blacklist <- blacklist
	Landscape@Colonies <- Colonies
	Landscape@totalbees <- totalbees
	Landscape@Traps <- Traps
	Landscape
	#new("Landscape", Foraging=Foraging, Nesting=Nesting, Traps=Traps, Colonies=Colonies, dim=dim, totalbees=totalbees, ncolonies=ncolonies, fpar=fpar, lpar=lpar, buffer=buffer, blacklist=blacklist)
}

trapBee <- function(Landscape){
	trapKernel = rep(0, Landscape@Traps@ntraps)
	for(i in 1:Landscape@ncolonies)
		trapKernel = trapKernel + Landscape@Colonies[[i]]@kernel * Landscape@Colonies[[i]]@bees / Landscape@totalbees
	whichTrap = sample(1:Landscape@Traps@ntraps, 1, prob=trapKernel)

	colKernel = rep(0, Landscape@ncolonies)
	for(i in 1:Landscape@ncolonies)
		colKernel[i] = Landscape@Colonies[[i]]@kernel[whichTrap] * Landscape@Colonies[[i]]@bees / Landscape@totalbees
	whichCol = sample(1:Landscape@ncolonies, 1, prob = colKernel)

	Landscape@Colonies[[whichCol]]@bees = Landscape@Colonies[[whichCol]]@bees - 1
	Landscape@totalbees = Landscape@totalbees - 1
	Landscape@Traps@trapcounts[whichTrap] = Landscape@Traps@trapcounts[whichTrap] + 1
	Landscape@Traps@trappedbees = rbind(Landscape@Traps@trappedbees, c( "beeid"=sum(Landscape@Traps@trapcounts), "colony"=whichCol, "trap"=whichTrap ))
	return(Landscape)
}

runSim <- function(Landscape){
	trapping = TRUE
	while(trapping & Landscape@totalbees > 0){
		Landscape <- trapBee(Landscape)
		trapping = !(Landscape@Traps@stoppingRule(Landscape@Traps@trapcounts))
	}
	return(Landscape)
}

plotKernel <- function(Landscape, colNum, bins=NULL, fpar=NULL, traps=NULL){
	library(ggplot2)

	if(!is.null(fpar))
		Landscape@fpar <- fpar

	## foraging kernel
	kernel = .cellDistance(Landscape@fpar["beta"], Landscape@fpar["theta"], Landscape@Colonies[[colNum]]@coords, Landscape@dim, Landscape@Foraging@data[,1])
	kernel = kernel[[2]]
	
	crd <- expand.grid(y=1:Landscape@dim[1], x=1:Landscape@dim[2])
	crd$pr = kernel
	plt <- ggplot(crd, aes(x=x, y=y)) + geom_tile(aes(fill=pr, color=pr)) 
	if(!is.null(bins))
		plt = plt + geom_contour(aes(z=pr), color="black", bins=bins, alpha=0.3, size=0.25)
	if(!is.null(traps)){
		colnames(traps) <- c("x", "y")
		traps <- data.frame(traps)		
		traps$pr <- kernel[Landscape@Traps@target]
		plt = plt + geom_point(data=traps, aes(x=x, y=y, size=pr), pch=21, alpha = 0.5, fill="dodgerblue", col="black")
	}
	plt + theme_minimal() +
		annotate(geom="point", pch=23, x=Landscape@Colonies[[colNum]]@coords[1], y=Landscape@Colonies[[colNum]]@coords[2], fill="firebrick", col="black", size=3) +
#		scale_fill_gradient(low="gray20", high="white") + scale_color_gradient(low="gray20", high="white") + 
		scale_fill_gradient(low="white", high="gray20") + scale_color_gradient(low="white", high="gray20") + 
		theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position="none", panel.border=element_rect(fill=NA)) +
		scale_x_continuous(expand=c(0,0), breaks = seq(50,250,100), labels = seq(0.5,2.5,1)) + scale_y_continuous(expand=c(0,0), breaks = seq(50,250,100), labels = seq(0.5,2.5,1)) + xlab("Easting") + ylab("Northing") + theme(axis.title=element_blank())
}

plotLikelihood <- function(Landscape, trapcounts, fpar=NULL, bins=NULL, logt=FALSE){
	library(ggplot2)
	
	if(length(trapcounts) != Landscape@Traps@ntraps)
		stop("trapcounts is wrong length")

	if(!is.null(fpar))
		Landscape@fpar <- fpar

	crd <- expand.grid(y=seq(1,Landscape@dim[1]), x=seq(1,Landscape@dim[2]))

	#
	traps <- Landscape@Traps@coords
	colnames(traps) <- c("x", "y")
	traps <- data.frame(traps)
	traps$cnt <- trapcounts
	traps$flo <- Landscape@Foraging@data[Landscape@Traps@target,1]
	colfunc <- colorRampPalette(c("white", "gold"))#colorRampPalette(c("gold", "white", "blue"))
	cols <- data.frame(collis = colfunc(100), val=seq(min(traps$flo), max(traps$flo), length.out=100))
	traps$floCol <- cols$collis[sapply(traps$flo, function(x) which.min( abs( x - cols$val ) ) )]

	# calculate likelihood here
	crd$ll <- .colLik(Landscape@fpar["beta"], Landscape@fpar["theta"], Landscape@dim, Landscape@Foraging@data[Landscape@Traps@target,1], Landscape@Traps@coords, trapcounts)
	fudge <- max(crd$ll)
	crd$ll <- crd$ll - fudge
	if(!logt){
		crd$l <- exp(crd$ll)
	} else {
		crd$l <- crd$ll
	}
	p1 <- ggplot( crd, aes(x=x, y=y)) + geom_tile(aes(fill=l, color=l)) + theme_minimal()

	## add contours
	if(!is.null(bins)) 
		p1 <- p1 + geom_contour(aes(z=l), bins=bins, color="black", alpha=0.5, size=0.5)

#	p1 <- p1 + geom_point(data=traps, aes(x=x, y=y, size=cnt), pch=21, alpha=0.5, fill="dodgerblue", col="black")
	p1 <- p1 + geom_point(data=traps, aes(x=x, y=y, size=cnt), pch=21, alpha=0.5, fill=traps$floCol, col="black")
	p1 <- p1 + scale_fill_gradient(low="white", high="gray20") + scale_color_gradient(low="white", high="gray20")
	p1 <- p1 + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), legend.position="none", panel.border=element_rect(fill=NA), panel.grid=element_blank())
	p1 <- p1 + scale_x_continuous(expand=c(0,0), breaks = seq(50,250,100), labels = seq(0.5,2.5,1)) + scale_y_continuous(expand=c(0,0), breaks = seq(50,250,100), labels = seq(0.5,2.5,1)) 
	p1 <- p1 + xlab("Easting") + ylab("Northing") + theme(axis.title=element_blank())

	p1
}


plotForaging <- function(Landscape){
	## foraging landscape
	image(matrix(Landscape@Foraging@data[,1], nrow=Landscape@dim[1], ncol=Landscape@dim[2]))
}

resetTrap <- function(Landscape){
	Landscape@totalbees <- 0
	for(i in 1:Landscape@ncolonies){
		Landscape@Colonies[[i]]@bees <- Landscape@Colonies[[i]]@beesStart
		Landscape@totalbees = Landscape@totalbees + Landscape@Colonies[[i]]@beesStart
		Landscape@Traps@trappedbees <- matrix(nrow=0,ncol=3)
		Landscape@Traps@trapcounts <- rep(0, Landscape@Traps@ntraps)
	}
	return(Landscape)
}

changeTrap <- function(Landscape, Trap, snap=NULL){
	Traps <- .trapCoords(Trap, Landscape, snap=snap)	
	Landscape <- resetTrap(Landscape)
	Landscape@Traps = Traps
	return(Landscape)
}

recalcKernel <- function(Landscape, fpar){
	for(i in 1:Landscape@ncolonies){
		Landscape@Colonies[[i]] <- Colony( exp(Landscape@Nesting@data[,1]*fpar["phi"]), Landscape@Foraging@data[,1], Landscape@dim, Landscape@Traps@target, Landscape@blacklist, fpar, coords = Landscape@Colonies[[i]]@coords, beesStart = Landscape@Colonies[[i]]@beesStart)
	}
	Landscape <- resetTrap(Landscape)
	Landscape@fpar = fpar
	return(Landscape)
}

.makeResp <- function(landscape, dropCol = TRUE){
	dt <- as.data.frame(landscape@Traps@trappedbees)
	dt$trap <- factor(dt$trap, levels=1:landscape@Traps@ntraps)
	dt$colony <- factor(dt$colony, levels=1:landscape@ncolonies)	
	y <- with(dt, table(colony, trap))[]
	if(dropCol){
		y <- y[apply(y,1,sum)>0,]
	}
	y
}

runModel <- function(landscape, outputname, models, chains=3, iter=2000, settings=NA){
	library(rstan)

	## make data here
	myData = list()
#	myData$y = table(landscape@Traps@trappedbees[,"colony"], landscape@Traps@trappedbees[,"trap"] )[]
	myData$y = .makeResp(landscape)
	myData$trap = landscape@Traps@coords # trap coords
	myData$C = nrow(myData$y)
	myData$K = ncol(myData$y)
#	myData$trap_distance = as.matrix(dist(landscape@Traps@coords))
	myData$lowerbound = 0
	myData$upperbound = max(landscape@dim)
	myData$floral = as.vector(landscape@Foraging@data[landscape@Traps@target,1])
	myData$priorVa = 10
	myData$priorCo = 10
	
	## fit models
	output <- list()
	flags <- rep(NA, length(models))
	for(i in 1:length(models)){
		output[[i]] <- try( sampling(models[[i]], myData, chains=chains, iter=iter), silent=TRUE )
		if( class(output[[i]]) == "stanfit" ){
			if( all( summary(output[[i]])$summary[,"Rhat"] < 1.1 ) ){
				flags[i] <- 0
			} else {
				flags[i] <- 1
			}
		}
	}
	output[["Landscape"]] <- landscape
	output[["Data"]] <- myData
	output[["filename"]] <- outputname
	output[["settings"]] <- settings
	output[["flags"]] <- flags

	save(output, file=outputname)
}

landscapeParSim <- function(nsim){
	# simulate landscape
	lpars <- runif(2, 0.5, 1.5)
	fpars <- matrix(nrow=nsim, ncol=13)
	colnames(fpars) <- c("seed", "id", "ncolonies", "lambda", "beta", "theta", "phi", "shape", "dens", "snap", "toc", "total", "countsPer")
	# simulate within landscape
	seeds <- sample.int(10e8, nsim)
	for(i in 1:nsim){
		ncolonies <- sample(50:500, 1)
		lambda <- sample(50:400, 1)
		beta <- runif(1, -0.03, 0)
		theta <- runif(1, 0, 3)
		phi <- runif(1, 0, 2)

		shape <- sample( 1:4, 1) ## c("grid", "transect", "cross", "random")
		dens <- sample(2:6, 1)^2
		snap <- sample( 1:0, 1) ## c("yes", "no")

		total.or.counts <- sample(1:2, 1)
		total <- sample(50:1400, 1)
		countsPer <- sample(5:40, 1)
		fpars[i,] <- c("seed"=seeds[i], "id"=i, "ncolonies"=ncolonies, "lambda"=lambda, "beta"=beta, "theta"=theta, "phi"=phi, "shape"=shape, "dens"=dens, "snap"=snap, "toc"=total.or.counts, "total"=total, "countsPer"=countsPer)
	}
	while( any(duplicated(fpars[,"seed"])) ){
		whichdup <- which(duplicated(fpars[,"seed"]))
		fpars[whichdup,"seed"] <- fpars[whichdup,"seed"] + sample(1:50000,length(whichdup),replace=TRUE)
	}
	out <- list()
	out[[1]] <- lpars
	out[[2]] <- fpars
	out
}

makeModels <- function(modelfile){
	library(rstan)
	source(modelfile)
	## compile
	out <- list()	
	for(i in 1:length(modelcode))
		out[[i]] = stan_model(model_code = modelcode[[i]])
	return(out)
}

makeSeeds <- function(nseeds){
	set.seed(as.integer(Sys.time()))
	seeds <- sample.int(10e8, nseeds)
	if(any(duplicated(seeds)))
		return("fail")
	else
		return(seeds)
}
