modelNDVI <- function(ndvi.values, year.int, multipleSeasons=FALSE, correction="bise", 
	method="LinIP", MARGIN=2, doParallel=FALSE, silent=TRUE, ...){
	if (!is.vector(ndvi.values) & !is.matrix(ndvi.values) & !is.array(ndvi.values)){
		stop("'ndvi.values' has to be of type 'vector', 'matrix' or 'array'")
	}
	if (tolower(correction)!="none" & tolower(correction)!="bise" &	tolower(correction)!="ravg"){
		stop("'correction' has to be 'none', 'bise' or 'ravg'")
	}
	if (tolower(method)!="linip" & tolower(method)!="spline" &
		tolower(method)!="dsig" & tolower(method)!="dsigc" &
		tolower(method)!="dlogistic" & tolower(method)!="gauss" &
		tolower(method)!="growth" & tolower(method)!="fft" &
		tolower(method)!="savgol" & tolower(method)!="none" & 
		tolower(method)!="gaussmix") {
			stop("method has to be one of the following: 
				'none', 'LinIP', 'Spline', 'DSig', 
				'DSigC', 'DLogistic', 'Gauss', 'GaussMix', 
				'Growth', 'FFT' or 'SavGol'")
	}
	if (!is.logical(doParallel)){
		stop("'doParallel' has to be of type 'logical'")
	}

	if (!is.logical(multipleSeasons)){
		stop("'multipleSeasons' has to be of type 'logical'")
	}
	
	if (is.null(dim(ndvi.values))){
		ndvi.values <- as.matrix(ndvi.values)
	}
	if ((length(dim(ndvi.values))-1) != length(MARGIN)){
		stop("'MARGIN' has to be of length 'dim(ndvi.values)-1'")
	}


	dnames <- dimnames(ndvi.values)
	if (is.character(MARGIN)) {
       	if (is.null(dnames))
			stop("'X' must have named dimnames")
        	MARGIN <- match(MARGIN, dnames)
       	if (any(is.na(MARGIN)))
			stop("not all elements of 'MARGIN' are names of dimensions")
	}

	dims <- dim(ndvi.values)
	dimnrs <- seq_len(length(dim(ndvi.values)))

	m.days <- dimnrs[-MARGIN]
    	m.it <- dimnrs[MARGIN]

	ndvi.values <- aperm(ndvi.values, c(m.days, m.it))
	values.length <- prod(dims[m.it])	

	arguments <- names(list(...))
	if (is.na(match("slidingperiod",arguments))){
		slidingperiod <- 40
	} else { slidingperiod <- list(...)$slidingperiod }
	if (is.na(match("growthFactorThreshold",arguments))){
		growthFactorThreshold <- 0.1
	} else { growthFactorThreshold <- list(...)$growthFactorThreshold }
	if (is.na(match("cycleValues",arguments))){
		cycleValues <- TRUE
	} else { cycleValues <- list(...)$cycleValues }
	if (is.na(match("window.ravg",arguments))){
		window.ravg <- 7
	} else { window.ravg <- list(...)$window.ravg }

	if (doParallel){
		# check if parallel backend is available
		if(!foreach::getDoParRegistered()) {
			if (!silent){ cat("No parallel backend detected! Problem will be solved sequential.\n",sep="") }
			foreach::registerDoSEQ()
		} else {
			if (!silent){ cat("Parallel backend detected.\n",sep="") }
		}
		
		ndvi.list <- foreach(i=1:values.length, .inorder=TRUE) %dopar% {
			position <- ((i-1)*dims[m.days]+1):(i*dims[m.days])
			ndvi.vec <- ndvi.values[position]
			ndvi.vec[which(ndvi.vec > 1 | ndvi.vec < -1)] <- NA
			ndvi <- new("NDVI", year=as.integer(year.int), values=ndvi.vec)
			if (multipleSeasons){ ndvi <- detectSeasons(ndvi) }
				
			if (correction=="bise"){ ndvi <- bise(ndvi, slidingperiod, 
							growthFactorThreshold, cycleValues) }
			if (correction=="ravg"){ ndvi <- runningAvg(ndvi, window.ravg) }

			if (method != "none"){ ndvi <- modelValues(ndvi, method=method, ...) }
			return(ndvi)
		}
	} else {
		ndvi.list <- list()
		for (i in seq_len(values.length)){
			position <- ((i-1)*dims[m.days]+1):(i*dims[m.days])
			ndvi.vec <- ndvi.values[position]
			ndvi.vec[which(ndvi.vec > 1 | ndvi.vec < -1)] <- NA
			ndvi <- new("NDVI", year=as.integer(year.int), values=ndvi.vec)
			if (multipleSeasons){ ndvi <- detectSeasons(ndvi) }

			if (correction=="bise"){ ndvi <- bise(ndvi, slidingperiod, 
							growthFactorThreshold, cycleValues) }
			if (correction=="ravg"){ ndvi <- runningAvg(ndvi, window.ravg) }

			if (method != "none"){ ndvi <- modelValues(ndvi, method=method, ...) }
			ndvi.list[[length(ndvi.list)+1]] <- ndvi
		}
	}

	return(ndvi.list)
}
