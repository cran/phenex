.initNDVI <-function() {
	setClass("OptionalInteger")
	setClassUnion("OptionalInteger", c("logical", "integer"))

	setClass("OptionalVector")
	setClassUnion("OptionalVector", c("NULL", "vector"))

	setClass("NDVI", representation(year="OptionalInteger", seasons="OptionalVector", values="vector", 					
		correctedValues="OptionalVector", modelledValues="OptionalVector"), 
		prototype(year=NA, seasons=NULL, values=rep(NA, 365), correctedValues=NULL, modelledValues=NULL))

	setGeneric("year", function(x) standardGeneric("year"))
	setMethod("year", "NDVI", function(x) x@year)	

	setGeneric("seasons", function(x) standardGeneric("seasons"))
	setMethod("seasons", "NDVI", function(x) x@seasons)

	setGeneric("values", function(x) standardGeneric("values"))
	setMethod("values", "NDVI", function(x) x@values)

	setGeneric("correctedValues", function(x) standardGeneric("correctedValues"))
	setMethod("correctedValues", "NDVI", function(x) x@correctedValues)

	setGeneric("modelledValues", function(x) standardGeneric("modelledValues"))
	setMethod("modelledValues", "NDVI", function(x) x@modelledValues)

	setGeneric("isLeapYear", function(x) standardGeneric("isLeapYear"))
	setMethod("isLeapYear", "NDVI", function(x) {
		leap <- FALSE
		y <- year(x)
		if ((y %% 4)==0) { leap <- TRUE }
		if ((y %% 100)==0) { leap <- FALSE }
		if ((y %% 400)==0) { leap <- TRUE }
		return(leap)
	})

	setValidity("NDVI", function(object) {
		if (!is.vector(values(object))){
			stop("'values' has to be a vector")
		}
		if (length(which(!is.na(values(object)))) > 0){
			if ( (min(values(object), na.rm=TRUE) < -1) | 
				(max(values(object), na.rm=TRUE) > 1)){
				stop("Values of 'values' have to be in interval (-1,1)")
			}
		}
		if ( !is.null(correctedValues(object)) & !is.vector(correctedValues(object)) ){
			stop("'correctedValues' has to be NULL or a vector")
		}
		if ( !is.null(modelledValues(object)) & !is.vector(modelledValues(object)) ){
			stop("'correctedValues' has to be NULL or a vector")
		}
		return(TRUE)
	})

	setGeneric("year<-", function(x, value) standardGeneric("year<-"))
	setReplaceMethod("year", "NDVI", function(x, value) {x@year <- value; validObject(x); x})

	setGeneric("seasons<-", function(x, value) standardGeneric("seasons<-"))
	setReplaceMethod("seasons", "NDVI", function(x, value) {x@seasons <- value; validObject(x); x})

	setGeneric("values<-", function(x, value) standardGeneric("values<-"))
	setReplaceMethod("values", "NDVI", function(x, value) {x@values <- value; validObject(x); x})

	setGeneric("runningAvg<-", function(x, value) standardGeneric("runningAvg<-"))
	setReplaceMethod("runningAvg", "NDVI", function(x, value) {x@runningAvg <- value; validObject(x); x})

	setGeneric("correctedValues<-", function(x, value) standardGeneric("correctedValues<-"))
	setReplaceMethod("correctedValues", "NDVI", function(x, value) {x@correctedValues <- value; validObject(x); x})

	setGeneric("modelledValues<-", function(x, value) standardGeneric("modelledValues<-"))
	setReplaceMethod("modelledValues", "NDVI", function(x, value) {x@modelledValues <- value; validObject(x); x})

	setGeneric("bise", function(x, slidingperiod, growthFactorThreshold, cycleValues) standardGeneric("bise"))
	setMethod("bise", "NDVI", function(x, slidingperiod, growthFactorThreshold, cycleValues) {
		if (missing("slidingperiod")){ slidingperiod <- 40 }
		if (missing("growthFactorThreshold")){ growthFactorThreshold <- 0.1 }
		if (missing("cycleValues")){ cycleValues <- TRUE }

		ndvi <- values(x)
		daysofyear <- length(ndvi)
		ndvi[which(ndvi <= 0)] <- NA

		if (cycleValues){ 
			days <- 3*daysofyear
			ndvi <- c(ndvi,ndvi,ndvi)
		} else {
			days <- daysofyear
		}
		
		corrected <- rep(NA,length=days)
		pos <- 1
		lastValidPos <- 0
		while (pos <= days){
			if (is.na(ndvi[pos])){
				pos <- pos+1
				next
			}

			if (lastValidPos==0){
				if ((ndvi[pos] > mean(ndvi, na.rm=TRUE)/5) && 
					(ndvi[pos] <= mean(ndvi, na.rm=TRUE))){
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
				}
				pos <- pos+1
				next
			}

			validIncrease <- (1+(growthFactorThreshold*(pos-lastValidPos)))*
						corrected[lastValidPos]

			validIncrease <- ifelse(validIncrease > 1, 1, validIncrease)

			if (ndvi[pos] >= corrected[lastValidPos]) {
				if ((ndvi[pos] <= validIncrease) || (ndvi[pos] < 0.2)) {
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
				}
				pos <- pos+1
			} else if (ndvi[pos] < corrected[lastValidPos]) {

				endOfPeriod <- pos+slidingperiod
				endOfPeriod <- ifelse(endOfPeriod>days,days,endOfPeriod)
				period <- pos:endOfPeriod
				nextValues <- ndvi[period]

				cloudValues <- nextValues-corrected[lastValidPos]
				cloudHop <- which(cloudValues>0)
				if (length(cloudHop)>0){
					pos <- period[cloudHop[1]]
					next
				}

				slopeThreshold <- ndvi[pos]+0.2*(corrected[lastValidPos]-ndvi[pos])
				values <- nextValues-slopeThreshold
				possibleValues <- which(values>0)
				if (length(possibleValues)>0){
					pos <- period[possibleValues[1]]
				} else {
					corrected[pos] <- ndvi[pos]
					lastValidPos <- pos
					pos <- pos+1
				}
			}
		}

		if (cycleValues){
			corrected <- corrected[(daysofyear+1):(2*daysofyear)]
		}

		correctedValues(x) <- corrected
		validObject(x)
		return(x)
	})

	setGeneric("runningAvg", function(x, window) standardGeneric("runningAvg"))
	setMethod("runningAvg", "NDVI", function(x, window) {
		if (missing("window")){ window <- 7 }
		ndvi <- values(x)
		days <- length(ndvi)

		res <- correctedndvi <- vector(mode="numeric",length=days);
		
		ndvi <- ifelse( is.na(ndvi), -1, ndvi)

		cndvi <- .C("runAVG", rdays = as.integer(days), 
				ndvi=as.numeric(ifelse(is.na(ndvi), -1, ndvi)),
				window=as.integer(7), cndvi=as.numeric(correctedndvi), 
				PACKAGE="phenex")$cndvi
		res <- ifelse(cndvi < 0, NA, cndvi)

		#evaluate res (not enough values)
		if (length(which(is.na(res) == FALSE)) < 5){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}
		if (res[order(res, decreasing=TRUE)[5]] < 0.1){
			correctedValues(x) <- as.vector(rep(NA, days))
			validObject(x)
			return(x)
		}

		correctedValues(x) <- as.vector(res)
		validObject(x)
		return(x)
	})

	setGeneric("modelValues", function(x, method, ...) standardGeneric("modelValues"))
	setMethod("modelValues", "NDVI", function(x, method, ...) {
		# ...: filter.threshold for fft
		# ...: asym for gauss
		# ...: window.sav, degree, smoothing for savgol
		# method can be 'LinIP', 'Spline', 'DSig', 'DSigC',
		# 'DLogistic', 'Gauss', 'Growth', 'FFT' or 'SavGol'
		if (tolower(method)!="linip" & tolower(method)!="spline" &
		tolower(method)!="dsig" & tolower(method)!="dsigc" &
		tolower(method)!="dlogistic" & tolower(method)!="gauss" &
		tolower(method)!="growth" & tolower(method)!="fft" &
		tolower(method)!="savgol" & tolower(method)!="gaussmix") {
			stop("method has to be one of the following: 
				'LinIP', 'Spline', 'DSig', 'DSigC', 
				'DLogistic', 'Gauss', 'GaussMix', 
				'Growth', 'FFT' or 'SavGol'")
		}
		arguments <- names(list(...))
		if (tolower(method)=="gauss"){
			if (is.na(match("asym",arguments))){
				asym <- FALSE
			} else { asym <- list(...)$asym }
		}
		if (tolower(method)=="gaussmix"){
			if (is.na(match("components",arguments))){
				components <- 3
			} else { components <- list(...)$components }
		}
		if (tolower(method)=="fft"){
			if (is.na(match("filter.threshold",arguments))){
				filter.threshold <- 3
			} else { filter.threshold <- list(...)$filter.threshold }
		}
		if (tolower(method)=="savgol"){
			if (is.na(match("window.sav",arguments))){
				window.sav <- 7
			} else { window.sav <- list(...)$window.sav }
			if (is.na(match("degree",arguments))){
				degree <- 2
			} else { degree <- list(...)$degree }
			if (is.na(match("smoothing",arguments))){
				smoothing <- 10
			} else { smoothing <- list(...)$smoothing }
		}
		
		
		if(is.null(correctedValues(x))){
			ndvi <- values(x)
		} else {
			ndvi <- correctedValues(x)
		}
		
		ndvi.mod.all <- rep(NA, length(ndvi))

		if (is.null(seasons(x))) { seasons(x) <- c(1,length(ndvi)) }

		if (tolower(method)=="linip" | tolower(method)=="spline" | 
			tolower(method)=="fft" | tolower(method)=="savgol"){
			if (length(which(!is.na(ndvi))) >= 5) {
				if (tolower(method)=="linip"){ ndvi.mod.all <- .linIP(ndvi) }
				if (tolower(method)=="spline"){ ndvi.mod.all <- .spline(ndvi) }
				if (tolower(method)=="fft"){ ndvi.mod.all <- .fftfilter(ndvi, filter.threshold) }
				if (tolower(method)=="savgol"){ ndvi.mod.all <- .savGol(ndvi, window.sav, degree, smoothing) }
			}
		} else {
			ndvi.all <- ndvi
			for (season in 1:(length(seasons(x))-1)){
				ndvi <- ndvi.all[seasons(x)[season]:seasons(x)[season+1]]
				ndvi.mod <- rep(NA, length(ndvi))
				if (length(which(!is.na(ndvi))) >= 5) {
					if (tolower(method)=="dsig"){ ndvi.mod <- .dSig(ndvi) }
					if (tolower(method)=="dsigc"){ ndvi.mod <-  .dSigC(ndvi) }
					if (tolower(method)=="dlogistic"){ ndvi.mod <- .dLogistic(ndvi) }
					if (tolower(method)=="gauss"){ ndvi.mod <- .gauss(ndvi, asym) }
					if (tolower(method)=="gaussmix"){ ndvi.mod <- .gaussMix(ndvi, components) }
					if (tolower(method)=="growth"){ ndvi.mod <- .growth(ndvi) }
			
				}
				ndvi.mod.all[seasons(x)[season]:seasons(x)[season+1]] <- ndvi.mod
			}
		}
		modelledValues(x) <- ndvi.mod.all
		validObject(x)
		return(x)
	})

	if (!isGeneric("plot"))
      		setGeneric("plot", function(x,y,...) standardGeneric("plot"))
		setMethod("plot", "NDVI", function(x,y=NULL,...) {
		plot(1:length(values(x)), values(x), type="p",
			ylim=c(0,1), xlab="Day of the Year",
			ylab="NDVI", col="black", ...)
		if(!is.null(seasons(x))){
			abline(v=seasons(x), 
				col="orange", lty=2, ...)
		}
		if(!is.null(correctedValues(x))){
			points(1:length(correctedValues(x)), correctedValues(x), 
				col="red", ...)
		}
		if (!is.null(modelledValues(x))){
			lines(1:length(modelledValues(x)), modelledValues(x), 
				col="blue", ...)
		}
	})

	setGeneric("phenoPhase", function(x, phase, method, threshold, n) standardGeneric("phenoPhase"))
	setMethod("phenoPhase", "NDVI", function(x, phase, method, threshold, n) {
		if (is.null(seasons(x))) { seasons(x) <- c(1,length(modelledValues(x))) }
		nrOfSeasons <- length(seasons(x))-1
		errorreturn <- list(mean=rep(NA, nrOfSeasons), sd=rep(NA, nrOfSeasons))
		if ( length(which(!is.na(modelledValues(x)))) < 5 ) { return(errorreturn) }
		if (missing("n")){ n <- 1000 }
		
		if(is.null(correctedValues(x))){
			observations <- values(x)
		} else {
			observations <- correctedValues(x)
		} 

		if (tolower(phase)!="max" & tolower(phase)!="min" & tolower(phase)!="maxval" & tolower(phase)!="minval" &
			tolower(phase)!="greenup" & tolower(phase)!="senescence"){
			stop("'phase' has to be 'min', 'minval', 'max', 'maxval', 'greenup' or 'senescence'")
		}

		if (tolower(phase)=="minval" | tolower(phase)=="maxval" | 
			tolower(phase)=="min" | tolower(phase)=="max"){
			vals <- c()
			sds <- c()
			doys <- c()
			for (season in 1:(length(seasons(x))-1)){
				days <- seasons(x)[season]:seasons(x)[season+1]
				model <- modelledValues(x)[days]
				obs <- observations[days]
				maxDOY <- order(model, decreasing=TRUE, na.last=TRUE)[1]
				if (tolower(phase)=="minval" | tolower(phase)=="min"){
					doy <- order(model[1:maxDOY], decreasing=FALSE, na.last=TRUE)[1]
				} else {
					doy <- maxDOY
				}

				if (tolower(phase)=="minval" | tolower(phase)=="maxval"){
					val <- model[doy]
					sdVal <- 0.02+0.02*val+sd(model-obs, na.rm=TRUE)
		
					vals <- c(vals, val)
					sds <- c(sds, sdVal)
				} else {
					doy <- seasons(x)[season]-1+doy
					doys <- c(doys, doy)
				}
			}
			if (tolower(phase)=="minval" | tolower(phase)=="maxval"){
				return(list(mean=vals, sd=sds))
			} else {
				return(list(mean=doys, sd=rep(NA, length(doys))))
			}
		}

		if (tolower(method)!="local" & tolower(method)!="global") {
			stop("'method' has to be 'local' or 'global'")
		}
		if (!is.numeric(threshold)){
			stop("'threshold' has to be a numeric value")
		}
		if (threshold < 0 | threshold > 1){
			stop("'threshold' has to be in interval (0,1)")
		}

		doys <- c()
		sds <- c()
		for (season in 1:(length(seasons(x))-1)){
			days <- seasons(x)[season]:seasons(x)[season+1]
			model <- modelledValues(x)[days]
			obs <- observations[days]

			if (tolower(phase)=="greenup"){ start <- 1 }
			if (tolower(phase)=="senescence"){ start <- length(model) }

			doy <- NA
			maxDoy <- order(model, decreasing=TRUE, na.last=TRUE)[1]
			if (method=="local"){				
				if (tolower(phase)=="senescence"){
					minmodel <- min(model[maxDoy:length(model)], na.rm=TRUE)
				} else {
					minmodel <- min(model[1:maxDoy], na.rm=TRUE)
				}
				maxmodel <- max(model, na.rm=TRUE)
				thresval <- ((maxmodel - minmodel) * threshold) + minmodel
			} else if (method=="global"){
				thresval <- threshold
			}

			if (n > 0){
				varThres <- 0.02+0.02*thresval + sd(model-obs, na.rm=TRUE)
				possibleThres <- rnorm(n=n, mean=thresval, sd=varThres)
			} else {
				possibleThres <- thresval
			}

			if (length(which(possibleThres > 1 | possibleThres < 0)) > 0){
				possibleThres <- possibleThres[-which(possibleThres > 1 | possibleThres < 0)]
			}

			modelhalf <- model[start:which(model == max(model, na.rm=TRUE))[1]]

			possibledoy <- c()
			for (thres in possibleThres){
				thresvec <- modelhalf-rep(thres, length(modelhalf))

				if (length(which(thresvec > 0)) > 0){
					doy <- which(thresvec > 0)[na.omit(match(which(thresvec <= 0)+1, which(thresvec > 0)))]
					if (tolower(phase)=="greenup"){
						doy <- doy[order(doy, decreasing=FALSE)[1]]
					} else if (tolower(phase)=="senescence"){ 
						doy <- doy[order(doy, decreasing=TRUE)[1]]
						doy <- length(model) + 1 - doy
						doy <- doy+1 #reverse correction
					} else {
						doy <- NA
					}
				} else {
					doy <- NA
				}
				possibledoy <- c(possibledoy, doy)
			}
			if (length(possibledoy)==0){
				doys <- c(doys, NA)
				sds <- c(sds, NA)
			} else if (length(which(!is.na(possibledoy)))>0){
				doys <- c(doys, days[round(mean(possibledoy, na.rm=TRUE),0)])
				sds <- c(sds, sd(possibledoy, na.rm=TRUE))
			} else {
				doys <- c(doys, NA)
				sds <- c(sds, NA)
			} 
		}
		return(list(mean=doys, sd=sds))
	})

	
	setGeneric("rsquare", function(x) standardGeneric("rsquare"))
	setMethod("rsquare", "NDVI", function(x) {
		if (is.null(modelledValues(x))){ return(NA) }
		if(is.null(correctedValues(x)) | length(correctedValues(x)) < 2){
			ndvi <- values(x)
		} else {
			ndvi <- correctedValues(x)
		}

		model <- modelledValues(x)

		check <- which(!is.na(ndvi) & !is.na(model))

		if ((length(check) > 0) && 
			(var(ndvi[check])!=0) && 
			(var(model[check])!=0)){
			r2 <- cor(ndvi[check], model[check])^2
		} else {
			r2 <- NA
		}
		return(r2)
	})

	setGeneric("integrateTimeserie", function(x, start, end, n) standardGeneric("integrateTimeserie"))
	setMethod("integrateTimeserie", "NDVI", function(x, start, end, n){
		if (missing("n")){ n <- 1000 }		
		if (is.null(seasons(x))) { seasons(x) <- c(1,length(modelledValues(x))) }
		nrOfSeasons <- length(seasons(x))-1
		errorreturn <-list(mean=rep(NA, nrOfSeasons), sd=rep(NA, nrOfSeasons))
		if(is.null(modelledValues(x))){ return(errorreturn) }
		if (length(which(!is.na(start$mean)))==0 | length(which(!is.na(end$mean)))==0){ 
			return(errorreturn) 
		}
		if ((n > 0) & (length(which(!is.na(start$sd)))==0 | length(which(!is.na(end$sd)))==0)){
			return(errorreturn) 
		}
		if (length(start$mean)!=length(end$mean) | length(start$mean)!=length(start$sd) | 
			length(end$sd)!=length(end$mean) ){ return(errorreturn) }

		model <- modelledValues(x)
		if(length(which(!is.na(model)))<2){ return(errorreturn) }
		modelfunc <- approxfun(x=(-length(model)+1):(2*length(model)), y=c(model,model,model))
	
		gsiviValue <- c()
		gsiviSd <- c()
		for (season in 1:length(start$mean)){
			values <- c()
			if (n > 0){
				starts <- rnorm(n=n, mean=start$mean, sd=start$sd)
				ends <- rnorm(n=n, mean=end$mean, sd=end$sd)
				nLength <- n
			} else {
				starts <- start$mean
				ends <- end$mean
				nLength <- 1

			}
			for (i in 1:nLength){		
				count <- 1
				success <- FALSE
				repeat {
					subdiv <- 100
					res <- try(integrate(f=modelfunc, lower=starts[i], 
							upper=ends[i], subdivisions=subdiv), 
						silent=TRUE)
					if (!inherits(res, "try-error")){
						success <- TRUE
						break;
					} else {
						subdiv <- subdiv*2
						count <- count+1
						if (count >= 20){
							break;
						}
					}
				}

				if (success){
					values <- c(values, res$value)
				} else {
					values <- NA
				}
			}
			if (length(which(!is.na(values)))==0){ 
				gsiviValue <- c(gsiviValue, NA)
				gsiviSd <- c(gsiviSd, NA)  
			} else {
				gsiviValue <- c(gsiviValue, mean(values, na.rm=TRUE))
				gsiviSd <- c(gsiviSd, sd(values, na.rm=TRUE))  
			}
		}		
		
		
		return(list(mean=gsiviValue, sd=gsiviSd))
	})

	setGeneric("detectSeasons", function(x, minValRange, ...) standardGeneric("detectSeasons"))
	setMethod("detectSeasons", "NDVI", function(x, minValRange, ...){
		if (missing("minValRange")){ minValRange <- 50 }
		if (is.null(values(x))){ return(x) }
		if (is.null(minValRange)){ return(x) }

		arguments <- names(list(...))
		if (is.na(match("slidingperiod",arguments))){
				slidingperiod <- 30
		} else { slidingperiod <- list(...)$slidingperiod }
		if (is.na(match("growthFactorThreshold",arguments))){
				growthFactorThreshold <- 0.1
		} else { growthFactorThreshold <- list(...)$growthFactorThreshold }
		if (is.na(match("cycleValues",arguments))){
				cycleValues <- TRUE
		} else { cycleValues <- list(...)$cycleValues }

		bisePoints <- correctedValues(bise(x, slidingperiod, growthFactorThreshold, cycleValues))

		fftvals <- rep(1, length(bisePoints))

		thres <- 100

		while (var(fftvals)==0){
			fftvals <- .fftfilter(bisePoints, filter.threshold=thres)
			thres <- thres-1
		}

		dFFT <- diff(fftvals)
		minima <- which(dFFT>0)[na.omit(match(which(dFFT<=0)+1,which(dFFT>0)))]
		maxima <- which(dFFT<=0)[na.omit(match(which(dFFT>0)+1,which(dFFT<=0)))]

		if (maxima[1] < minima[1]){ minima <- c(1, minima) }
		if (tail(maxima,n=1) > tail(minima,n=1)){ minima <- c(minima, length(bisePoints)) }

		bisePointsOrdered <- order(bisePoints, decreasing=FALSE)
		seasonStarts <- c()
		for (minimum in minima){
			availableMinima <- which(abs(bisePointsOrdered-minimum) < minValRange)
			if (length(availableMinima) > 0){
				bisePos <- bisePointsOrdered[availableMinima[1]]
				seasonStarts <- c(seasonStarts, bisePos)
			}
		}

		seasonStarts <- unique(seasonStarts)
		seasons(x) <- seasonStarts
		validObject(x)
		return(x)
	})
}
