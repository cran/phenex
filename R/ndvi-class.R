.initNDVI <-function() {
	setClass("OptionalInteger")
	setClassUnion("OptionalInteger", c("logical", "integer"))

	setClass("OptionalVector")
	setClassUnion("OptionalVector", c("NULL", "vector"))

	setClass("NDVI", representation(year="OptionalInteger", values="vector", correctedValues="OptionalVector", modelledValues="OptionalVector"), 
		prototype(year=NA, values=rep(NA, 365), correctedValues=NULL, modelledValues=NULL))

	setGeneric("year", function(x) standardGeneric("year"))
	setMethod("year", "NDVI", function(x) x@year)

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
		if (!is.na(year(object)) & !is.integer(year(object))){
			stop("Year has to be 'integer' or 'NA'")
		}
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

	setGeneric("values<-", function(x, value) standardGeneric("values<-"))
	setReplaceMethod("values", "NDVI", function(x, value) {x@values <- value; validObject(x); x})

	setGeneric("runningAvg<-", function(x, value) standardGeneric("runningAvg<-"))
	setReplaceMethod("runningAvg", "NDVI", function(x, value) {x@runningAvg <- value; validObject(x); x})

	setGeneric("correctedValues<-", function(x, value) standardGeneric("correctedValues<-"))
	setReplaceMethod("correctedValues", "NDVI", function(x, value) {x@correctedValues <- value; validObject(x); x})

	setGeneric("modelledValues<-", function(x, value) standardGeneric("modelledValues<-"))
	setReplaceMethod("modelledValues", "NDVI", function(x, value) {x@modelledValues <- value; validObject(x); x})

	setGeneric("bise", function(x, slidingperiod) standardGeneric("bise"))
	setMethod("bise", "NDVI", function(x, slidingperiod) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("slidingperiod")){ slidingperiod <- 40 }
		days <- ifelse(isLeapYear(x), 366, 365)

		res <- correctedndvi <- vector(mode="numeric",length=days);

		ndvi <- values(x)
		ndvi <- ifelse( is.na(ndvi), 0, ndvi)

		res <- .C("bise", rdays=as.integer(days), ndvi=as.numeric(ndvi), 
			rslidperiod=as.integer(slidingperiod), cndvi=as.numeric(correctedndvi), 
			PACKAGE="phenex")$cndvi

		res <- ifelse ( res <= 0, NA, res)

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

		peaks <- c()
		threshold <- 1.5
		check <- which(is.na(res) == FALSE)
		for (i in 1:length(check)){
			if (i==1) {
				if (res[check[i]] > (threshold * mean(res[c(check[length(check)], check[i+1])]))){
					if (res[check[i]] > 0.3) {
						peaks <- c(peaks,check[i])
					}
				}
			} else {
				if (i==length(check)){
					if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[1])]))){
						if (res[check[i]] > 0.3) {
							peaks <- c(peaks,check[i])
						}
					}
				} else {
					if (res[check[i]] > (threshold * mean(res[c(check[i-1], check[i+1])]))){
						if (res[check[i]] > 0.3) {
							peaks <- c(peaks,check[i])
						}
					}
				}
			}
		}
		res[peaks] <- NA
		correctedValues(x) <- as.vector(res)
		validObject(x)
		return(x)
	})

	setGeneric("runningAvg", function(x, window) standardGeneric("runningAvg"))
	setMethod("runningAvg", "NDVI", function(x, window) {
		if (is.na(year(x))){ stop("'year' has to be set") }
		if (missing("window")){ window <- 7 }
		days <- ifelse(isLeapYear(x), 366, 365)

		res <- correctedndvi <- vector(mode="numeric",length=days);

		ndvi <- values(x)
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
		tolower(method)!="savgol") {
			stop("method has to be one of the following: 
				'LinIP', 'Spline', 'DSig', 'DSigC', 
				'DLogistic', 'Gauss', 'Growth', 
				'FFT' or 'SavGol'")
		}

		if (tolower(method)=="gauss"){
			if (!exists("asym")){ asym <- FALSE }
		}
		if (tolower(method)=="fft"){
			if (!exists("filter.threshold")){ filter.threshold <- 3 }
		}
		if (tolower(method)=="savgol"){
			if (!exists("window.sav")){ window.sav <- 7 }
			if (!exists("degree")){ degree <- 2 }
			if (!exists("smoothing")){ smoothing <- 10 }
		}

		if (is.na(year(x))){ stop("'year' has to be set") }
		
		ndvi.mod <- rep(NA, ifelse(isLeapYear(x),366,365))
		if(is.null(correctedValues(x))){
			ndvi <- values(x)
		} else {
			ndvi <- correctedValues(x)
		}
		
		if (length(which(!is.na(ndvi))) >= 5) {
			if (tolower(method)=="linip"){ ndvi.mod <- phenex:::.linIP(ndvi) }
			if (tolower(method)=="spline"){ ndvi.mod <- phenex:::.spline(ndvi) }
			if (tolower(method)=="dsig"){ ndvi.mod <- phenex:::.dSig(ndvi) }
			if (tolower(method)=="dsigc"){ ndvi.mod <-  phenex:::.dSigC(ndvi) }
			if (tolower(method)=="dlogistic"){ ndvi.mod <- phenex:::.dLogistic(ndvi) }
			if (tolower(method)=="gauss"){ ndvi.mod <- phenex:::.gauss(ndvi, asym) }
			if (tolower(method)=="growth"){ ndvi.mod <- phenex:::.growth(ndvi) }
			if (tolower(method)=="fft"){ ndvi.mod <- phenex:::.fftfilter(ndvi, filter.threshold) }
			if (tolower(method)=="savgol"){ ndvi.mod <- phenex:::.savGol(ndvi, window.sav, degree, smoothing) }
		}
		modelledValues(x) <- ndvi.mod
		validObject(x)
		return(x)
	})

	if (!isGeneric("plot"))
      		setGeneric("plot", function(x,y,...) standardGeneric("plot"))
	setMethod("plot", "NDVI", function(x,y=NULL,...) {
		plot(1:ifelse(isLeapYear(x),366,365), values(x), type="p",
			xlim=c(0,365), ylim=c(0,1), xlab="Day of the Year",
			ylab="NDVI", col="black", ...)
		if(!is.null(correctedValues(x))){
			points(1:ifelse(isLeapYear(x),366,365), correctedValues(x), 
				col="red", ...)
		}
		if (!is.null(modelledValues(x))){
			lines(1:ifelse(isLeapYear(x),366,365), modelledValues(x), 
				col="blue", ...)
		}
	})

	setGeneric("phenoPhase", function(x, phase, method, threshold) standardGeneric("phenoPhase"))
	setMethod("phenoPhase", "NDVI", function(x, phase, method, threshold) {
		if ( length(which(!is.na(modelledValues(x)))) < 5 ) { return(NA) }
		if (tolower(phase)!="max" & tolower(phase)!="min" & 
			tolower(phase)!="greenup" & tolower(phase)!="senescence"){
			stop("'phase' has to be 'min', 'max', 'greenup' or 'senescence'")
		}
		
		if (tolower(phase)=="max"){
			return(which(modelledValues(x)==max(modelledValues(x), na.rm=TRUE))[1])
		}
		if (tolower(phase)=="min"){
			return(which(modelledValues(x)==min(modelledValues(x), na.rm=TRUE))[1])
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

		if (tolower(phase)=="greenup"){ start <- 1 }
		if (tolower(phase)=="senescence"){ start <- length(modelledValues(x)) }

		doy <- NA
		model <- modelledValues(x) 

		if (method=="local"){
			minmodel <- min(model, na.rm=TRUE)
			maxmodel <- max(model, na.rm=TRUE)
			threshold <- ((maxmodel - minmodel) * threshold) + minmodel
		}

		modelhalf <- model[start:which(model == max(model, na.rm=TRUE))[1]]
		thresvec <- modelhalf-rep(threshold, length(modelhalf))
		
		doy <-  which(abs(thresvec) == min(abs(thresvec), na.rm=TRUE))
		doy <- doy[order(doy, decreasing=FALSE)[1]]
	
		if (tolower(phase)=="senescence"){ doy <- length(modelledValues(x)) - doy }	

		return(doy)
	})

	
	setGeneric("rsquare", function(x) standardGeneric("rsquare"))
	setMethod("rsquare", "NDVI", function(x) {
		if (is.null(modelledValues(x))){ return(NA) }
		if(is.null(correctedValues(x))){
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
}
