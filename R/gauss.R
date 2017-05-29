.gauss <-
function(ndvi, asym=FALSE){
	if (!is.logical(asym)){ stop("'asym' should be of type 'logical'") }
	
	days <- length(ndvi)
	time <- which(!is.na(ndvi))	
	ndvi <- ndvi[time]

	#calculate gaussian function
	if (asym){
		Ga <- function(mu, sig1, sig2, scal, base, time, days){
			mu <- round(mu)
			sig <- rep(NA, days)
			sig <- c(rep(sig1, length(1:mu)), rep(sig2, length((mu+1):days)))
			sig <- sig[time]

			erg <- ((scal / (sig * sqrt(2 * pi)))*
					exp(-0.5*((((time/days) - 
					(mu/days)) / sig)^2)))+base

			if ( length(which(time<=mu)) > 0 & length(which(time>mu)) > 0){
				firstmax <- max(erg[which(time <= mu)], na.rm=TRUE)
				secondmax <-max(erg[which(time > mu)], na.rm=TRUE)

				if ((firstmax > secondmax | secondmax > 1) & firstmax < 1){
					erg[which(time > mu)] <- erg[which(time > mu)]*firstmax/secondmax
				} else {
					erg[which(time <= mu)] <- erg[which(time <= mu)]*secondmax/firstmax
				}
			}
			return(erg)
		} 

		Gax <- function(x, time, days, ndvi){
			erg <- sum((Ga(mu=x[1], sig1=x[2], sig2=x[3], scal=x[4], base=x[5], time, days)-ndvi)^2)
			return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
		}
		
		lower <- c(0,0,0,0,0)
		upper <- c(days,1,1,2,1)

		model <- DEoptim(fn=Gax, time=time, days=days, ndvi=ndvi,
		lower=lower, upper=upper,
		control=list(VTR=0, strategy=1, NP=200, 
			itermax = 200, trace=FALSE, CR=0.9))

	} else {
		G <- function(mu, sig, scal, base, time, days){
			erg <- ((scal / (sig * sqrt(2 * pi)))*
					exp(-0.5*((((time/days) - 
					(mu/days)) / sig)^2)))+base
			return(erg)
		} 

		Gx <- function(x, time, days, ndvi){
			erg <- sum((G(mu=x[1], sig=x[2], scal=x[3], base=x[4], time, days)-ndvi)^2)
			return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
		}

		lower <- c(0,0,0,0)
		upper <- c(days,1,2,1)

		model <- DEoptim(fn=Gx, time=time, days=days, ndvi=ndvi,
		lower=lower, upper=upper,
		control=list(VTR=0, strategy=1, NP=200, 
			itermax = 200, trace=FALSE, CR=0.9))

	}


	if (asym){
		mu <- model$optim$bestmem[1]
		sig1 <- model$optim$bestmem[2]
		sig2 <- model$optim$bestmem[3]
		scal <- model$optim$bestmem[4]
		base <- model$optim$bestmem[5]
	
		model.interpol <- Ga(mu, sig1, sig2, scal, base, time=1:days, days)
	} else {
		mu <- model$optim$bestmem[1]
		sig <- model$optim$bestmem[2]
		scal <- model$optim$bestmem[3]
		base <- model$optim$bestmem[4]

		model.interpol <- G(mu, sig, scal, base, time=1:days, days)
	}

	return(model.interpol)
}
