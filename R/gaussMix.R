.gaussMix <- function(ndvi, components=2:4){
	days <- length(ndvi)
	ndvi.interpol <- .linIP(ndvi)
	time <- which(!is.na(ndvi))	
	ndvi <- ndvi[time]

	#calculate gaussian function
	Gm <- function(mu, sig, scal, base, time, days){
		erg <- ((scal / (sig * sqrt(2 * pi)))*
				exp(-0.5*((((time/days) - 
				(mu/days)) / sig)^2)))+base
		return(erg)
	} 

	Gm2 <- function(params, component, time, days){
		erg <- rep(0, length(time))
		for (comp in 1:component){
			erg <- erg + Gm(mu=params[(comp-1)*4+1], 
					sig=params[(comp-1)*4+2], 
					scal=params[(comp-1)*4+3], 
					base=params[(comp-1)*4+4], time, days)
		}
		return(erg)
	} 

	Gmx <- function(x, time, days, ndvi, component){
		erg <- sum((Gm2(params=x, component=component, time, days)-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}
	
	params <- list()
	lsq <- c()
	for (componentNr in 1:length(components)){
		component <- components[componentNr]
		lower <- rep(c(0,0,0,0), component)
		upper <- rep(c(days,2,2,1), component)

		NP <- 200*component

		model <- DEoptim(fn=Gmx, time=time, days=days, ndvi=ndvi, component=component,
			lower=lower, upper=upper,
			control=list(VTR=0, strategy=1, NP=NP, 
				itermax = 400, trace=FALSE, CR=0.9))

		lsq <- c(lsq, Gmx(x=model$optim$bestmem, time=1:days, days=days, ndvi=ndvi.interpol, component=component))
		params[[componentNr]] <- model$optim$bestmem
	}

	componentNr <- which(lsq==min(lsq,na.rm=TRUE)[1])
	component <- components[componentNr]
	parameterset <- params[[componentNr]]

	model.interpol <- Gm2(parameterset, component=component, time=1:days, days)
	
	return(model.interpol)
}
