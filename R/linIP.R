.linIP <- function(ndvi){
	days <- length(ndvi)
	f <- approxfun(x=1:(3*days), y=c(ndvi,ndvi,ndvi))
	ndvi.interpol <- f((days+1):(2*days))
	return(ndvi.interpol)
}
