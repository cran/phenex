date2doy <- 
function(date) {
	#Calculates the Julian Day out of the integer date (YYMMDD)
	tmp <- 0
	doy <- .C("dateToDoy", as.integer(date), doy=as.integer(tmp), PACKAGE="phenex")$doy
	return(doy)
}
