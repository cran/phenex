analyzeBits <-
function(value,mode=1){
	if (is.na(value)){
		return(NA)
	}
	
	mask <- 0

	if (mode == 1)
		mask <- .C("getLandWater", number=as.integer(value), lw=as.integer(mask), PACKAGE="phenex")$lw;
	if (mode == 2)
		mask <- .C("getCloudmask", number=as.integer(value), cm=as.integer(mask), PACKAGE="phenex")$cm;
	if (mode == 3)
		mask <- .C("getDayNo", number=as.integer(value), day=as.integer(mask), PACKAGE="phenex")$day;
	
	return(mask)
}