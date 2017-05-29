.dSig <-
function(ndvi){
	days <- length(ndvi)
	
	time <- which(!is.na(ndvi))	
	ndvi <- ndvi[time]

	DS <- function(pos1, width1, pos2, width2, time){
		erg <- 0.5*(tanh((time-pos1)/width1)-tanh((time-pos2)/width2))
		return(erg)
	}

	DSx <- function(x, time, ndvi){
		erg <- sum((DS(pos1=x[1],width1=x[3],pos2=x[2],width2=x[4],time)-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}

	model <- DEoptim(fn=DSx, time=time, ndvi=ndvi,
		lower=rep(0,4), upper=rep(days,4),
		control=list(VTR=0, strategy=1, NP=100, 
			itermax = 200, trace=FALSE, CR=0.9))

	
	pos1 <- model$optim$bestmem[1]
	pos2 <- model$optim$bestmem[2]
	width1 <- model$optim$bestmem[3]
	width2 <- model$optim$bestmem[4]  

	model.interpol <- DS(pos1, width1, pos2, width2, time=1:days) 

	return(model.interpol)
}
