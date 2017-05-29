.dLogistic <- 
function(ndvi){
	days <- length(ndvi)
	
	time <- which(!is.na(ndvi))	
	ndvi <- ndvi[time]

	Wt <- function(vb,ve,k,p,d,c,q, time) {
   		erg <- vb + (k/(1+exp(-c*(time-p)))) - ((k+vb-ve)/(1+exp(-d*(time-q))))
		return(erg)
  	}

	Wx <- function(x, time, ndvi){
		erg <- sum((Wt(vb=x[1], ve=x[2], k=x[3], p=x[4], d=x[5], c=x[6], q=x[7], time)-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}

	
	model <- DEoptim(fn=Wx, time=time, ndvi=ndvi,
		lower=c(0,0,0,1,0,0,1),
		upper=c(1,1,1,days,1,1,days),
		control=list(VTR=0, strategy=1, NP=200, itermax = 200, trace=FALSE, CR=0.9))

	vb <- model$optim$bestmem[1]
	ve <- model$optim$bestmem[2]
	k <- model$optim$bestmem[3]
	p <- model$optim$bestmem[4]
	d <- model$optim$bestmem[5]
	c <- model$optim$bestmem[6]
	q <- model$optim$bestmem[7]
	model.interpol <- Wt(vb,ve,k,p,d,c,q, time=1:days) 

	return(model.interpol)
}
