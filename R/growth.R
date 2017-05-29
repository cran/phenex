.growth <- 
function(ndvi){
	days <- length(ndvi)
	ndvi.interpol <- .linIP(ndvi)
	
	time <- which(!is.na(ndvi))	
	ndvi <- ndvi[time]

	Wt <- function(W, a, r, p, my, time){
		p <- log(abs(p))
		erg <- (W * ((a+1)^(r/exp(p))) * exp(((-1) * my) * time)) / 
			((1+a*exp(-exp(p)*time))^(r/exp(p)))
		return(erg)
	}

	Wx <- function(x, time, ndvi){
		erg <- sum((Wt(W=x[1], a=x[2], r=x[3], p=x[4], my=x[5], time)-ndvi)^2)
		return(ifelse((is.infinite(erg)||is.nan(erg)),1e50,erg))
	}
	

	model <- DEoptim(fn=Wx, time=time, ndvi=ndvi,
		lower=c(0,0,0,0,0), upper=c(0.8,1e4,0.3,0.1,0.5),
		control=list(VTR=0, strategy=1, NP=200, 
			itermax = 200, trace=FALSE, CR=0.9))

	
	W <- model$optim$bestmem[1]
	a <- model$optim$bestmem[2]
	r <- model$optim$bestmem[3]
	p <- model$optim$bestmem[4]
	my <- model$optim$bestmem[5]

	model.interpol <- Wt(W,a,r,p,my,time=1:days)

	return(model.interpol)
}
