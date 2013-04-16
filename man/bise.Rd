\name{bise}
\alias{bise}

\title{Best index slope extraction}
\description{Reduces noise in NDVI time-series.}
\usage{bise(x, slidingperiod)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing raw NDVI values.}
	\item{slidingperiod}{Sliding Period of the BISE-algorithm, default value is 40.}
}
\details{Knowledge regarding the phenological cycle in temperate climates and its temporal evolution 
	is used to detect and eliminate cloud contaminated observations. As the algorithm is searching 
	forward within daily NDVI observations over 1 year, decreases are only accepted if no higher 
	value is found within a so called sliding period. A period of 40 days proved best for our 
	study area but might have to be modified when study areas in different climates are 
	investigated.}
\value{An object of class \sQuote{NDVI} containing raw and corrected NDVI values.}
\seealso{\code{\link{modelNDVI}}, \code{\linkS4class{NDVI}}}
\author{Daniel Doktor, Maximilian Lange}
\references{
Viovy, N., Arino, O. and Belward, A.S. (1992). The Best Index Slope Extraction (BISE) - a method for reducing noise in NDVI time-series. \emph{International Journal of Remote Sensing}, \bold{13}, 1585-1590.
}
\examples{
	# load data
	data(avhrr)

	# create NDVI object
	ndvi <- new("NDVI", values=avhrr.ndvi/10000, year=as.integer(1995))

	# correct values (bise)
	ndvi.bise <- bise(ndvi, slidingperiod=40)

	#plot
	plot(ndvi.bise)
}