\name{integrateTimeserie}
\alias{integrateTimeserie}

\title{Integration of timeserie}
\description{Extracts the integral of the vegetation index between start and end date.}
\usage{integrateTimeserie(x, start, end)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing modelled NDVI values.}
	\item{start}{Starting date for integration.} 
	\item{end}{End date for integration.} 
}
\value{An object of class \sQuote{integrate}, containing \dQuote{value}, 
	\dQuote{abs.error}, \dQuote{subdivisions}, \dQuote{method} and \dQuote{call}.}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{integrate}}}
\author{Daniel Doktor, Maximilian Lange}
\examples{
	# load data
	data(avhrr)

	# create NDVI object, correct and model NDVI values
	ndvi <- modelNDVI(ndvi.values=avhrr.ndvi/10000, year.int=1995, 
				correction="bise", method="LinIP", MARGIN=2, 
				doParallel=FALSE, slidingperiod=40)[[1]]

	# extract greenup DOY
	greenup <- phenoPhase(ndvi, phase="greenup", method="local", threshold=0.55)
	senesc <- phenoPhase(ndvi, phase="senescence", method="local", threshold=0.55)

	# extract green season integrated vegetation index
	gsivi <- integrateTimeserie(ndvi, start=greenup, end=senesc)
}
