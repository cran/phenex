\name{integrateTimeserie}
\alias{integrateTimeserie}

\title{Integration of timeserie}
\description{Extracts the integral of the vegetation index between start and end date.}
\usage{integrateTimeserie(x, start, end, n)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing modelled NDVI values.}
	\item{start}{A list containing the starting date(s) for integration as \sQuote{mean} and its standard deviation(s) as \sQuote{sd}. 
			Use a list with multiple entries (as vector) for \sQuote{mean} and \sQuote{sd} if the NDVI object contains multiple seasons. } 
	\item{end}{A list containing the end date(s) for integration as \sQuote{mean} and its standard deviation(s) as \sQuote{sd}. 
		Use a list with multiple entries (as vector) for \sQuote{mean} and \sQuote{sd} if the NDVI object contains multiple seasons.} 
	\item{n}{The number \sQuote{n} of normal distributed values to create around start and end date.}
}
\value{A list containing the integral(s) as \sQuote{mean} and a standard devation(s) \sQuote{sd}.}
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
	greenup <- phenoPhase(ndvi, phase="greenup", method="local", threshold=0.55, n=1000)
	senesc <- phenoPhase(ndvi, phase="senescence", method="local", threshold=0.55, n=1000)

	# extract green season integrated vegetation index
	gsivi <- integrateTimeserie(ndvi, start=greenup, end=senesc, n=1000)
}
