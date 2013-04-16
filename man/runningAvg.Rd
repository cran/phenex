\name{runningAvg}
\alias{runningAvg}

\title{Running Average}
\description{Reduces noise in NDVI time-series through running averaging.}
\usage{runningAvg(x, window)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing raw NDVI values.}
	\item{window}{Window size of the running averaging algorithm. Default is 7.}
}
\value{An object of class \sQuote{NDVI} containing raw and corrected NDVI values.}
\seealso{\code{\link{modelNDVI}},\code{\linkS4class{NDVI}}}
\author{Daniel Doktor, Maximilian Lange}
\examples{
	# load data
	data(avhrr)

	# create NDVI object
	ndvi <- new("NDVI", values=avhrr.ndvi/10000, year=as.integer(1995))

	# correct values (bise)
	ndvi.ravg <- runningAvg(ndvi, window=7)

	#plot
	plot(ndvi.ravg)
}