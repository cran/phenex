\name{correctedValues}
\alias{correctedValues}

\title{Corrected Value Accessor}
\description{Access to corrected values of NDVI object.}
\usage{correctedValues(x)}
\arguments{
  \item{x}{An object of class \sQuote{NDVI} containing raw and corrected NDVI values.}
}
\value{Returns a vector containing corrected NDVI values.}
\author{Lange, Maximilian and Doktor, Daniel}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{modelNDVI}}, \code{\link{bise}}, \code{\link{runningAvg}}}
\examples{
	# load data
	data(avhrr)

	# create NDVI object, correct and model values
	ndvi.list <- modelNDVI(ndvi.values=avhrr.ndvi/10000, year.int=1995, 
			correction="bise", method="LinIP", MARGIN=2, 
			doParallel=FALSE, slidingperiod=40)
	ndvi <- ndvi.list[[1]]

	#get modelled values
	biseValues <- correctedValues(ndvi)
}

