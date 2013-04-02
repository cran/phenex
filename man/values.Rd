\name{values}
\alias{values}

\title{Raw Value Accessor}
\description{Access to raw values of NDVI object.}
\usage{values(x)}
\arguments{
  \item{x}{An object of class \sQuote{NDVI} containing raw NDVI values.}
}
\value{Returns a vector containing raw NDVI values.}
\author{Lange, Maximilian and Doktor, Daniel}
\seealso{\code{\linkS4class{NDVI}}}
\examples{
	# load data 
	data(avhrr)

	# create NDVI object
	ndvi <- new("NDVI", values=avhrr.ndvi/10000, year=as.integer(1995))

	#extract raw values
	rawValues <- values(ndvi)
}

