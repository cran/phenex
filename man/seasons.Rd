\name{seasons}
\alias{seasons}

\title{Seasons Accessor}
\description{Access to detected seasons.}
\usage{seasons(x)}
\arguments{
  \item{x}{An object of class \sQuote{NDVI} containing raw NDVI values and detected seasons.}
}
\value{Returns a vector containing position of seasons.}
\author{Lange, Maximilian and Doktor, Daniel}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{detectSeasons}}}
\examples{
	# load data
	data(avhrr)

	# create NDVI object, correct and model values
	ndvi.list <- modelNDVI(ndvi.values=c(avhrr.ndvi/10000,avhrr.ndvi/10000), 
			year.int=1995, detectSeasons=TRUE,
			correction="bise", method="LinIP", MARGIN=2, 
			doParallel=FALSE, slidingperiod=40)
	ndvi <- ndvi.list[[1]]

	#get seasons
	seasondates <- seasons(ndvi)
}

