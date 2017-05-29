\name{detectSeasons}
\alias{detectSeasons}

\title{Season Detection}
\description{Detects seasons in timeseries.}
\usage{detectSeasons(x, minValRange, ...)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing raw NDVI values.}
	\item{minValRange}{Range in which to search for lowest bise value around detected season start.} 
	\item{...}{Optional parameters passed to internal bise call.} 
}
\value{An object of class \sQuote{NDVI}, containing \dQuote{seasons}.}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{seasons}}, \code{\link{bise}}}
\author{Daniel Doktor, Maximilian Lange}
\examples{
	# load data
	data(avhrr)

	# generate ndvi object
	ndvi <- new("NDVI", values=rep(avhrr.ndvi/10000,5), year=NA)

	# extract seasons
	ndvi <- detectSeasons(ndvi)
	seasons(ndvi)
}
