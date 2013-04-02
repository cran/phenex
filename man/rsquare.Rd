\name{rsquare}
\alias{rsquare}

\title{Squared Pearson correlation coefficient}
\description{Calculates the squared Pearson correlation coefficient of corrected and modelled timeserie.}
\usage{rsquare(x)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing corrected and modelled values. 
		If corrected timeserie is not available, the raw one out of \sQuote{values} will be used}
}
\value{The squared Pearson correlation coefficient as numeric value.}
\seealso{\code{\linkS4class{NDVI}}}
\author{Daniel Doktor, Maximilian Lange}
\examples{
	# load data
	data(avhrr)

	# create NDVI object, correct and model values
	ndvi.list <- modelNDVI(ndvi.values=avhrr.ndvi/10000, 
			year.int=1995, correction="bise", method="Growth", MARGIN=2, 
			doParallel=FALSE, slidingperiod=40)

	#plot
	plot(ndvi.list[[1]])

	# squared Pearson correlation coefficient
	rsquare(ndvi.list[[1]])
}