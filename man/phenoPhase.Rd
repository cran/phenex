\name{phenoPhase}
\alias{phenoPhase}

\title{Phenological Phase Extraction}
\description{Extracts phenological metrics based on modelled NDVI values.}
\usage{phenoPhase(x, phase, method, threshold, n)}
\arguments{
	\item{x}{An object of class \sQuote{NDVI} containing modelled NDVI values.}
	\item{phase}{Determines which phase will be extracted:
			
			\dQuote{\bold{max}}: Day of the year with highest NDVI value is returned.

			\dQuote{\bold{maxval}}: The highest modelled NDVI value and its standard deviation is returned.

			\dQuote{\bold{min}}: Day of the year (before day with maximum NDVI value) 
						with lowest NDVI value is returned. 
		
			\dQuote{\bold{minval}}: The lowest modelled NDVI value and its standard deviation is returned.

			\dQuote{\bold{greenup}}: Day of the year at which greenup takes place and its standard deviation is returned. 
				This day refers to the point where the function of modelled NDVI values 
				exceeds a certain threshold.
			
			\dQuote{\bold{senescence}}: Day of the year at which senescence takes place and its standard deviation is returned. 
				This day refers to the point where the function of modelled NDVI values exceeds 
				a certain threshold.}
	\item{method}{Determines whether a global or local threshold is used for
					greenup and senescence extraction.

			\dQuote{\bold{global}} threshold: The day of the year is returned, 
				where NDVI values are first equal or higher as the value of \sQuote{threshold}.
				If the threshold is higher than the values of the timeseries, \sQuote{-1} will be returned.
					
			\dQuote{\bold{local}} threshold: The day of the year is returned, for which NDVI values first reach the value of
			\sQuote{threshold} (interpreted as percentage) between lowest and highest NDVI value of timeserie. 
					The lowest NDVI value is chosen depending on phase selected. For \dQuote{greenup}, 
					the lowest value before day of maximum NDVI value is used. For \dQuote{senescence}, 
					the lowest value after day of maximum NDVI value is used.}

	\item{threshold}{Threshold for local or global greenup/senescence-extraction method. 
		Should have numerical value between \sQuote{0} and \sQuote{1}. Global thresholds 
		refer to a fix NDVI value independent of actual NDVI magnitudes. 
		Local thresholds are based on (max - min) * e.g. 0.55. Hence, the respective 
		NDVI value will vary dependent on e.g. land-cover.}
	
	\item{n}{The number \sQuote{n} of normal distributed values around the threshold for estimation of \sQuote{sd}. The normal distribution uses the threshold as mean and a standaed deviation consisting of 			satellite error and standard deviation of fit residuals. The satellite error is considered as 0.02+0.02*value.}
}
\value{A list containing the julian day of the year at which the phenologial phase occurs as list entry \sQuote{mean} and its standard deviation as list entry \sQuote{sd}. The list contains vectors in \sQuote{mean} and \sQuote{sd} if multiple seasons are available in \sQuote{NDVI} object.}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{detectSeasons}}, \code{\link{seasons}}}
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
}
