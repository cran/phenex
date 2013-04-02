\name{modelValues}
\alias{modelValues}

\title{NDVI modelling}
\description{Models NDVI values.}
\usage{modelValues(x, method, ...)}
\arguments{
  \item{x}{An object of class \sQuote{NDVI} containing raw and/or corrected NDVI values.}
  \item{method}{Determines which model will be fitted to the corrected NDVI-timeserie (if the 
				corrected timeserie is not available, the raw one out of \sQuote{values} will be used). 

				\dQuote{\bold{LinIP}}: A linear interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1).
				Applied in \emph{Badeck et. al} (2004) and \emph{Doktor et. al} (2009).

				\dQuote{\bold{Spline}}: A spline interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1).

				\dQuote{\bold{DSig}}: Fits a double sigmoidal function to NDVI values (according to \emph{Zhang et. al} (2003)).

				\dQuote{\bold{DSigC}}: Fits another double sigmoidal function (own C implementation) to NDVI values.

				\dQuote{\bold{DLogistic}}: Fits a double logistic function after \emph{Fischer, Alberte} (1994) to NDVI values.

				\dQuote{\bold{Gauss}}: Fits a symmetric or asymmetric (determined by boolean parameter \sQuote{asym}) gaussian 
					function (own C implementation after \emph{Press, W.H.} (1992)) to NDVI values.

				\dQuote{\bold{Growth}}: Fits a growth model after \emph{Richter et. al} (1991) to NDVI values.

				\dQuote{\bold{FFT}}: Smoothes the corrected or raw NDVI values with fast fourier transfusion (implemented in R). 
					The smoothing intensity can be controlled with parameter \sQuote{filter.threshold} with default to 3.

				\dQuote{\bold{SavGol}}: Smoothes the corrected or raw NDVI values with a Savitzky-Golay 
					filter (own C implementation after \emph{Press, W.H.} (1992)). 
					The smoothing algorithm can be modified with parameters 
					\sQuote{window} (window size of filter, default to 7),
					\sQuote{degree} (degree of fitting polynomial, default to 2) 
					and \sQuote{smoothing} (repetition quantity, default to 10).}
	\item{...}{Other parameters passed to modelling function. 
					\sQuote{asym} for method \dQuote{Gauss},
					\sQuote{filter.threshold} for method \dQuote{FFT} and
					\sQuote{degree}, \sQuote{window} and \sQuote{smoothing} for method \dQuote{SavGol}.}

}
\details{Returns an object of type \sQuote{NDVI} containing raw data and/or corrected NDVI values and modelled NDVI values.}
\author{Lange, Maximilian and Doktor, Daniel}
\references{
Badeck, F.W., Bondeau, A., Boettcher, K., Doktor, D., Lucht, W., Schaber, J. and Sitch, S. (2004). Responses of spring phenology to climate change. \emph{New Phytologist}, \bold{162}, 295-309. 

Doktor, D., Bondeau, A., Koslowski, D. and Badeck, F.W. (2009). Influence of heterogeneous landscapes on computed green-up dates based on daily AVHRR NDVI observations. \emph{Remote Sensing of Environment}, \bold{113}, 2618-2632

Fischer, Alberte (1994). A Model for the Seasonal Variations of Vegetation Indices in Coarse Resolution Data and Its Inversion to Extract Crop Parameters. \emph{Remote Sensing of Environment}, \bold{48}, 220-230.

Press, W.H. (1992). Numerical recipes in C: The Art of Scientific Computing, vol. 1. Cambridge University Press, Cambridge, 2nd edn.

Richter, O., Spickermann, U. and Lenz, F. (1991). A new model for plant-growth. \emph{Gartenbauwissenschaft}, \bold{56}, 99-106.

Viovy, N., Arino, O. and Belward, A.S. (1992). The Best Index Slope Extraction (BISE) - a method for reducing noise in NDVI time-series. \emph{International Journal of Remote Sensing}, \bold{13}, 1585-1590.

Zhang, X.Y., Friedl, M.A., Schaaf, C.B., Strahler, A.H., Hodges, J.C.F., Gao, F., Reed, B.C. and Huete, A. (2003). Monitoring vegetation phenology using MODIS. \emph{Remote Sensing of Environment}, \bold{84}, 471-475.
}
\seealso{\code{\linkS4class{NDVI}}, \code{\link{modelNDVI}}}
\examples{
	# load data
	data(avhrr)

	# create NDVI object
	ndvi <- new("NDVI", values=avhrr.ndvi/10000, year=as.integer(1995))

	# correct values (bise)
	ndvi <- bise(ndvi, slidingperiod=40)

	#model values
	ndvi <- modelValues(ndvi, method="LinIP")

	# plot
	plot(ndvi)
}

