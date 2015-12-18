\name{modelNDVI}
\alias{modelNDVI}

\title{Main function for NDVI correction and modelling}
\description{Fits a suite of functions/models to raw NDVI or selected NDVI observations (after BISE).}
\usage{modelNDVI(ndvi.values, year.int, correction="bise", 
	method="LinIP", MARGIN=2, doParallel=FALSE, silent=TRUE, ...)}
\arguments{
  \item{ndvi.values}{A vector or matrix containing NDVI observations with values in the interval (-1,1).}
  \item{year.int}{Observation year}
  \item{correction}{A character string determining which method will be used for correction of NDVI values. 
			Should be either:

			\dQuote{\bold{none}}: no correction will be performed.
			
			\dQuote{\bold{ravg}}: Smoothing with running average. Default window size is 7 
				and can be modified by parameter \sQuote{window.ravg}. See also \code{\link{runningAvg}}.

			\dQuote{\bold{bise1}}: Best index slope extraction after \emph{Viovy et. al} (1992).
			Can be modified with parameter \sQuote{slidingperiod}. Default to 40, optimised for 
			the area of Central Europe. This might has to be modified 
			when vegetation dynamics of other climates/habitats are analysed. 
			See also \code{\link{bise1}} and \code{\link{bise2}}.

			\dQuote{\bold{bise2}}: Best index slope extraction after \emph{Viovy et. al} (1992).
			Second interpretation of bise algorithm.
			Can be modified with parameter \sQuote{slidingperiod}. Default to 40, optimised for 
			the area of Central Europe. This might has to be modified 
			when vegetation dynamics of other climates/habitats are analysed. 
			An maximum increase threshold is defined via parameter \sQuote{growthFactorThreshold} with 
			default to 0.1 (10 percent increase per day allowed).
			The timeserie is cycled via parameter \sQuote{cycleValues}, which is a 
			boolean value determining whether the end of the ndvi timeserie 
			is combined with its beginning or not (default value is true). If false, gaps or low ndvi 
			values at the beginning of timeserie influence the result.
			See also \code{\link{bise1}} and \code{\link{bise2}}.
			}
  \item{method}{Determines which model will be fitted to the corrected NDVI-timeserie (if the 
				corrected timeserie is not available, the raw one out of \sQuote{values} will be used). 

				\dQuote{\bold{LinIP}}: A linear interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1). 
				Applied in \emph{Badeck et. al} (2004) and \emph{Doktor et. al} (2009).

				\dQuote{\bold{Spline}}: A spline interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1).

				\dQuote{\bold{DSig}}: Fits a double sigmoidal function to NDVI values (according to Zhang et. al (2003)).

				\dQuote{\bold{DSigC}}: Fits another double sigmoidal function (own C implementation) to NDVI values.

				\dQuote{\bold{DLogistic}}: Fits a double logistic function after \emph{Fischer, Alberte} (1994) to NDVI values.

				\dQuote{\bold{Gauss}}: Fits a symmetric or asymmetric (determined by boolean parameter \sQuote{asym}) gaussian 
					function to NDVI values (own C implementation after \emph{Press, W.H.} (1992)).

				\dQuote{\bold{Growth}}: Fits a plant growth model after \emph{Richter et. al} (1991) to NDVI values.

				\dQuote{\bold{FFT}}: Smoothes the corrected or raw NDVI values with fast fourier transformation (implemented in R). 
					The smoothing intensity can be controlled with parameter \sQuote{filter.threshold} with default to 3.

				\dQuote{\bold{SavGol}}: Smoothes the corrected or raw NDVI values with a Savitzky-Golay 
					filter (own C implementation after \emph{Press, W.H.} (1992)). 
					The smoothing algorithm can be modified with parameters 
					\sQuote{window.sav} (window size of filter, default to 7),
					\sQuote{degree} (degree of fitting polynomial, default to 2) 
					and \sQuote{smoothing} (repetition quantity, default to 10).}
	\item{MARGIN}{A vector giving the subscripts which the function will be applied over.  
			E.g., for a matrix \sQuote{1} indicates rows, \sQuote{2} indicates columns, 
			\sQuote{c(1, 2)} indicates rows and columns. Where \sQuote{X} has named dimnames, 
			it can be a character vector selecting dimension names.

			Should be \sQuote{2} if \sQuote{ndvi.values} is a vector instead of a matrix/array.}
	\item{doParallel}{This method uses 'foreach'. If a parallel backend is registered, setting \sQuote{do.parallel} to \sQuote{TRUE} enables parallel processing.}
	\item{silent}{A boolean flag determining whether debug information is shown.}
	\item{...}{Other parameters passed to correction or modelling function. 
			These are: \sQuote{slidingperiod} for correction \dQuote{bise},
					\sQuote{window.ravg} for correction \dQuote{ravg},
					\sQuote{asym} for method \dQuote{Gauss},
					\sQuote{filter.threshold} for method \dQuote{FFT} and
					\sQuote{degree}, \sQuote{window.sav} and \sQuote{smoothing} for method \dQuote{SavGol}.}

}
\value{Returns an object of type \sQuote{NDVI} containing raw data, corrected NDVI values and modelled NDVI values.}
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
\seealso{\code{\link{bise}}, \code{\link{runningAvg}}, \code{\linkS4class{NDVI}}, \code{\link{phenoPhase}}}
\examples{
data(avhrr)
data(modis)

# create NDVI object, correct and model values
ndvi.list1 <- modelNDVI(ndvi.values=cbind(avhrr.ndvi/10000, modis.ndvi/10000), 
			year.int=1995, correction="bise", method="LinIP", MARGIN=2, 
			doParallel=FALSE, slidingperiod=40)
ndvi.list2 <- modelNDVI(ndvi.values=cbind(avhrr.ndvi/10000, modis.ndvi/10000), 
			year.int=1995, correction="ravg", method="FFT", MARGIN=2, 
			doParallel=FALSE, filter.threshold=7)

# plot Values
for (ndvi.ob in ndvi.list1){ plot(ndvi.ob) } 
for (ndvi.ob in ndvi.list2){ plot(ndvi.ob) }
}

