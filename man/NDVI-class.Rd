\name{NDVI-class}
\Rdversion{1.1}
\docType{class}
\alias{NDVI-class}
\alias{bise,NDVI-method}
\alias{bise1,NDVI-method}
\alias{bise2,NDVI-method}
\alias{checkLength,NDVI-method}
\alias{correctedValues<-,NDVI-method}
\alias{correctedValues,NDVI-method}
\alias{isLeapYear,NDVI-method}
\alias{modelledValues<-,NDVI-method}
\alias{modelledValues,NDVI-method}
\alias{modelValues,NDVI-method}
\alias{phenoPhase,NDVI-method}
\alias{plot,NDVI-method}
\alias{runningAvg<-,NDVI-method}
\alias{runningAvg,NDVI-method}
\alias{values<-,NDVI-method}
\alias{values,NDVI-method}
\alias{year<-,NDVI-method}
\alias{year,NDVI-method}
\alias{rsquare,NDVI-method}
\alias{integrateTimeserie,NDVI-method}

\title{Class \code{"NDVI"}}
\description{
	Class \sQuote{NDVI} provides functions to smooth NDVI (Normalized Difference Vegetation Index) 
	timeseries obtained from satellite observations.
}
\section{Objects from the Class}{
	NDVI-Objects can be created by calls of the form \code{modelNDVI(ndvi.values, ...)}.
	These objects contain the timeserie as given by \sQuote{ndvi.values}, smoothed values due 
	to the performed correction (see also \code{\link{modelNDVI}}, parameter 'correction') 
	and modelled values due to the chosen method (see also \code{\link{modelNDVI}}, parameter \sQuote{method}).
	Furthermore, the year of the timeserie can be stored as integer value.
}

\section{Slots}{
  \describe{
    \item{\code{year}:}{The year of the timeserie as integer value or \sQuote{NA}.}
    \item{\code{values}:}{A vector of NDVI-values (between minus one and one as numeric) of length 365 or 
	366 (respective to length of the year).}
    \item{\code{correctedValues}:}{A vector of corrected NDVI-values. Is set after use of correction methods \sQuote{bise} or \sQuote{runningAvg}.}
    \item{\code{modelledValues}:}{A vector of corrected NDVI-values. Is set after use of method \sQuote{modelValues}.}
  }
}
\section{Methods}{
  \describe{
    \item{bise1}{\code{signature(x = "NDVI")}: Best index slope extraction after \emph{Viovy et. al} (1992). 
			See also \code{\link{bise1}} and \code{\link{bise2}}.

		This routine tries to restore the temporal NDVI profile, 
		i.e. separate true observations from noise. The sliding period default of BISE-algorithm 
		is 40 days, optimised for the area of Central Europe. This might has to be modified 
		when vegetation dynamics of other climates/habitats are analysed.

		This method is used by function \code{\link{modelNDVI}} if parameter correction is set to \dQuote{bise1}.}
    \item{bise2}{\code{signature(x = "NDVI")}: Best index slope extraction after \emph{Viovy et. al} (1992). 
			Second interpretation of algorithm. See also \code{\link{bise1}} and \code{\link{bise2}}.

		This routine tries to restore the temporal NDVI profile, 
		i.e. separate true observations from noise. The sliding period default of BISE-algorithm 
		is 40 days, optimised for the area of Central Europe. This might has to be modified 
		when vegetation dynamics of other climates/habitats are analysed.
		An maximum increase threshold is defined via parameter \sQuote{growthFactorThreshold} with 
		default to 0.1 (10 percent increase per day allowed).
		The timeserie is cycled via parameter \sQuote{cycleValues}, which is a 
		boolean value determining whether the end of the ndvi timeserie 
		is combined with its beginning or not (default value is true). If false, gaps or low ndvi 
		values at the beginning of timeserie influence the result.

		This method is used by function \code{\link{modelNDVI}} if parameter correction is set to \dQuote{bise2}.}
    \item{bise}{\code{signature(x = "NDVI")}: Calls method \code{\link{bise2}}.}
    \item{checkLength}{\code{signature(x = "NDVI")}: Checks the length of the timeserie respective 
			to the length of the year (365 days or 366 for leap years).}
    \item{correctedValues<-}{\code{signature(x = "NDVI")}: Replacement method for slot \code{correctedValues}.}
    \item{correctedValues}{\code{signature(x = "NDVI")}: Acessor method for slot \code{correctedValues}. }
    \item{isLeapYear}{\code{signature(x = "NDVI")}: Checks whether the year of the NDVI-object is a leap year.}
    \item{modelledValues<-}{\code{signature(x = "NDVI")}: Replacement method for slot \code{modelledValues}. }
    \item{modelledValues}{\code{signature(x = "NDVI")}: Acessor method for slot \code{modelledValues}. }
    \item{modelValues}{\code{signature(x = "NDVI")}: Fits a suite of functions/models to corrected NDVI-timeserie (if the 
				corrected timeserie is not available, the raw one out of \sQuote{values} will be used). 
				Parameter \sQuote{method} determines the used model:

				\dQuote{\bold{LinIP}}: A linear interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1). 
				Applied in \emph{Badeck et. al} (2004) and \emph{Doktor et. al} (2009).

				\dQuote{\bold{Spline}}: A spline interpolation is performed. For interpolation, the end of 
				timeserie is connected to the beginning (e.g. after day 365 follows day 1).

				\dQuote{\bold{DSig}}: Fits a double sigmoidal function to NDVI values (according to \emph{Zhang et. al} (2003)).

				\dQuote{\bold{DSigC}}: Fits another double sigmoidal function (own C implementation) to NDVI values.

				\dQuote{\bold{DLogistic}}: Fits a double logistic function after \emph{Fischer, Alberte} (1994) to NDVI values.

				\dQuote{\bold{Gauss}}: Fits a symmetric or asymmetric (determined by boolean parameter \sQuote{asym}) gaussian 
					function to NDVI values (own C implementation after \emph{Press, W.H.} (1992)).

				\dQuote{\bold{Growth}}: Fits a growth model after \emph{Richter et. al} (1991) to NDVI values.

				\dQuote{\bold{FFT}}: Smoothes the corrected or raw NDVI values with fast fourier transfusion (implemented in R). 
					The smoothing intensity can be controlled with parameter \sQuote{filter.threshold} with default to 3.

				\dQuote{\bold{SavGol}}: Smoothes the corrected or raw NDVI values with a Savitzky-Golay 
					filter (own C implementation after \emph{Press, W.H.} (1992)). 
					The smoothing algorithm can be modified with parameters 
					\sQuote{window.sav} (window size of filter, default to 7),
					\sQuote{degree} (degree of fitting polynomial, default to 2) 
					and \sQuote{smoothing} (repetition quantity, default to 10).

				This method is used by function \code{\link{modelNDVI}} with modelling method respective to parameter \sQuote{method}.}
    \item{phenoPhase}{\code{signature(x = "NDVI")}: Extracts the start of phenological phases out of the modelled NDVI values. 
					Parameter \sQuote{phase} determines which phase is extracted:

					\dQuote{\bold{max}}: Day of the year with highest NDVI value is returned.
	
					\dQuote{\bold{min}}: Day of the year with lowest NDVI value is returned. 

					\dQuote{\bold{greenup}}: Day of the year at which greenup takes place is returned. 
						This day refers to the point where the function of modelled NDVI 
						values exceeds a certain threshold.
			
					\dQuote{\bold{senescence}}: Day of the year at which senescence takes place is returned. 
						This day refers to the point where the function of modelled NDVI 
						values exceeds a certain threshold.

					Parameter \sQuote{method} determines whether a global or local threshold is used for
					greenup and senescence extraction:
		
					With \dQuote{\bold{global}} threshold, the day of the year is returned, 
					where NDVI values are first equal or higher as the value of \sQuote{threshold}. 
					If the threshold is higher than the values of the timeseries, \sQuote{-1} will be returned.
					
					With \dQuote{\bold{local}} threshold, the day of the year is returned, for which NDVI values first reach the value of
					\sQuote{threshold} (interpreted as percentage) between lowest and highest NDVI value of timeserie. 
					The lowest NDVI value is chosen depending on phase selected. For \dQuote{greenup}, 
					the lowest value before day of maximum NDVI value is used. For \dQuote{senescence}, 
					the lowest value after day of maximum NDVI value is used.

					The \sQuote{threshold} for local or global greenup/senescence-extraction method
					should have numerical value between \sQuote{0} and \sQuote{1}. Global thresholds 
					refer to a fix NDVI value independent of actual NDVI magnitudes. 
					Local thresholds are based on (max - min) * e.g. 0.55. Hence, the respective 
					NDVI value will vary dependent on e.g. land-cover.
					}
    \item{plot}{\code{signature(x = "NDVI")}: Plots raw data as black dots (slot 'values'), 
			corrected values as red dots (slot \sQuote{correctedValues}) and 
			modelled values as blue line (slot \sQuote{modelledValues}).}
    \item{runningAvg}{\code{signature(x = "NDVI")}: This routine performes an averaging with a 
			running window on NDVI values. Default window size is 7 and can be modified by parameter \sQuote{window}.

			This method is used by function \code{\link{modelNDVI}} if parameter correction is set to \dQuote{ravg}.
			See also \code{\link{runningAvg}}.}
    \item{values<-}{\code{signature(x = "NDVI")}: Replacement method for slot \code{values}. }
    \item{values}{\code{signature(x = "NDVI")}: Acessor method for slot \code{values}. }
    \item{year<-}{\code{signature(x = "NDVI")}: Replacement method for slot \code{year}. }
    \item{year}{\code{signature(x = "NDVI")}: Acessor method for slot \code{year}. }
    \item{rsquare}{\code{signature(x = "NDVI")}: Calculates the squared Pearson correlation coefficient out of 
				corrected (if the corrected timeserie is not available, the raw one out of 
				\sQuote{values} will be used) and modelled timeserie. }
    \item{integrateTimeserie}{\code{signature(x = "NDVI")}: Calculates the integral of the timeserie between the bounds \sQuote{start} and \sQuote{end}.}
    }
}

\seealso{\code{\link{bise}}, \code{\link{runningAvg}}, \code{\link{modelNDVI}}, \code{\link{phenoPhase}}}
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

\examples{
	 ##first example
	# load data
	data(avhrr)
	data(modis)

	# create NDVI object
	ndvi1 <- new("NDVI", values=avhrr.ndvi/10000, year=as.integer(1995))
	ndvi2 <- new("NDVI", values=modis.ndvi/10000, year=as.integer(1995))

	# correct values (bise)
	ndvi1.bise <- bise(ndvi1)
	ndvi2.bise <- bise(ndvi2)

	# or running Average
	ndvi1.ravg <- runningAvg(ndvi1)
	ndvi2.ravg <- runningAvg(ndvi2)

	# model Values
	res1 <- modelValues(ndvi1.bise, method="LinIP")
	res2 <- modelValues(ndvi1.ravg, method="FFT", filter.threshold=7)
	res3 <- modelValues(ndvi2.bise, method="LinIP")
	res4 <- modelValues(ndvi2.bise, method="FFT", filter.threshold=7)

	# plot Values
	plot(res1)
	plot(res2)
	plot(res3)
	plot(res4)

	# extract greenup date
	greenup <- phenoPhase(res1, phase="greenup", method="local", threshold=0.55)
	# extract date with highest ndvi
	max.ndvi <- phenoPhase(res1, phase="max")

	## or simpler
	data(avhrr)
	data(modis)
	
	# create NDVI objects, correct and model values
	ndvi.list1 <- modelNDVI(ndvi.values=cbind(avhrr.ndvi/10000, modis.ndvi/10000),
			year.int=1995, correction="bise", method="LinIP", MARGIN=2,
			doParallel=FALSE, slidingperiod=40)
	ndvi.list2 <- modelNDVI(ndvi.values=cbind(avhrr.ndvi/10000, modis.ndvi/10000),
			year.int=1995, correction="ravg", method="FFT", MARGIN=2,
			doParallel=FALSE, filter.threshold=7)

	# plot Values
	for (ndvi.ob in ndvi.list1){ plot(ndvi.ob) }
	for (ndvi.ob in ndvi.list2){ plot(ndvi.ob) }
	
	# extract greenup date
	greenup <- phenoPhase(ndvi.list1[[1]], phase="greenup", 
			method="local", threshold=0.55)
	# extract senescence date
	senescence <- phenoPhase(ndvi.list1[[1]], phase="senescence", 
			method="local", threshold=0.55)
	# extract date with highest NDVI
	max.ndvi <- phenoPhase(ndvi.list1[[1]], phase="max")

	# calculate green season integrated vegetation index
	gsivi <- integrateTimeserie(ndvi.list1[[1]], start=greenup, end=senescence)
}
