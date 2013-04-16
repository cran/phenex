\name{date2doy}
\alias{date2doy}

\title{Date to Julian Day Converter}
\description{Converts a date into a Julian day}
\usage{date2doy(date)}
\arguments{
	\item{date}{Date YYMMDD as integer}
}
\value{The Julian Day (day of year) of the date}
\author{Daniel Doktor, Maximilian Lange}

\examples{
	date <- 891208
	doy <- 0
	doy <- date2doy(date)
	doy
}