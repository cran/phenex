\name{yearlength}
\alias{yearlength}

\title{Number of Days}
\description{Returns the number of days of the given year.}
\usage{yearlength(year)}
\arguments{
	\item{year}{A vector of years as integer values.}
}
\value{The number of days of the year in the date}
\author{Daniel Doktor, Maximilian Lange}

\examples{
	year <- c(1995, 2000, 2005, 2010)
	days <- yearlength(year)
	days
}