\name{leapYears}
\alias{leapYears}

\title{Leap Year Check}
\description{Checks whether the given years are leap years.}
\usage{leapYears(year)}
\arguments{
	\item{year}{A vector of years as integer values.}
}
\value{A vector of boolean values determining whether the given years are leap years.}
\author{Daniel Doktor, Maximilian Lange}

\examples{
	year <- c(1995, 2000, 2005, 2010)
	leaps <- leapYears(year)
	leaps
}