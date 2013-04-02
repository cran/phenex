\name{analyzeBits}
\alias{analyzeBits}

\title{Bit Analysis}
\description{Analyses bits of a 16-bit integer}
\usage{analyzeBits(value, mode)}
\arguments{
  \item{value}{A 16-bit integer value}
  \item{mode}{If mode is 1 (default), the algorithm returns the most significant bit, 
if mode is 2 it returns the second significant bit and if mode is 3 the return value is the integer value of the last four bits}
}
\details{This routine analyses 16-bit integer values to get the indicators of MEDOKADS data}
\value{An integer value respective to parameter \sQuote{mode}.}
\author{Daniel Doktor, Maximilian Lange}
\references{
Koslowsky, D., Billing, H. and Friedrich, K. (2005): MEDOKADS: A long-term data set for detection and monitoring for desertification risks in the mediterranean. In \emph{RGLDD Conference}.
}
\examples{
value <- -32768
res <- analyzeBits(value, mode=3)
res
}

