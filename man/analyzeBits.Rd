\name{analyzeBits}
\alias{analyzeBits}

\title{Bit Analysis}
\description{Analyses bits of a 16-bit integer}
\usage{analyzeBits(value, mode=1, bitpos=0)}
\arguments{
  \item{value}{A 16-bit integer value}
  \item{mode}{If mode is 0, the algorithm returns the bit on position \sQuote{bitpos}. If mode is 1 (default), the algorithm returns the most significant bit, 
if mode is 2 it returns the second significant bit and if mode is 3 the return value is the integer value of the last four bits}
  \item{bitpos}{An integer value between 0 and 15 determining the position of the bit to read when mode is 0.}
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

