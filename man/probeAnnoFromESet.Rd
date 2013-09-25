\name{probeAnnoFromESet}
\alias{probeAnnoFromESet}

\title{Build a probeAnno from ExpressionSet}

\description{Build a probeAnno from ExpressionSet.  }

\usage{
  probeAnnoFromESet(eSet, probeLength)
}

\arguments{
  \item{eSet}{object of class ExpressionSet}
  \item{probeLength}{integer; specifying probe lenght on array}
}


\details{
  Creates a probeAnno object from ExpressionSet object. The function uses the
  PROBE_ID, CHROMOSOME and POSITION information store in the ExpressionSet object.  
}


\value{
An object of class probeAnno holding the mapping between probes and 
genomic positions.
}

\author{Bettina Fischer}


\seealso{
  \code{\link[Ringo]{probeAnno-class}}
}


\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  probeAnno <- probeAnnoFromESet(SGR, probeLength=50)
}
