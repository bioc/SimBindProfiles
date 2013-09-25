\name{probeLengthPlot}
\alias{probeLengthPlot}

\title{Plot the probe length frequency}

\description{Plot the probe length frequency.  }


\usage{
  probeLengthPlot(xSet, sgrset = 1, chr = NULL, bound.cutoff, probe.max.spacing = 200, xlim.max = 25)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{integer; specifying which data set to use from the ExpressionSet}
  \item{chr}{probes from which chromosome to be used, default NULL means all probes}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{probe.max.spacing}{integer; maximum amount of base pairs at which bound probes are condensed into one region.}
  \item{xlim.max}{integer; maximum number of probes plotted along the x-axsis}
}


\details{Plot the probe length frequency per regions based on a bound.cutoff value.  }

\value{Plot of frequency of the number of probes within bound regions}

\author{Bettina Fischer}

\seealso{
  \code{\link{plot}}
}

\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  probeLengthPlot(SGR, sgrset=1, chr=NULL, bound.cutoff=1.2, probe.max.spacing=200, xlim.max=50)
}
