\name{eSetScatterPlot}
\alias{eSetScatterPlot}

\title{Scatterplot of ExpressionSet object}

\description{Create a smooth scatterplot with correlation  }

\usage{
    eSetScatterPlot(xSet)
}

\arguments{
    \item{xSet}{object of class ExpressionSet}
}


\details{Scatterplot and correlations of all data sets in ExpressionSet object}

\value{No useful return. The function is called for its side effect to produce the scatterplot.}

\author{Bettina Fischer}

\seealso{
    \code{\link[graphics]{plot}}, \code{\link[graphics]{smoothScatter}}
}

\examples{
    dataPath <- system.file("data", package = "SimBindProfiles")
    load(file.path(dataPath, "SGR.RData"))
    eSetScatterPlot(SGR)
}
