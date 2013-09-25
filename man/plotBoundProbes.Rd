\name{plotBoundProbes}
\alias{plotBoundProbes}

\title{Plot coloured bound probes}

\description{Scatterplot of bound probes in colour based on the method.  }

\usage{
plotBoundProbes(xSet, sgrset, method=c("pairwise" , "compensation" , "increasedBinding"), 
                 bound.cutoff, diff.cutoff, cols=NULL, pcex=2)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{integer; specifying which data set to use from the ExpressionSet}
  \item{method}{a character string equal to \dQuote{pairwise} , \dQuote{increasedBinding} , \dQuote{compensation}}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{diff.cutoff}{numeric; difference threshold to determine if object 1 and object 2 are uniquely bound}
  \item{cols}{vector of colours to highlight probes, otherwise colours are set by default}
  \item{pcex}{a numerical vector giving the amount by which probe symbols should be scaled relative to the default}
}

\details{Scatterplot of the bound probes in colour based on the selected method}

\value{Coloured scatter plot}

\author{Bettina Fischer}

\seealso{
  \code{\link[graphics]{plot}},
  \code{\link{pairwiseRegions}}, 
  \code{\link{compensationRegions}}, 
  \code{\link{increasedBindingRegions}}
}

\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  plotBoundProbes(SGR, sgrset=c(1,2), method="pairwise", bound.cutoff=1.86, diff.cutoff=1.4)
}
