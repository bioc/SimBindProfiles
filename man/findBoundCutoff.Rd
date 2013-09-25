\name{findBoundCutoff}
\alias{findBoundCutoff}

\title{Find the bound.cutoff  }

\description{
  Define the bound.cutoff using either the normalNull or the twoGaussiansNull method.  
}

\usage{
  findBoundCutoff(xSet, method = c("normalNull", "twoGaussiansNull"), mean.method = "mode", 
    pvalue = FALSE, fdr = FALSE, pvalPlot = FALSE)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{method}{a character string equal to \dQuote{normalNull} or \dQuote{twoGaussiansNull}}
  \item{mean.method}{a characyer string equal to \dQuote{mode} or \dQuote{zero}. This is used in combination with the \dQuote{normalNull} method.}
  \item{pvalue}{decimal specifying the p-value cutoff (either pvalue or fdr can be used not both!)}
  \item{fdr}{decimal specifying the fdr cutoff (either pvalue or fdr can be used not both!)}
  \item{pvalPlot}{logical, if TRUE the pvalue histogram is written}
}

\details{
  We implemented two methods to set the bound.cutoff, probes above this threshold are considered \dQuote{bound}. 
  The twoGaussiansNull method established in the Ringo package (Toedling et al., 2007), by which the data 
  is assumed to follow a mixture of two Gaussian distributions. The one Gaussian with the lower mean value 
  is assumed to be the null distribution and probe levels are assigned p-values based on this null 
  distribution. Alternatively the user can select the normalNull method instead which assumes the null 
  distribution is normal and symmetrical around the mode (or zero). For both methods the user can decide if the 
  resulting p-values are to be adjusted for multiple testing (fdr) or selected by a p-value threshold.
  The function also provides QC plots for the twoGaussiansNull and an optional p-value histogram.  
}

\value{Returns a numeric}

\references{
  Toedling J., Skylar O., Krueger T, Fischer J.J., Sperling S., Huber W. 2007 Ringo - an
  R/Bioconductor package for analyzing ChIP-chip readouts. BMC Bioinformatics, 8:221
}

\note{
  Please note that the use of the package \dQuote{mclust} is only free for strict academic 
  use (see the license of \dQuote{mclust} here: http://www.stat.washington.edu/mclust/license.txt). 
  The alternative function normalNull does not have this restriction.  
}

\author{Bettina Fischer}

\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  bound.cutoff <- findBoundCutoff(SGR, method="twoGaussiansNull", fdr=0.25)
}
