\name{SGR}
\alias{SGR}
\alias{precomputed}
\alias{bound.cutoff}
\alias{compR}
\alias{increasedR}
\docType{data}
\title{ A Sample Data object }
\description{ 
  The SGR datset consists of three data sets stored in an ExpressionSet object.
  These are SoxNDam, SoxN-DDam and DDam with probes located on chromosome X 
  between 1-6000000 bp. The precomputed dataset contains precomputed objects 
  for the bound.cutoff, compR and increasedR to save time for building the vignette.  
}

\usage{
  data(SGR)
  data(precomputed)
}

\format{
  An ExpressionSet object, with three data sets.
}

\examples{
  data(SGR)
  data(precomputed)
}

\keyword{datasets}
