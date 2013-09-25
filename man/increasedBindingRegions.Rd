\name{increasedBindingRegions}
\alias{increasedBindingRegions}

\title{Classify two Binding Profiles - Increased Binding}

\description{
  This function is used to classify two Binding Profiles into regions which 
  are more bound in one data set than the other.  
}


\usage{
increasedBindingRegions(xSet, sgrset = c(1, 2), bound.cutoff, diff.cutoff, probes, 
  probe.max.spacing, writeBedFile = TRUE)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{vector of lenght 2; specifying which data sets to compare from the ExpressionSet}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{diff.cutoff}{numeric; difference threshold to determine if object 1 and object 2 are uniquely bound}
  \item{probes}{integer; minimum number of probes in a valid region}
  \item{probe.max.spacing}{integer; maximum number of base pairs in a gap before splitting a region into 2 regions}
  \item{writeBedFile}{logical; should a bed file be written}
}


\details{
Probe signal values above the bound.cutoff threshold for both data sets are 
compared where set 1 is above the diff.cutoff of set 2. These probes are 
then filtered into regions using the probes and probe.max.spacing details. The
score is calculated as mean (probes in region set 1 minus set 2). Optional bed 
file formated result files are written using the choosen options in the file names.  
}


\value{
data.frame with the following columns:
  \item{name}{name(s) of data set to which region belongs}
  \item{class.group}{class group, in this case only 1}
  \item{chr}{chromsome}
  \item{start}{start position of region}
  \item{end}{end position of region}
  \item{scrore}{score of region}
  \item{nProbes}{number of probes in region}
}

\author{Bettina Fischer, Robert Stojnic}

\seealso{
  \code{\link{pairwiseRegions}}, 
  \code{\link{compensationRegions}}, 
  \code{\link{threewayRegions}}, 
  \code{\link{plotBoundProbes}}
}

\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  overcompBC <- increasedBindingRegions(SGR, sgrset=c(2,3), bound.cutoff=1.86, diff.cutoff=1.4, 
                probes=10, probe.max.spacing=200)
}
