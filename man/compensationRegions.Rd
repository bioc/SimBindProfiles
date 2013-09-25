\name{compensationRegions}
\alias{compensationRegions}

\title{Classify two Binding Profiles - Compensation}

\description{
  This function is used to classify three Binding Profiles into regions 
  which are bound in data set one and two but not three.  
}

\usage{
compensationRegions(xSet, sgrset = c(1, 2, 3), bound.cutoff, diff.cutoff, 
  probes, probe.max.spacing, writeBedFile = TRUE)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{vector of lenght 3; specifying which data sets to compare from the ExpressionSet}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{diff.cutoff}{numeric; difference threshold to determine if objects uniquely bound}
  \item{probes}{integer; minimum number of probes in a valid region}
  \item{probe.max.spacing}{integer; maximum number of base pairs in a gap before splitting a region into 2 regions}
  \item{writeBedFile}{logical; should bed file be written}
}


\details{
  Select probes with a signal above the bound.cutoff in data set 1 and 2 and below in set 3,
  and for which the average signal of set 1 and set2 is above the diff.cutoff compared
  to set 3. These probes are then filtered into regions using the probes and 
  probe.max.spacing details. The score is calculated as mean (probes in region
  average(set 1 and set 2) minus set 3). Optional bed file formated result files 
  are written using the choosen options in the file names.  
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
  \code{\link{threewayRegions}}, 
  \code{\link{increasedBindingRegions}}
}

\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  transcompABC <- compensationRegions(SGR, sgrset=c(1,2,3), bound.cutoff=1.86, 
                  diff.cutoff=1.4, probes=10, probe.max.spacing=200)
}
