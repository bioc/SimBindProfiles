\name{threewayRegions}
\alias{threewayRegions}

\title{Classify three Binding Profiles}

\description{
  Classify three Binding Profiles into unique and common binding regions and 
  write results to bed files.  
}

\usage{
threewayRegions(xSet, sgrset = c(1, 2, 3), bound.cutoff, diff.cutoff, 
  probes, probe.max.spacing, writeBedFile=TRUE)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{vector of lenght 3; specifying which data sets to compare from the ExpressionSet}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{diff.cutoff}{numeric; difference threshold to determine if objects are unqiquely bound}
  \item{probes}{integer; minimum number of probes in a valid region}
  \item{probe.max.spacing}{integer; maximum number of base pairs in a gap before splitting a region into 2 regions}
  \item{writeBedFile}{logical; should a bed file be written}
}

\details{
  Signal values are flagged as bound and unbound, and then the bound probes are classified into 7 groups:
  class.group=1: unique probes in object 1;
  class.group=2: unique probes in object 2;
  class.group=3: unique probes in object 3;
  class.group=4: common probes in object 1+2;
  class.group=5: common probes in object 2+3;
  class.group=6: common probes in object 1+3;
  class.group=7: common probes in object 1+2+3.
  The classified probes are then filtered into regions using the probes and probe.max.spacing details.
  For scores calculations please refer to the pairwiseRegions function. Optional bed file 
  formated result files are written using the choosen options in the file names.  
}


\value{
  \item{name}{name(s) of data set to which region belongs}
  \item{class.group}{class group; 1 to 7, see details above}
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
  \code{\link{increasedBindingRegions}}
}


\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  threewayRegions(SGR, sgrset=c(1,2,3), bound.cutoff=1.86, diff.cutoff=1.4, 
                  probes=10, probe.max.spacing=200)
}
