\name{pairwiseRegions}
\alias{pairwiseRegions}

\title{Classify two Binding Profiles}

\description{Classify two Binding Profiles into unique and common binding regions and write results to bed files.  }

\usage{
pairwiseRegions(xSet, sgrset = c(1, 2), bound.cutoff, diff.cutoff, probes, probe.max.spacing, writeBedFile=TRUE)
}

\arguments{
  \item{xSet}{object of class ExpressionSet}
  \item{sgrset}{vector of lenght 2; specifying which data sets to compare from the ExpressionSet}
  \item{bound.cutoff}{numeric; threshold above probes are considered \dQuote{bound}}
  \item{diff.cutoff}{numeric; difference threshold to determine if object 1 and object 2 are unqiquely bound}
  \item{probes}{integer; minimum number of probes in a valid region}
  \item{probe.max.spacing}{integer; maximum number of base pairs in a gap before splitting a region into 2 regions}
  \item{writeBedFile}{logical; should bed file be written}
}


\details{
Probe signal values above the bound.cutoff in both data are classified as common 
bound. Probes which are above the bound.cutoff and in one data and higher than the
diff.cutoff to the other data are called unique. Then these probes are then filtered 
into regions using the probes and probe.max.spacing details. The score for the unique 
regions is calculated as mean (probes in region set 1 minus set 2), or vise versa. 
The score for the common region is the mean (probes in region (set 1 plus set 2)/2).
Optional bed file formated result files are written using the choosen options in the file names.  
}


\value{
data.frame with the following columns:
  \item{name}{name(s) of data set to which region belongs}
  \item{class.group}{class group; 1, 2 or 3 for common regions between both sets}
  \item{chr}{chromsome}
  \item{start}{start position of region}
  \item{end}{end position of region}
  \item{scrore}{score of region}
  \item{nProbes}{number of probes in region}
}

\author{Bettina Fischer, Robert Stojnic}


\seealso{
  \code{\link{compensationRegions}},
  \code{\link{increasedBindingRegions}},
  \code{\link{threewayRegions}}
}


\examples{
  dataPath <- system.file("data",package="SimBindProfiles")
  load(file.path(dataPath,"SGR.RData"))
  pairAB <- pairwiseRegions(SGR, sgrset=c(1,2), bound.cutoff=1.86, diff.cutoff=1.4, 
            probes=10, probe.max.spacing=200)
}
