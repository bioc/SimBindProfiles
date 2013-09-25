\name{readSgrFiles}
\alias{readSgrFiles}

\title{Read sgr data}

\description{Function to read sgr files into a list object, quantile normalise signal.  }

\usage{readSgrFiles(X, dataPath = getwd(), fileExt=".txt", normalise = TRUE)}

\arguments{
  \item{X}{object of class ExpressionSet}
  \item{dataPath}{path to the directory holding the data sgr data files}
  \item{normalise}{logical, should data be quantil normalised}
  \item{fileExt}{character specifying the filename extension}
}

\details{
The function reads the data files in sgr file tab delimited format: chr, position, signal. 
All data sets must be from the same array platform and have the same chromosome
names and positions. The data is read and then quantile normalised and stored as an 
ExpressionSet.  
}

\value{Returns normalized, transformed values as an object of class ExpressionSet}

\author{Bettina Fischer}

\seealso{
  \code{\link[Biobase]{ExpressionSet}}
}

\examples{
  sgrfiles <- c("SoxNDam_trunc","SoxN-DDam_trunc","DDam_trunc")
  dataPath <- system.file("extdata",package="SimBindProfiles")
  readTestSGR <- readSgrFiles(X=sgrfiles, dataPath)
}
