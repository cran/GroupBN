\name{discretize.dens}
\alias{discretize.dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{discretize.dens
%%  ~~function to do ... ~~
}
\description{
density approximative discretization. Significant peaks in the density are determined and used as starting points for k-means based discretization. If only one peak is present, distribution quartiles are used for binning.
}
\usage{
discretize.dens(data, graph=F, title="Density-approxmative Discretization",
rename.level=F, return.all=T, cluster=F, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a vector containing the data that may be discretized}
  \item{graph}{a boolean value, if TRUE, the density and the determined binning are plotted}
  \item{title}{a title for the plot}
  \item{rename.level}{a boolean value, if TRUE, factor levels are replaced by integers 1:n}
  \item{return.all}{a boolean value, if FALSE, only the discretized data are returned.}
  \item{cluster}{a boolean value, if data is a cluster variable and may already be discrete or not}
  \item{seed}{a random seed number}
}
\value{
\item{discretized}{the discretized data}
\item{levels}{the factor levels}
\item{optima}{the x and y coordinates of the determined peaks}
}
\author{Ann-Kristin Becker
%%  ~~who you are~~
}
\examples{
testdata = c(rnorm(100,-3,1), rnorm(100,3,1))
d<-discretize.dens(testdata, graph=TRUE)
summary(d$discretized)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
