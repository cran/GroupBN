\name{plot.groupbn}
\alias{plot.groupbn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.groupbn
%%  ~~function to do ... ~~
}
\description{
generic plot function for class groupbn
}
\usage{
\method{plot}{groupbn}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class groupbn}
  \item{...}{further arguments}
}
\details{
Plot the group bayesian network structure
}
\author{Ann-Kristin Becker
%%  ~~who you are~~
}
\value{
No return value, called for plotting
}

\seealso{
\code{\link{groupbn}}
}
\examples{
data("wine.groupbn.refined")
plot(wine.groupbn.refined)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

