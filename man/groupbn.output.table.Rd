\name{groupbn.output.table}
\alias{groupbn.output.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{groupbn.output.table
%%  ~~function to do ... ~~
}
\description{Create an output table with clusters and included variables with similarity scores
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
groupbn.output.table(res, with.scores=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{gn object}
  \item{with.scores}{if TRUE, similarity scores of every cluster member to the cluster center are added to the table}
%%     ~~Describe \code{x} here~~
}
\value{a table with one column per group,
similarity scores to cluster centers are calculated for each variable
}
\author{Ann-Kristin Becker
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{groupbn}}
\code{\link{groupbn_refinement}}
}
\examples{
data("wine.groupbn.refined")
df<-groupbn.output.table(wine.groupbn.refined)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

