\name{groupbn.vis.html.plot}
\alias{groupbn.vis.html.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{groupbn.vis.html.plot
%%  ~~function to do ... ~~
}
\description{
Create an interactive html network object with visNet (displaying similarity scores and number of variables in a score)
}
\usage{
groupbn.vis.html.plot(res, df=NULL, save.file=TRUE, save.name=NULL,
hierarchical=FALSE, nodecolor.all="#E0F3F8",
nodecolor.special="cornflowerblue", main=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{a groupbn object}
  \item{df}{output from output.table if already calculated, otherwise the same table is calculated internally}
  \item{save.file}{boolean; if TRUE a html file is produced}
  \item{save.name}{name for saving html object, date is additionally used}
  \item{hierarchical}{boolean; if TRUE the network is plotted with a hierarchical layout}
  \item{nodecolor.all}{a color for "normal" nodes}
  \item{nodecolor.special}{a color for the target variable and all separated nodes, if any.}
  \item{main}{optionally a title for the plot}
}
\value{
an html widget of class visNetwork
}
\details{Plots an interactive network plot using visNetwork package
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
groupbn.vis.html.plot(wine.groupbn.refined, hierarchical=TRUE, save.file=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
