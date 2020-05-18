\name{groupbn_refinement}
\alias{groupbn_refinement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{groupbn_refinement
%%  ~~function to do ... ~~
}
\description{Adaptive Refinement of a group Bayesian Network
using hierarchical Clustering
}
\usage{
groupbn_refinement(res, hierarchy, refinement.part="mb", restart=0, perturb=1,
max.step=10, max.min=Inf, R=100, return.all=FALSE, debug=FALSE, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{an object of class groupbn}
  \item{hierarchy}{a cluster object from ClustOfVar}
  \item{refinement.part}{"mb", "mb2" or "all", selects if the refinement steps should be done only within the markov blanket of the target variable, within the second-order markov blanket or within all clusters. Default: "mb"}
  \item{restart}{a positive integer number, the number of restarts}
  \item{perturb}{a positive integer number, the number of perturbations (splits) in each restart}
  \item{max.step}{a positive integer number, the maximal number of refinement steps, default is 10}
   \item{max.min}{a positive integer number, the maximal run time in minutes, default is unlimited}
  \item{R}{number of bootstrap replicates for model averaging, default is 100}
  \item{return.all}{a boolean, if TRUE, the output is a whole list of group models, if FALSE, the output is only the best-scoring model.}
  \item{debug}{a boolean,
  if TRUE, debugging messages are printed}
  \item{seed}{a random seed number}
}
\details{Based on a variable grouping, data are aggregated and a Bayesian network is learned. The target variable is kept separated during this procedure, so that the resulting network model can be used for risk prediction and classification. Starting from a coarse group network, groups are iteratively refined to smaller groups. The heuristic refinement happens downwards along the dendrogram, and stops, if it no longer improves the predictive performance of the model. The refinement part is implemented using a hill-climbing procedure.
}
\value{
returns an object of class groupbn
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Ann-Kristin Becker
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{groupbn}}
\code{\link{groupbn.output.table}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#load example data
data(wine)
wine.test<-wine[wine$Soil\%in\%c("Reference", "Env1"),1:29]
wine.test$Soil<-factor(wine.test$Soil)
levels(wine.test$Soil)<-c("0", "1")

#cluster data
hierarchy<-hclustvar(X.quanti=wine.test[,3:29], X.quali=wine.test[,1:2])

#Learn group network among 5 clusters with "Soil" as target variable
wine.groupbn<-groupbn(hierarchy, k=5, target="Soil", separate=NULL,
X.quanti=wine.test[,3:29], X.quali=wine.test[,1:2], seed=321)

#Do one refinement step
#Set max.step higher to optimize completely
wine.groupbn.refined<-groupbn_refinement(wine.groupbn, hierarchy,
refinement.part="mb", max.step = 1, seed=321)

#Plot refined network
plot(wine.groupbn.refined)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

