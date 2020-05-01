\name{groupbn}
\alias{groupbn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{groupbn
%%  ~~function to do ... ~~
}
\description{creates groupbn object
(determines an initial clustering based on a hierarchy with target variable and 'separated' variables separated, learns a Bayesian network from grouped data and saves discretization and pca parameters)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
groupbn(hierarchy, k, target, separate, X.quanti, X.quali,
debug=FALSE, R=100, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hierarchy}{a cluster object from ClustOfVar.}
  \item{k}{a positive integer number, the number of initial clusters.}
  \item{target}{a string, the name of the target variable.}
  \item{separate}{a vector of strings, names of variables that should be separated from the groups, such as age, sex,...}
  \item{X.quanti}{a numeric matrix of data,
  or an object that can be coerced to such a matrix
  (such as a numeric vector or a data frame with all numeric columns).}
  \item{X.quali}{a categorical matrix of data,
  or an object that can be coerced to such a matrix
  (such as a character vector, a factor or a data frame with all factor columns).}
  \item{debug}{a boolean,
  if TRUE, debugging messages are printed}
  \item{R}{number of bootstrap replicates for model averaging, default is 100}
  \item{seed}{a random seed number}
}
\value{an object of class groupbn
  \item{bn}{a Bayesian Network structure of bn class from bnlearn.}
  \item{fit}{a Bayesian Network with fitted parameters of bn.fit class from bnlearn.}
  \item{X.quanti}{a data.frame containing only the quantitative variables.}
  \item{X.quali}{a data.frame containing only the qualitative variables.}
  \item{grouping}{a vector of positive integers, giving the cluster assignment.}
  \item{k}{the number of clusters.}
  \item{group.data}{a data.frame containing the cluster representants.}
  \item{target}{a string, the name of the target variable.}
  \item{separate}{a vector of strings, names of variables that should be separated from the groups.}
  \item{pca.param}{the PCAmix used to determine the cluster representants.}
  \item{disc.param}{the cutpoints used to discretize the cluster representants.}
  \item{score}{Different prediction scores for the target variable using the fitted network.}
}
\author{
Ann-Kristin Becker
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{groupbn_refinement}}
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

#Plot network
plot(wine.groupbn)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
