\name{groupbn_refine_manually}
\alias{groupbn_refine_manually}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{groupbn_refine_manually
%%  ~~function to do ... ~~
}
\description{Based on a GroupBN, a cluster can be selected manually, that is split and the refined model is learned.
}
\usage{
groupbn_refine_manually(res, hierarchy, refine, arc.thresh=NULL,
R=100, debug=FALSE, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{an object of class groupbn}
  \item{hierarchy}{a cluster object from ClustOfVar}
  \item{refine}{name of group to be refined}
  \item{arc.thresh}{threshold for bootstrap arcs}
  \item{R}{number of bootstrap replicates for model averaging, default is 100}
  \item{debug}{a boolean,
  if TRUE, debugging messages are printed}
  \item{seed}{a random seed number}
}
\value{
returns an object of class groupbn
}
\author{
Ann-Kristin Becker
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

#Refine cluster 2
wine.groupbn.refined<-groupbn_refine_manually(wine.groupbn, hierarchy,
refine = "cl2", seed=321)

#Plot refined network
plot(wine.groupbn.refined)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
