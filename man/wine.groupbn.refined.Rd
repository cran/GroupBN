\name{wine.groupbn.refined}
\alias{wine.groupbn.refined}
\docType{data}
\title{wine.groupbn.refined
%%   ~~ data name/kind ... ~~
}
\description{
A refined group Bayesian network with 8 groups learned from dataset 'wine'.
}
\usage{data("wine.groupbn.refined")}
\format{group Bayesian network (class 'groupbn')

name of target variable: Soil
number of groups:  8
achieved scoring: F1:  0.92 ; Precision:  1 ; Recall:  0.86 ; AUC-PR:  1 ; AUC-ROC:  1 ; cross-entr.:  1.43; BIC (netw.): -77.21

 name          description
 "$bn"         "Bayesian network structure"
 "$fit"        "fitted Bayesian network (multinomial)"
 "$arc.confid" "arc confidence"
 "$X.quali"    "qualitative variables in a data.frame"
 "$X.quanti"   "quantitative variables in a data.frame"
 "$grouping"   "group memberships"
 "$k"          "number of groups of initial grouping"
 "$group.data" "group representatives used for network inference"
 "$target"     "name of target variable"
 "$separate"   "name of any other separated variables"
 "$pca.param"  "pca parameters of each group"
 "$disc.param" "discretization intervals of each group"
 "$score"      "cross entropy and additional scoring information"
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(wine.groupbn.refined)
}
\keyword{datasets}
