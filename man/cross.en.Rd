\name{cross.en}
\alias{cross.en}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{cross.en
%%  ~~function to do ... ~~
}
\description{
Calculates the weighted cross entropy / log-loss for a vector of observations and predicted probabilities (weighted by class proportions)
}
\usage{
cross.en(pred, obs, sdpred=NULL, weighted=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{a numeric vector, the predicted probabilities of the reference class}
  \item{obs}{the vector of observations, a categorical variable with 2-4 levels}
  \item{sdpred}{either NULL or a vector containing the standard deviations of every estimate}
  \item{weighted}{a boolean, if FALSE, the unweighted logloss is calculated. By default, the weighted cross entropy is calculated.}
}
\details{
if sdpred contains the standard deviations for each estimated probability, then a lower bound of the log loss is returned.
}
\value{a numeric value: cross entropy / log loss for comparison of classifiers. The smaller, the better.
}
\author{ Ann-Kristin Becker
%%  ~~who you are~~
}
\examples{
#observations
obs<-as.factor(c("A","A","B"))
#correct prediction
pred1<-c(1,1,0)
#wrong prediction
pred2<-c(0,0,1)

cross.en(pred=pred1, obs=obs) #small
cross.en(pred=pred2, obs=obs) #large

#prediction of only majority class
pred3<-c(1,1,1)
#prediction of only minority class
pred4<-c(0,0,0)

cross.en(pred=pred3, obs=obs, weighted=TRUE)
cross.en(pred=pred4, obs=obs, weighted=TRUE)
#both equal (as weighted)

cross.en(pred=pred3, obs=obs, weighted=FALSE)
cross.en(pred=pred4, obs=obs, weighted=FALSE)
#unweighted, majority class is favored
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

