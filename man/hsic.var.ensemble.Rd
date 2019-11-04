\name{hsic.var.ensemble}
\alias{hsic.var.ensemble}
\title{
Error Variance Estimation in Genomic Prediction
}
\description{
Estimation of error variance using Ensemble method which combines bootstraping and sampling with srswor in HSIC LASSO.
}
\usage{
hsic.var.ensemble(x,y,b,d)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{b}{number of bootstrap samples.}
  \item{d}{number of variables to be selected from x.}
}
\details{
In this method, both bootstrapping and simple random sampling without replacement are combined to estimate error variance. Variables are selected using HSIC LASSO from the original datasets and all possible samples of a particular size are taken from the selected variables set with simple random sampling without replacement. With these selected samples error variance is estimated from bootstrap samples of the original datasets using least squared regression method. Finally the average of all the estimated variances is considered as the final estimate of the error variance.
}
\value{
\item{Error variance}{}
}
\references{
Yamada, M., Jitkrittum, W., Sigal, L., Xing, E. P. and Sugiyama, M. (2014). High-Dimensional Feature Selection by Feature-Wise Kernelized Lasso. \emph{Neural Computation}, 26(1):185-207. doi:10.1162/NECO_a_00537
}
\author{
Sayanti Guha Majumdar <\email{sayanti23gm@gmail.com}>, Anil Rai, Dwijesh Chandra Mishra
}
\examples{
library(GSelection)
data(GS)
x_trn <- GS[1:40,1:110]
y_trn <- GS[1:40,111]
x_tst <- GS[41:60,1:110]
y_tst <- GS[41:60,111]
hsic_var <- hsic.var.ensemble(x_trn,y_trn,2,10)
}
\keyword{Ensemble}
\keyword{HSIC LASSO}
