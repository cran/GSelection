\name{spam.var.ensemble}
\alias{spam.var.ensemble}
\title{
Error Variance Estimation in Genomic Prediction
}
\description{
Estimation of error variance using Ensemble method which combines bootstraping and sampling with srswor in Sparse Additive Models.
}
\usage{
spam.var.ensemble(x,y,b,d)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{b}{number of bootstrap samples}
  \item{d}{number of variables to be selected from x.}
}

\details{
In this method, both bootstrapping and simple random sampling without replacement are combined to estimate error variance. Variables are selected using Sparse Additive Models (SpAM) from the original datasets and all possible samples of a particular size are taken from the selected variables set with simple random sampling without replacement. With these selected samples error variance is estimated from bootstrap samples of the original datasets using least squared regression method. Finally the average of all the estimated variances is considered as the final estimate of the error variance.
}
\value{
\item{Error variance}{}
}
\references{
Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030. doi:10.1111/j.1467-9868.2009.00718.x
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
spam_var <- spam.var.ensemble(x_trn,y_trn,2,10)
}
\keyword{Ensemble}
\keyword{SpAM}
