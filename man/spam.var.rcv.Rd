\name{spam.var.rcv}
\alias{spam.var.rcv}
\title{
Error Variance Estimation in Genomic Prediction
}
\description{
Estimation of error variance using Refitted cross validation in Sparse Additive Models.
}
\usage{
spam.var.rcv(x,y,d)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{d}{number of variables to be selected from x.}
}
\details{
Refitted cross validation method (RCV) which is a two step method, is used to get the estimate of the error variance. In first step, dataset is divided into two sub-datasets and with the help of Sparse Additive Models (SpAM) most significant markers(variables) are selected from the two sub-datasets. This results in two small sets of selected variables. Then using the set selected from 1st sub-dataset error variance is estimated from the 2nd sub-dataset with ordinary least square method and using the set selected from the 2nd sub-dataset error variance is estimated from the 1st sub-dataset with ordinary least square method. Finally the average of those two error variances are taken as the final estimator of error variance with RCV method.
}
\value{
\item{Error variance}{}
}
\references{Fan, J., Guo, S., Hao, N. (2012).Variance estimation using refitted cross-validation in ultrahigh dimensional regression. \emph{Journal of the Royal Statistical Society}, 74(1), 37-65.
\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030. doi:10.1111/j.1467-9868.2009.00718.x}
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
spam_var <- spam.var.rcv(x_trn,y_trn,10)
}

\keyword{RCV}
\keyword{SpAM}
