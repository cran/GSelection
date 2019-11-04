\name{hsic.var.rcv}
\alias{hsic.var.rcv}
\title{
Error Variance Estimation in Genomic Prediction
}
\description{
Estimation of error variance using Refitted Cross Validation in HSIC LASSO.
}
\usage{
hsic.var.rcv(x,y,d)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{d}{number of variables to be selected from x.}
}
\details{
Refitted cross validation method (RCV) which is a two step method, is used to get the estimate of the error variance. In first step, dataset is divided into two sub-datasets and with the help of HSIC LASSO most significant markers(variables) are selected from the two sub-datasets. This results in two small sets of selected variables. Then using the set selected from 1st sub-dataset error variance is estimated from the 2nd sub-dataset with ordinary least square method and using the set selected from the 2nd sub-dataset error variance is estimated from the 1st sub-dataset with ordinary least square method. Finally the average of those two error variances are taken as the final estimator of error variance with RCV method.
}
\value{
\item{Error variance}{}
}
\references{Fan, J., Guo, S., Hao, N. (2012). Variance estimation using refitted cross-validation in ultrahigh dimensional regression. \emph{Journal of the Royal Statistical Society}, 74(1), 37-65.
\cr Yamada, M., Jitkrittum, W., Sigal, L., Xing, E. P. and Sugiyama, M. (2014). High-Dimensional Feature Selection by Feature-Wise Kernelized Lasso. \emph{Neural Computation}, 26(1):185-207. doi:10.1162/NECO_a_00537
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
hsic_var <- hsic.var.rcv(x_trn,y_trn,10)
}
\keyword{RCV}
\keyword{HSIC LASSO}
