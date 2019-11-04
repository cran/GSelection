\name{feature.selection}
\alias{feature.selection}
\title{
Genomic Feature Selection
}
\description{
Feature (marker) selection in case of genomic prediction with integrated model framework using both additive (Sparse Additive Models) and non-additive (HSIC LASSO) statistical models.
}
\usage{
feature.selection(x,y,d)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{y}{a column vector of response variable.}
  \item{d}{number of variables to be selected from x.}
}
\details{
Integrated model framework was developed by combining one additive model (Sparse Additive Model) and one non-additive model (HSIC LASSO) for selection of important markers from whole genome marker data.
}
\value{
Returns a LIST containing
  \item{spam_selected_feature_index}{returns index of selected markers from x using Sparse Additive Model}
  \item{coefficient.spam}{returns coefficient values of selected markers using Sparse Additive Model.}
  \item{hsic_selected_feature_index}{returns index of selected markers from x using HSIC LASSO.}
  \item{coefficient.hsic}{returns coefficient values of selected markers using HSIC LASSO.}
  \item{integrated_selected_feature_index}{returns index of selected markers from x using integrated model framework.}
}
\references{Guha Majumdar, S., Rai, A. and Mishra, D. C. (2019). Integrated framework for selection of additive and non-additive genetic markers for genomic selection. \emph{Journal of Computational Biology}. doi:10.1089/cmb.2019.0223
\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030. doi:10.1111/j.1467-9868.2009.00718.x
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
fit <- feature.selection(x_trn,y_trn,d=10)
}
\keyword{SpAM}
\keyword{HSIC LASSO}
\keyword{Integrated Model}
