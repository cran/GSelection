\name{genomic.prediction}
\alias{genomic.prediction}
\title{
Genomic Prediction
}
\description{
Prediction of phenotypic values based on selected markers with integrated model framework using both additive (Sparse Additive Models) and non-additive (HSIC LASSO) statistical models.
}
\usage{
genomic.prediction(x,spam_error_var,hsic_error_var,
spam_selected_feature_index,hsic_selected_feature_index,
coefficient.spam,coefficient.hsic)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables for which phenotype will be predicted. Each column contains one marker and each row represents an individual.}
  \item{spam_error_var}{estimated error variance of genomic prediction by Sparse Additive Model.}
  \item{hsic_error_var}{estimated error variance of genomic prediction by HSIC LASSO.}
  \item{spam_selected_feature_index}{index of selected markers from x using Sparse Additive Model}
  \item{hsic_selected_feature_index}{index of selected markers from x using HSIC LASSO.}
  \item{coefficient.spam}{coefficient values of selected markers using Sparse Additive Model.}
  \item{coefficient.hsic}{coefficient values of selected markers using HSIC LASSO.}
}
\details{
Phenotypic values will be predicted for given genotype of markers by using previously fitted model object. Integrated model framework is used for this purpose which is developed by combining selected features from SpAm and HSIC LASSO.
}
\value{
\item{Integrated_y}{returns predicted phenotype}
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

## estimate spam_var from function spam.var.ensemble or spam.var.rcv
spam_var <- 2.681972
## estimate hsic_var from function hsic.var.ensemble or hsic.var.rcv
hsic_var <- 10.36974

fit <- feature.selection(x_trn,y_trn,d=10)
pred_y <- genomic.prediction(x_tst,spam_var,hsic_var,
fit$spam_selected_feature_index,fit$hsic_selected_feature_index,
fit$coefficient.spam,fit$coefficient.hsic)
}
\keyword{SpAM}
\keyword{HSIC LASSO}
\keyword{Integrated Model}
\keyword{Genomic Prediction}
