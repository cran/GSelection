\name{RED}
\alias{RED}
\title{
Redundancy Rate
}
\description{
Calculate the redundancy rate of the selected features(markers). Value will be high if many redundant features are selected.
}
\usage{
RED(x,spam_selected_feature_index,hsic_selected_feature_index,
integrated_selected_feature_index)
}
\arguments{
  \item{x}{a matrix of markers or explanatory variables, each column contains one marker and each row represents an individual.}
  \item{spam_selected_feature_index}{index of selected markers from x using Sparse Additive Model.}
  \item{hsic_selected_feature_index}{index of selected markers from x using HSIC LASSO.}
  \item{integrated_selected_feature_index}{index of selected markers from x using integrated model framework}
}
\details{
The RED score (Zhao et al., 2010) is determined by average of the correlation between each pair of selected markers. A large RED score signifies that selected features are more strongly correlated to each other which means many redundant features are selected. Thus, a small redundancy rate is preferable for feature selection.
}
\value{
Returns a LIST containing
  \item{RED_spam}{returns redundancy rate of features selected by using Sparse Additive Model.}
  \item{RED_hsic}{returns redundancy rate of features selected by using HSIC LASSO.}
  \item{RED_I}{returns redundancy rate of features selected by using integrated model framework.}
}
\references{Guha Majumdar, S., Rai, A. and Mishra, D. C. (2019). Integrated framework for selection of additive and non-additive genetic markers for genomic selection. \emph{Journal of Computational Biology}. doi:10.1089/cmb.2019.0223
\cr Ravikumar, P., Lafferty, J., Liu, H. and Wasserman, L. (2009). Sparse additive models. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 71(5), 1009-1030. doi:10.1111/j.1467-9868.2009.00718.x
\cr Yamada, M., Jitkrittum, W., Sigal, L., Xing, E. P. and Sugiyama, M. (2014). High-Dimensional Feature Selection by Feature-Wise Kernelized Lasso. \emph{Neural Computation}, 26(1):185-207. doi:10.1162/NECO_a_00537
\cr Zhao, Z., Wang, L. and Li, H. (2010). Efficient spectral feature selection with minimum redundancy. \emph{In AAAI Conference on Artificial Intelligence (AAAI)}, pp 673-678.
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
red <- RED(x_trn,fit$spam_selected_feature_index,fit$hsic_selected_feature_index,
fit$integrated_selected_feature_index)
}
\keyword{ Redundancy rate }
\keyword{ Integrated Model }
