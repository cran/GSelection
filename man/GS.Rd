\name{GS}
\alias{GS}
\docType{data}
    \title{Genotypic and phenotypic simulated dataset}
    \format{A data frame with 60 rows as genotypes with 111 columns (i.e. contains information of genotyped markers and phenotypic traits).}
    \details{
It has total 60 rows which represents 200 individuals genotypes and a total of 111 of columns, in which first 110 columns contain information of genotyped markers and last column represents value of phenotypic trait associated with genotype under study.
}
    \source{
Yandell, B. S., Mehta, T., Banerjee, S., Shriner, D., Venkataraman, R. et al. (2007). R/qtlbim: QTL with Bayesian Interval Mapping in experimental crosses. \emph{Bioinformatics}, 23, 641-643.
\cr Yandell, B. S., Nengjun, Y., Mehta, T., Banerjee, S., Shriner, D. et al. (2012). qtlbim: QTL Bayesian Interval Mapping. R package version 2.0.5. http://CRAN.R-project.org/package=qtlbim
}
    \description{
     This dataset is simulated with the R package "qtlbim" where 10 are true features associated with the trait of study and remaining 100 are random markers. we consider 10 chromosomes each containing 10 markers. Each chromosome have 1 qtl which is the true feature.
}

    \usage{data("GS")}

\examples{
library(GSelection)
data(GS)
X<-GS[,1:110]   ## Extracting Genotype
Y<-GS[,111] ## Extracting Phenotype
}
