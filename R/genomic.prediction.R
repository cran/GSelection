#' @title Genomic Prediction
#'
#' @description Prediction of phenotypic values based on selected markers with integrated model framework using both additive (Sparse Additive Models) and non-additive (HSIC LASSO) statistical models.
#'
#' @param x
#'
#' @param spam_error_var
#'
#' @param hsic_error_var
#'
#' @param spam_selected_feature_index
#'
#' @param hsic_selected_feature_index
#'
#' @param coefficient.spam
#'
#' @param coefficient.hsic
#'
#' @return Selected Features(markers)
#'
#' @examples
#'
#' @export

genomic.prediction <- function(x,spam_error_var,hsic_error_var,spam_selected_feature_index,hsic_selected_feature_index,coefficient.spam,coefficient.hsic){
y_spam <- x[,spam_selected_feature_index] %*% coefficient.spam
y_hsic <- x[,hsic_selected_feature_index] %*% coefficient.hsic
Integrated_y <- ((hsic_error_var/(hsic_error_var+spam_error_var))*y_spam)+((spam_error_var/(hsic_error_var+spam_error_var))*y_hsic)
return(Integrated_y)
}
