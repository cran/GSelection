#' @title Error Variance Estimation in Genomic Prediction
#'
#' @description Estimation of error variance using Refitted cross validation in Sparse Additive Models.
#'
#' @param x
#'
#' @param y
#'
#' @param d
#'
#' @return Error variance
#'
#' @examples
#'
#' @export

spam.var.rcv <- function(x,y,d){
  p<- ncol(x)
  n<- nrow(x)

  k <- floor(n/2)
  x1 <- x[1:k, ]
  y1 <- y[1:k]
  x2 <- x[(k + 1):n, ]
  y2 <- y[(k + 1):n]
  n1 <- k
  n2 <- n-k

  requireNamespace("SAM")
  spam_fit_n1 <- samQL(x1,y1,p=1)
  w_n1 <- row(as.matrix(spam_fit_n1$w[,30]))[which(spam_fit_n1$w[,30] != 0)]
  w_order_n1 <- head(order(spam_fit_n1$w[,30],decreasing = TRUE),d)
  w_value_n1 <- as.matrix(spam_fit_n1$w[,30])[w_order_n1,]
  spam_selected_feature_n1<- w_order_n1

  M1 <- length(spam_selected_feature_n1)
  selected_x2 <- x2[,spam_selected_feature_n1]

  fit_x2 <- lm(y2 ~ selected_x2 - 1)
  var1 <- sum((fit_x2$resid)^2)/(n - k - M1)



  spam_fit_n2 <- samQL(x2,y2,p=1)
  w_n2 <- row(as.matrix(spam_fit_n2$w[,30]))[which(spam_fit_n2$w[,30] != 0)]
  w_order_n2 <- head(order(spam_fit_n2$w[,30],decreasing = TRUE),d)
  w_value_n2 <- as.matrix(spam_fit_n2$w[,30])[w_order_n2,]
  spam_selected_feature_n2<- w_order_n2

  M2 <- length(spam_selected_feature_n2)
  selected_x1 <- x1[,spam_selected_feature_n2]

  fit_x1 <- lm(y1 ~ selected_x1 - 1)
  var2 <- sum((fit_x1$resid)^2)/(k - M2)
  var_rcv <- (var1 + var2)/2
  return(var_rcv)
}
