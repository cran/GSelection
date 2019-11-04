#' @title Error Variance Estimation in Genomic Prediction
#'
#' @description Estimation of error variance using Ensemble method in Sparse Additive Models.
#'
#' @param x
#'
#' @param y
#'
#' @param b
#'
#' @param d
#'
#' @return Error variance
#'
#' @examples
#'
#' @export

spam.var.ensemble <- function(x,y,b,d){
p<- ncol(x)
n<- nrow(x)


requireNamespace("SAM")
spam_fit_n1 <- samQL(x,y,p=1)
w_n1 <- row(as.matrix(spam_fit_n1$w[,30]))[which(spam_fit_n1$w[,30] != 0)]
w_order_n1 <- head(order(spam_fit_n1$w[,30],decreasing = TRUE),d)
w_value_n1 <- as.matrix(spam_fit_n1$w[,30])[w_order_n1,]
spam_selected_feature_n1<- w_order_n1
srswor_index <- combn(spam_selected_feature_n1,floor(d/2),FUN = NULL, simplify = FALSE)
s <- dim(combn(d,floor(d/2)))[2]
M1 <- dim(combn(d,floor(d/2)))[1]


final_var <- vector("list",b)
bs<- vector("list",b)
    xb <- vector("list",b)
    for (i in 1:b){
      bs[[i]] <- sample(1:n,n, replace = TRUE)
      xb[[i]] <- x[bs[[i]],]
      var <- numeric()
      selected_x <- vector("list",s)
      for (j in 1:s){
        selected_x[[j]] <- xb[[i]][,srswor_index[[j]]]

        fit_x2 <- lm(y ~ selected_x[[j]] - 1)
        var[[j]] <- sum((fit_x2$resid)^2)/(n - M1)
      }
      final_var[[i]] <- var
    }

    final_spam_var2 <- numeric()
    final_spam_var1<- vector("list",b)
    for (i in 1:b){
      final_spam_var<- numeric()
      for (j in 1:b){
        final_spam_var[[j]]<- mean(final_var[[i]][[j]])
      }
      final_spam_var1[[i]] <- final_spam_var
      final_spam_var2[[i]] <- mean(final_spam_var1[[i]])
    }

    final_spam_var3 <- mean(final_spam_var2)

    return(final_spam_var3)
  }
