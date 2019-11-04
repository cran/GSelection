#' @title Redundancy Rate
#'
#' @description Calculate the redundancy rate of the selected features(markers). Value will be high if many redundant features are selected.
#'
#' @param x
#'
#' @param spam_selected_feature_index
#'
#' @param hsic_selected_feature_index
#'
#' @param integrated_selected_feature_index
#'
#' @return Redundancy score
#'
#' @examples
#'
#' @export

RED <- function(x,spam_selected_feature_index,hsic_selected_feature_index, integrated_selected_feature_index){

spam_selected_feature = spam_selected_feature_index
select_hsicbeta = hsic_selected_feature_index
indc = integrated_selected_feature_index
requireNamespace("gdata")

cor_spam <- cor(x[,spam_selected_feature],x[,spam_selected_feature], method="spearman")
cor_spam_final <- upperTriangle(cor_spam, diag=FALSE)
m_spam<- length(spam_selected_feature)
RED_spam<- ((sum(abs(cor_spam_final)))/(m_spam*(m_spam-1)))

cor_hsic <- cor(x[,select_hsicbeta], x[,select_hsicbeta],method="spearman")
cor_hsic_final <- upperTriangle(cor_hsic, diag=FALSE)
m_hsic <- length(select_hsicbeta)
RED_hsic<- ((sum(abs(cor_hsic_final)))/(m_hsic*(m_hsic-1)))

cor_selectedfeature_I <- cor(x[,indc],x[,indc],method="spearman")
cor_selectedfeature_I_final <- upperTriangle(cor_selectedfeature_I, diag=FALSE)
m_I<- length(indc)
RED_I<- ((sum(abs(cor_selectedfeature_I_final)))/(m_I*(m_I-1)))
result <- list(RED_spam=RED_spam,RED_hsic=RED_hsic,RED_I=RED_I)
return(result)
}
