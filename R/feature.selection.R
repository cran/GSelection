#' @title Genomic Feature Selection
#'
#' @description Feature (marker) selection in case of genomic prediction with integrated model framework using both additive (Sparse Additive Models) and non-additive (HSIC LASSO) statistical models.
#'
#' @param x
#'
#' @param y
#'
#' @param d
#'
#' @return Selected Features(markers)
#'
#' @examples
#'
#' @export

feature.selection <- function(x,y,d){

##SPAM
requireNamespace("SAM")
spam_fit <- samQL(x,y,p=1)
w <- row(as.matrix(spam_fit$w[,30]))[which(spam_fit$w[,30] != 0)]
w_order <- head(order(spam_fit$w[,30],decreasing = TRUE),d)
w_value <- as.matrix(spam_fit$w[,30])[w_order,]
spam_selected_feature<- w_order



##HSIC
kernelLasso<-function(xin,res,lambda){
xin<-as.matrix(xin)
yin<-as.matrix(res)
requireNamespace("penalized")
kernelGausian=function(x,c,sigma){
x=as.matrix(x)
c=as.matrix(c)
d=nrow(x)
nx=ncol(x)
nc=ncol(c)
x2=colSums(x^2,1)
x2=as.matrix(x2)
c2=colSums(c^2,1)
c2=as.matrix(c2)
repmat = function(y,m,n){

my = dim(y)[1]
ny = dim(y)[2]
matrix(t(matrix(y,my,ny*n)),my*m,ny*n,byrow=T)
}
distance2=repmat(t(c2),nx,1)+repmat(x2,1,nc)-2*t(x)%*%c
X=exp(-distance2/(2*sigma^2))
}
d<-nrow(xin)
n<-ncol(xin)

#Normalization
x=xin/((as.matrix(apply(xin,1,sd)))%*%(rep(1,ncol(xin))) +.Machine$double.eps)

#Centering matrix
H = diag(n) - 1/n*rep(1,n);

#Transformation of input
KH = matrix(0,n^2,d);
for (ii in 1:d){
    Kx = kernelGausian(x[ii,,drop=FALSE],x[ii,,drop=FALSE],1.0);
    tmp = H%*%Kx%*%H;
    KH[,ii] = c(tmp);
}
KH
#Transformation of output
    y = yin;
    L = kernelGausian(y,y,1.0);
tmp = H%*%L%*%H;
LH = c(tmp);
pen<-penalized(LH,KH,lambda1=lambda,positive=TRUE);
return(pen)
}

new_x_hsic <- t(x)
new_y_hsic <- t(y)
hsic_result <- kernelLasso(new_x_hsic,new_y_hsic,1)
coef_hsic <- coef(hsic_result, "all")
coef_hsic_final <- coef_hsic[-1]
coef_hsic_final <- as.matrix(coef_hsic_final)
beta_hsic <- row(coef_hsic_final)[which(!coef_hsic_final==0)]
select_hsicbeta <- head(order(abs(coef_hsic_final),decreasing = TRUE),d)
coefficient.hsic <- coef_hsic_final[select_hsicbeta,]


##Integrated model

##scaling of coefficient
sum_hsicbeta<- sum(coef_hsic_final[select_hsicbeta,])
scale_hsicbeta <- coef_hsic_final[select_hsicbeta,]/sum_hsicbeta
sum_spambeta <- sum(w_value)
scale_spambeta <- w_value/sum_spambeta

##combined selected feature for hsic and spam

I_selected_beta <- c(scale_spambeta,scale_hsicbeta)
I_index <-c(spam_selected_feature,select_hsicbeta)
I_mat <- cbind(I_index,I_selected_beta)
indc <- unique(I_mat[,1])
integrated_selected_feature_index <- indc


spam_selected_feature_index <- spam_selected_feature
hsic_selected_feature_index <- select_hsicbeta
result <- list(spam_selected_feature_index = spam_selected_feature_index, coefficient.spam = w_value,hsic_selected_feature_index = hsic_selected_feature_index,coefficient.hsic=coefficient.hsic,integrated_selected_feature_index = integrated_selected_feature_index)
return(result)
}



