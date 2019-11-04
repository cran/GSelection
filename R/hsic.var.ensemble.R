#' @title Error Variance Estimation in Genomic Prediction
#'
#' @description Estimation of error variance using Ensemble method in HSIC LASSO.
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

hsic.var.ensemble <- function(x,y,b,d){
p<- ncol(x)
n<- nrow(x)

new_x <- t(x)
new_y <- t(y)


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

##kernelized y
kernelLasso_y<-function(xin,res,lambda){
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
##R equivalent of repmat (matlab)
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


#Transformation of output
    y = yin;
    L = kernelGausian(y,y,1.0);
tmp = H%*%L%*%H;
LH = c(tmp);
return(LH)
}


##kernelized x
kernelLasso_x<-function(xin,res,lambda){
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
##R equivalent of repmat (matlab)
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
return(KH)
}


hsic_result_n1 <- kernelLasso(new_x,new_y,1)
coef_hsic_n1 <- coef(hsic_result_n1, "all")
coef_hsic_final_n1 <- coef_hsic_n1[-1]
beta_hsic_n1<-row(as.matrix(coef_hsic_final_n1))[which(!coef_hsic_final_n1==0)]
select_hsicbeta_n1 <- head(order(abs(coef_hsic_final_n1),decreasing = TRUE),d)
srswor_index_hsic <- combn(select_hsicbeta_n1,floor(d/2),FUN = NULL, simplify = FALSE)
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
selected_x[[j]] <- xb[[i]][,srswor_index_hsic[[j]]]
LH_2 <- kernelLasso_y(t(selected_x[[j]]),t(y),1)
KH_2 <- kernelLasso_x(t(selected_x[[j]]),t(y),1)
fit_x2 <- lm(LH_2 ~ KH_2 - 1)
var[[j]] <- sum((fit_x2$resid)^2)/(n - M1)

}
final_var[[i]] <- var
}


final_hsic_var2 <- numeric()
final_hsic_var1<- vector("list",b)
for (i in 1:b){
final_hsic_var<- numeric()
for (j in 1:s){
final_hsic_var[[j]]<- mean(final_var[[i]][[j]])
}
final_hsic_var1[[i]] <- final_hsic_var
final_hsic_var2[[i]] <- mean(final_hsic_var1[[i]])
}
final_hsic_var3 <- mean(final_hsic_var2)
return(final_hsic_var3)
}
