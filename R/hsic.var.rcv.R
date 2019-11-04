#' @title Error Variance Estimation in Genomic Prediction
#'
#' @description Estimation of error variance using Refitted cross validation in HSIC LASSO.
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

hsic.var.rcv <- function(x,y,d){
p<- ncol(x)
n<- nrow(x)

k <- floor(n/2)
x1 <- x[1:k, ]
y1 <- y[1:k]
x2 <- x[(k + 1):n, ]
y2 <- y[(k + 1):n]
n1 <- k
n2 <- n-k


new_x1 <- t(x1)
new_x2 <- t(x2)
new_y1 <- t(y1)
new_y2 <- t(y2)



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
return(KH)
}


hsic_result_n1 <- kernelLasso(new_x1,new_y1,1)
coef_hsic_n1 <- coef(hsic_result_n1, "all")
coef_hsic_final_n1 <- coef_hsic_n1[-1]
beta_hsic_n1<-row(as.matrix(coef_hsic_final_n1))[which(!coef_hsic_final_n1==0)]
select_hsicbeta_n1 <- head(order(abs(coef_hsic_final_n1),decreasing = TRUE),d)
M1 <- length(select_hsicbeta_n1)
selected_x2 <- x2[,select_hsicbeta_n1]
new_selected_x2 <- t(selected_x2)


LH_2 <- kernelLasso_y(new_selected_x2,new_y2,1)
KH_2 <- kernelLasso_x(new_selected_x2,new_y2,1)
fit_x2 <- lm(LH_2 ~ KH_2 - 1)
var1 <- sum((fit_x2$resid)^2)/(n - k - M1)



hsic_result_n2 <- kernelLasso(new_x2,new_y2,1)
coef_hsic_n2 <- coef(hsic_result_n2, "all")
coef_hsic_final_n2 <- coef_hsic_n2[-1]
beta_hsic_n2<-row(as.matrix(coef_hsic_final_n2))[which(!coef_hsic_final_n2==0)]
select_hsicbeta_n2 <- head(order(abs(coef_hsic_final_n2),decreasing = TRUE),d)
M2 <- length(select_hsicbeta_n2)
selected_x1 <- x1[,select_hsicbeta_n2]
new_selected_x1 <- t(selected_x1)


LH_1 <- kernelLasso_y(new_selected_x1,new_y1,1)
KH_1 <- kernelLasso_x(new_selected_x1,new_y1,1)
fit_x1 <- lm(LH_1 ~ KH_1 - 1)
var2 <- sum((fit_x1$resid)^2)/(k - M2)


var_rcv <- (var1 + var2)/2
return(var_rcv)
}
