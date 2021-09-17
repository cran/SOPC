#' The sparse principal component can obtain sparse solutions of the eigenmatrix to better explain the relationship between principal components and original variables.
#'
#' @param data is a set of highly correlated variables 
#' @param m is the number of principal component
#' @param gamma is a sparse parameter
#'
#' @return As,Ds
#' @export
#'
#' @examples
#' SPC(data=PSA,m=3,gamma=0.03)
SPC<-function(data,m,gamma){
X<-scale(data)
R<-cor(X)
S<-R
eig<-eigen(S)
p<-nrow(S)
n<-nrow(X) 
paras<-rep(gamma,1*m,m) 

spc<-spca(R,K=m,type="Gram",sparse="penalty",trace=FALSE,para=paras)
lambda<-eig$values[1:m]
V<-spc$loadings[,1:m]
As<-matrix(0,nrow=p,ncol=m)
for (j in 1:m){
  As[,j]<-sqrt(lambda[j])*V[,j]
} 
As;table(As==0)
h2<-diag(As%*%t(As))
Ds<-diag(S-h2)
return(list(As=As,Ds=Ds))
}


