#' The traditional principal component method. This method can estimate the eigen space of the data set.
#'
#' @param data is a set of highly correlated variables
#' @param m is the number of principal component
#'
#' @return Ahat, Dhat
#' @export
#'
#' @examples 
#' PC(data=PSA,m=3)
PC<-function(data,m=m){
X<-scale(data)
R<-cor(X)
S<-R
eig<-eigen(S)
p<-nrow(S)
diag_S<-diag(S)
sum_rank<-sum(diag_S)
rowname<-paste("X",1:p,sep="")
colname<-paste("Factor",1:m,sep="")
Ahat<-matrix(0,nrow=p,ncol=m,
dimnames=list(rowname,colname))
rowname<-c("SS loadings","Proportion Var","Cumulative Var")
for (i in 1:m){
  Ahat[,i]<-sqrt(eig$values[i])*eig$vectors[,i]
  } 
h2<-diag(Ahat%*%t(Ahat))      
Dhat<-diag(S-h2)
return(list(Ahat=Ahat,Dhat=Dhat))
}

