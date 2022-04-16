#' The perturbation principal component can handle online data sets with highly correlated.
#'
#' @param data is a highly correlated online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#'
#' @return Ap,Dp
#' @export
#'
#' @examples
#' PPC(data=PSA,m=3,eta=0.8) 
PPC<-function(data,m,eta){
X<-as.matrix(scale(data))
S<-cov(X)   
n<-nrow(X)          
n0<-round(eta*n)  
p<-ncol(X)           
Xbar<-colMeans(X[1:n0,])      
eig1<-eigen(cov(X[1:n0,]-Xbar)) 
lambda<-eig1$values    
V<-eig1$vectors    
V1<-V[,1:m]

for (i in (n0+1):n) {
f<-1/i 
Xcenter<-t(X[i,])                    
lambda<-(1-f)*lambda       
Q<-sqrt(f)*t((X[i,]-Xbar)%*%V)
Q2<-Q*Q                      
num<-tcrossprod(Q)              
den<-matrix(lambda+Q2,p,p,byrow=T)-matrix(Q2+lambda,p,p)  # 
U<-num/den       
diag(U)<-1        
V<-V%*%U      
sigma2<-.colSums(V*V,p,p)
lambda<-lambda*sigma2
V<-V*rep.int(1/sqrt(sigma2),rep.int(p,p))
Xbar<-((i-1)/i)*Xbar+(1/i)*X[i,]   
}
ind<-order(lambda,decreasing=T)    
lambda<-lambda[ind]
V<-V[,ind]
V2<-V[,1:m]

Ap<-matrix(0,nrow=p,ncol=m)
for (j in 1:m){
  Ap[,j]<-sqrt(lambda[j])*V2[,j]
  } 
h2<-diag(Ap%*%t(Ap))
Dp<-diag(S-h2)  
return(list(Ap=Ap,Dp=Dp))
}


