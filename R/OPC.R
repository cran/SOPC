#' The online principal components can handle data sets that are updated in real time and streaming data.
#'
#' @param data is a highly correlated online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#'
#' @return Ao,Do
#' @export
#'
#' @examples
#' OPC(data=PSA,m=3,eta=0.8) 
OPC<-function(data,m,eta){
X<-as.matrix(scale(data))
S<-cov(X)   
n<-nrow(X)          
n0<-round(eta*n)  
p<-ncol(X)           
Xbar<-colMeans(X[1:n0,])      
eig1<-eigen(cov(X[1:n0,]-Xbar)) 
lambda<-eig1$values[1:m]       
V<-eig1$vectors[,1:m]      
V1<-V
T<-matrix(rep(0,(m+1)*(m+1)),nrow=(m+1))

for (i in (n0+1):n) {
Xcenter<-t(X[i,]-Xbar)       
g<-t(V)%*%t(Xcenter)  
Xhat<-t(V%*%g)+Xbar          
h<-t(X[i,]-Xhat)            
hmao<-norm(h,"2")            
gamma<-as.numeric(t(h/hmao)%*%t(Xcenter))   
T[1:m,]<-cbind(((i-1)/i)*diag(lambda)+((i-1)^2/i^3)*g%*%t(g),((i-1)^2/i^3)*gamma*g)
T[(m+1),]<-cbind(((i-1)^2/i^3)*gamma*t(g),((i-1)^2/i^3)*gamma^2)               
eig2<-eigen(T)                
lambda<-eig2$values[1:m]               
V<-(cbind(V,h/hmao)%*%eig2$vectors)[,1:m]
Xbar<-((i-1)/i)*Xbar+(1/i)*X[i,]   
}
V2<-V[,1:m]
Ao<-matrix(0,nrow=p,ncol=m)
for (j in 1:m){
  Ao[,j]<-sqrt(lambda[j])*V2[,j]
  } 
h2<-diag(Ao%*%t(Ao))
Do<-diag(S-h2) 
return(list(Ao=Ao,Do=Do))
}


