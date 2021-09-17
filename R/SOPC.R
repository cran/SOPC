#' The sparse online principal component can not only process the real-time updated data set and stream data set, but also obtain the sparse solution of the updated data set.
#'
#' @param data is a highly correlated online data set 
#' @param m is the number of principal component 
#' @param gamma is a sparse parameter 
#' @param eta is the proportion of online data to total data 
#'
#' @return Aso,Dso
#' @export
#'
#' @examples
#' SOPC(data=PSA,m=3,gamma=0.03,eta=0.6)
SOPC<-function(data,m,gamma,eta){
X<-scale(data)             
S<-cov(X)
eig<-eigen(S)
p<-nrow(S)
n<-nrow(X)                                                
n0<-round(eta*n)          
Xbar<-colMeans(X[1:n0,])   
S0<-cov(X[1:n0,])                             
lambda=eigen(S0)$values[1:m]  
V<-eigen(S0)$vectors[,1:m]                    
paras<-rep(gamma,1*m,m)       
Sd<-S0

iter1=0
for (i in (n0+1):n) {
iter1=iter1+1
Xcenter<-t(X[i,])           
Sd<-((i-1)/i)*Sd+(1/i)*t(Xcenter)%*%Xcenter 
lambda<-eigen(Sd)$values  
spc<-spca(Sd,K=m,type="Gram",max.iter=5,sparse="penalty",trace=FALSE,para=paras) 
V<-spc$loadings    
}
lambda2<-lambda[1:m]
V2<-V[,1:m]
Aso<-matrix(0,nrow=p,ncol=m)
for (j in 1:m){
  Aso[,j]<-sqrt(lambda2[j])*V2[,j]
  } 
Aso;table(Aso==0)
h2<-diag(Aso%*%t(Aso))
Dso<-diag(S-h2)
return(list(Aso=Aso,Dso=Dso))
}


