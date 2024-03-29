% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPC.R
\name{PPC}
\alias{PPC}
\title{The perturbation principal component can handle online data sets with highly correlated.}
\usage{
PPC(data, m, eta)
}
\arguments{
\item{data}{is a highly correlated online data set}

\item{m}{is the number of principal component}

\item{eta}{is the proportion of online data to total data}
}
\value{
Ap,Dp
}
\description{
The perturbation principal component can handle online data sets with highly correlated.
}
\examples{
PPC(data=PSA,m=3,eta=0.8) 
}
