\name{Hugging}
\alias{Hugging}
\docType{data}
\title{
Hugging
}
\description{
The EMG Physical Action-Hugging data set.
}
\usage{data("Hugging")}
\format{
A data frame with 9752 observations on the following 8 variables. 
  \describe{
    \item{\code{A}}{a numeric vector}
    \item{\code{B}}{a numeric vector}
    \item{\code{C}}{a numeric vector}
    \item{\code{D}}{a numeric vector}
    \item{\code{E}}{a numeric vector}
    \item{\code{F}}{a numeric vector}
    \item{\code{G}}{a numeric vector}
    \item{\code{H}}{a numeric vector}
  }
}
\details{
The data set is a body movement data set, including 10 normal and 10 aggressive body movements. The data frame with 9752 observations on the following 8 variables.
}
\source{
The Hugging data set comes from the UCI database.
}
\references{
Demir et al. (2019). Surface emg signals and deep transfer learning-based physical action classification. Neural Computing and Applications.
}
\examples{
data(Hugging)
## maybe str(Hugging) ; plot(Hugging) ...
}
\keyword{datasets}