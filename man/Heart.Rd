\name{Heart}
\alias{Heart}
\docType{data}
\title{
Heart failure
}
\description{
Heart failure
}
\usage{data("Heart")}
\format{
  A data frame with 299 observations on the following 13 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{anaemia}}{a numeric vector}
    \item{\code{creatinine_phosphokinase}}{a numeric vector}
    \item{\code{diabetes}}{a numeric vector}
    \item{\code{ejection_fraction}}{a numeric vector}
    \item{\code{high_blood_pressure}}{a numeric vector}
    \item{\code{platelets}}{a numeric vector}
    \item{\code{serum_creatinine}}{a numeric vector}
    \item{\code{serum_sodium}}{a numeric vector}
    \item{\code{sex}}{a numeric vector}
    \item{\code{smoking}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{DEATH_EVENT}}{a numeric vector}
  }
}
\details{
This dataset contains the medical records of 299 patients who had heart failure, collected during their follow-up period, where each patient profile has 13 clinical features.
}
\source{
The Heart failure data set comes from the UCI database.
}
\references{
Davide Chicco, Giuseppe Jurman. (2020). Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. BMC Medical Informatics and Decision Making.
}
\examples{
data(Heart)
## maybe str(Heart) ; plot(Heart) ...
}
\keyword{datasets}
