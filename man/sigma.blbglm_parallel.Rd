% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{sigma.blbglm_parallel}
\alias{sigma.blbglm_parallel}
\title{MSE of a blbglm_parallel object}
\usage{
\method{sigma}{blbglm_parallel}(object, confidence = FALSE, level = 0.95, cluster = 1, ...)
}
\arguments{
\item{object}{a fitted blblm_parallel class}

\item{confidence}{bool indicating if you want a confidence interval of sigma}

\item{level}{confidence level}

\item{cluster}{number of clusters to use for parallel computing. when = 1 will be running on single thread.}

\item{...}{others}
}
\description{
coefficients method of blbglm_parallel class returning coefficients
}
\examples{
#proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#sigma(proj, confidence = TRUE, level = 0.99, cluster = 2)
}
