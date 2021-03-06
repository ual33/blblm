% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{confint.blbglm_parallel}
\alias{confint.blbglm_parallel}
\title{confint method of blbglm_parallel}
\usage{
\method{confint}{blbglm_parallel}(object, parm = NULL, level = 0.95, cluster = 1, ...)
}
\arguments{
\item{object}{a blblm_parallel fitted class}

\item{parm}{specify parameter, otherwise all X variables in the formula}

\item{level}{confidence level}

\item{cluster}{for parallel computing if cluster > 1}

\item{...}{others}
}
\description{
confidence interval of blbglm_parallel
}
\examples{
#proj = blblm::blbglm_parallel(vs~mpg+disp, mtcars, m = 3, B = 100, cluster = 2)
#confint(proj, level = 0.95, parm = NULL, cluster = 2)
}
