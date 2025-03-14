% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Helpers.R
\name{CATE_score}
\alias{CATE_score}
\title{Evaluates the Neyman-Orthogonal Score Function for the CATE
at a given point}
\usage{
CATE_score(X, mu_0, mu_1, Y, W, pi)
}
\arguments{
\item{X}{Covariates}

\item{mu_0}{(Estimated) value of the untreated regression function at X}

\item{mu_1}{(Estimated) value of the treated regression function at X}

\item{Y}{Response}

\item{W}{Treatment Status (coded as 0 for untreated and 1 for treated)}

\item{pi}{(Estimated) Propensity Score}
}
\description{
Evaluates the Neyman-Orthogonal Score Function for the CATE
at a given point
}
