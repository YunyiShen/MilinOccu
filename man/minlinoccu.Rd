% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MilinOccu.R
\name{minlinoccu}
\alias{minlinoccu}
\title{Single Season Occupancy Model with Min-Linear Logistic Regression}
\usage{
minlinoccu(Y, Designs_occu, Designs_det, ...)
}
\arguments{
\item{Y}{The detection history matrix, with row as sites and column as repeats, missing data should be NA}

\item{Designs_occu}{A list with each element is a design matrix without intercept, different element can have overlap columns but should not be identical}

\item{Designs_det}{A list of design matrix list. Each element represent the design matrix list of repeat, and the design matrix list should have the same structure like Design_occu. Structure should be consistant among repeats.}

\item{...}{Extra arguments send to optim}
}
\value{
A list,`$formatted` is the formatted parameter estimation and `$optim` is the raw output of optim.
}
\description{
Fit a single season occupancy model based on min-linear Logistic regression.
}
\details{
Instead of logistic regression, this model fits a min-linear logistic based occupancy, i.e. logit(psi) and logit(p) are not linear to predictors but is the minimum among several linear predictors, which represent the worst environment.
}
