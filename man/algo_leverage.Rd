% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{algo_leverage}
\alias{algo_leverage}
\title{Implements algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows}
\usage{
algo_leverage(y, X, r, plot = F, method = "unif")
}
\arguments{
\item{y}{the response variable matrix}

\item{X}{the predictor variable matrix}

\item{r}{the number of subsamples to choose from}

\item{plot}{if true, a plot of points selected will be shown}

\item{method}{either "unif" or "lev", determines wheter uniform or leverage score will be used}
}
\value{
the linear regression coefficient on the subsample
}
\description{
Implements algorithmic leveraging for linear regression using uniform and leverage score based subsampling of rows
}
\examples{
X = matrix(rt(500, 6), nrow=500)
ep = matrix(rnorm(500, 0, 1), nrow=500)
y = -X + ep
r = 100
algo_leverage(y, X, r, method="unif")
algo_leverage(y, X, r, method="lev")
}
