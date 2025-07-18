% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barycenter.R
\name{gaussian_wd}
\alias{gaussian_wd}
\title{Compute the Wasserstain distance between two GMMs' components on the variance part}
\usage{
gaussian_wd(V1, V2)
}
\arguments{
\item{V1}{First GMM variance part}

\item{V2}{Second GMM variance part}
}
\value{
Returns matrix that contains the Wasserstain distance between two GMMs' components
(V1 -> V2 with rows correspond to components in the first GMM
and columns correspond to components in the second GMM)
}
\description{
Compute the Wasserstain distance between two GMMs' components on the variance part
}
