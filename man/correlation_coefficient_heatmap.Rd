% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation_coefficient_heatmap.R
\name{correlation_coefficient_heatmap}
\alias{correlation_coefficient_heatmap}
\title{Create a correlation coefficient heatmap figure.}
\usage{
correlation_coefficient_heatmap(x)
}
\arguments{
\item{x}{A matrix of expression data.}
}
\value{
A ComplexHeatmap of correlation.
}
\description{
This function creates a heatmap using ComplexHeatmap to show the correlation
among samples in a matrix (assuming samples are columns in the matrix). This
is a quick way to show the relationship between samples.
}
\details{
Of note, the color map is created such that the colors are from white to red, scaled
to the range of correlations seen. This is because when doing sample-sample correlations,
it's often the case that they are all highly correlated. Rather than over-emphasize
differences that aren't there, this is a cleaner approach.
}
