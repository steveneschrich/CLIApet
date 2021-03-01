#' Create a correlation coefficient heatmap figure.
#'
#' This function creates a heatmap using ComplexHeatmap to show the correlation
#' among samples in a matrix (assuming samples are columns in the matrix). This
#' is a quick way to show the relationship between samples.
#'
#' Of note, the color map is created such that the colors are from white to red, scaled
#' to the range of correlations seen. This is because when doing sample-sample correlations,
#' it's often the case that they are all highly correlated. Rather than over-emphasize
#' differences that aren't there, this is a cleaner approach.
#'
#' @param x A matrix of expression data.
#'
#' @return A ComplexHeatmap of correlation.
#' @export
correlation_coefficient_heatmap<-function(x) {
  assertthat::assert_that(methods::is(x, "matrix"))
  correlation_matrix<-stats::cor(x)
  col_fun<-circlize::colorRamp2(c(min(correlation_matrix), max(correlation_matrix)), c("white", "red"))

  ComplexHeatmap::Heatmap(correlation_matrix,
                          heatmap_legend_param = list(title = "Correlation\nCoefficient"),
                          col=col_fun
  )

}
