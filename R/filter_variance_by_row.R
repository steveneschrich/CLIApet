#' Filter rows to those having high enough variance.
#'
#' The ExpressionSet is filtered such that only genes
#' with a variance larger than a specific value (given
#' by threshold) are kept.
#'
#' @param x An ExpressionSet.
#' @param threshold Variance threshold to use.
#'
#' @return A matrix with filtered variance is returned.
#' @export
#'
filter_variance_by_row<-function(x, threshold=0) {
  # Filter out genes (rows) by low variance
  x_rows_with_variance<-which(apply(Biobase::exprs(x),1,stats::var)>0)
  x_filtered<-x

  Biobase::exprs(x)[x_rows_with_variance,]
}
