
#' Create a PCA ggplot of gene expression data.
#'
#' @param x An ExpressionSet to plot.
#' @param main A title for the graph.
#'
#' @return A ggplot representing a PCA plot.
#' @export
#'
#' @importFrom rlang .data
pca_plot<-function(x, main="PCA Plot of samples") {
  x_matrix_filtered<-filter_variance_by_row(x)

  pca<-stats::prcomp(t(x_matrix_filtered),scale.=T)
  percent_var<-100*pca$sdev^2/sum(pca$sdev^2)
  plot.data<-data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],SampleNames=Biobase::sampleNames(x))
  ggplot2::ggplot(plot.data, ggplot2::aes(x=.data$PC1, y=.data$PC2, label=.data$SampleNames)) +
    ggplot2::geom_point(size=5,shape=21,fill="red") +
    ggrepel::geom_text_repel(ggplot2::aes(label=.data$SampleNames),hjust=-0.2) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(main) +
    ggplot2::xlab(sprintf("PC1 (%2.2f%%)",percent_var[1])) +
    ggplot2::ylab(sprintf("PC2 (%2.2f%%)",percent_var[2]))
}
