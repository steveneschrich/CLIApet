#' Createa a ggplot-based MDS plot of gene expression data.
#'
#' @param x An ExpressionSet to plot.
#' @param main A title for the plot.
#'
#' @return A ggplot representing an MDS plot.
#' @export
#'
#' @importFrom rlang .data
#'
mds_plot<-function(x,main="MDS Plot of samples") {


  x_matrix_filtered<-filter_variance_by_row(x)

  mds<-stats::cmdscale(1-stats::cor(x_matrix_filtered))
  plot.data<-data.frame("MDS1"=mds[,1],"MDS2"=mds[,2],"SampleNames"=Biobase::sampleNames(x))
  ggplot2::ggplot(plot.data, ggplot2::aes(x=.data$MDS1, y=.data$MDS2, label=.data$SampleNames)) +
    ggplot2::geom_point(size=5,shape=21,fill="red") +
    ggrepel::geom_text_repel(ggplot2::aes(label=.data$SampleNames),hjust=-0.2, direction="both") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(main)

}
