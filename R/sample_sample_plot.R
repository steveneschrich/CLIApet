#' Pair a sample (n) with all the other samples in the matrix m.
#'
#' This function takes an integer n and a matrix m. The goal is to
#' generate all pairings of n with members of l (except n itself).
#' This is done in index space (e.g., if n=1, then it would provide
#' (1,2), (1,3), (1,4) for a 4-column matrix).
#'
#' I suspect there is a function in R to do this, but I couldn't find
#' it quickly.
#'
#' @param n The index of sample to pair everything against.
#' @param m A matrix of samples (one column is one sample).
#'
#' @return A matrix of columns representing pairings to consider.
#' @export
#'
pair_with<-function(n, m) {
  assertthat::assert_that(n>0, ncol(m)>0, n <= ncol(m))
  # Sequence is everything 1:l, except n
  s<-1:ncol(m)
  s<-s[-match(n, s)]

  matrix(data=c(rep(n, length(s)), s), byrow=T,ncol=length(s))
}



#' Create a sample-sample scatter plot.
#'
#' This is a plot type that does a scatter plot of two samples, that
#' are expected to be highly similar. Fold-change lines are included,
#' as are basic statistics on the degree of differences between the
#' two samples and degree of similarity (correlation).
#'
#' The assumption is that the plots are fairly similar, so there
#' should be a place to include text for metrics in the lower
#' right quadrant of the graph.
#'
#' @param x Sample 1
#' @param y Sample 2
#' @param fc Fold change lines
#' @param main Main title
#' @param xlab X label
#' @param ylab Y label
#'
#' @return ggplot of scatter plot.
#' @export
#'
sample_sample_plot<-function(x, y, fc=2.0,
                             main="Sample-Sample comparisons",
                             xlab="Sample 1",
                             ylab="Sample 2") {

  df <- data.frame(
    x=x, y=y
  )

  # Stats to compute
  cor_r<-stats::cor(x,y)
  diffs<-x-y
  diffs_1.5<-length(which(abs(diffs)>log2(1.5)))
  diffs_2<-length(which(abs(diffs)>log2(2)))

  xrange<-range(x)

  splot<-ggplot2::ggplot(df, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_point(fill="red",  shape=21) +
    ggplot2::theme_bw() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(main) +

    # Ensure the plot is square
    ggplot2::coord_fixed(ratio=1) +
    #  geom_text(x=range[2]+1,hjust=0,size=8)
    #theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin



    ggplot2::geom_abline(slope=1, intercept=0, color="blue", linetype=2) +
    ggplot2::geom_abline(slope=1, intercept=1, color="red", linetype=2) +
    ggplot2::geom_abline(slope=1, intercept=-1, color="red", linetype=2)

  #  stat_density_2d(aes(fill = ..density..,color="pink"), geom = 'raster', contour = FALSE)


  # Generate the text table next to the graph
  tplot<-ggpubr::ggtexttable(data.frame(labels=c("R","# diff at 1.5-fold","# diff at 2.0-fold"),
                                values=c(sprintf("%.5f",cor_r),
                                         sprintf("%d",diffs_1.5),
                                         sprintf("%d",diffs_2))
  ),
  rows=NULL,
  cols=NULL,
  theme = ggpubr::ttheme(
    base_size = 8,
    tbody.style = ggpubr::tbody_style(hjust=1, x=0.9))
  )

  # Arrange the scatter plot and text plot side-by-side, in a 3:1 ratio.
  ggpubr::ggarrange(splot,tplot,ncol=2, widths=c(3,1))
}



#' Generate a list of sample/sample density plots.
#'
#' This function is suitable for generating a bunch of sample-sample
#' plots that can then be ggarrange'd to plot out in bulk.
#'
#' @param x A matrix of samples (columns)
#' @param comparisons A matrix of sample pairings (columns)
#' @param main A title for the plot(s)
#'
#' @return A list of plots.
#' @export
#'
generate_sample_sample_density_plot_list<-function(x, comparisons=NULL, main="Sample-Sample comparisons") {
  assertthat::assert_that(ncol(x)>1, max(comparisons) <= ncol(x),
                          min(comparisons) > 0, !is.null(comparisons),
                          is.matrix(x))
  plot_list<-apply(comparisons, 2, function(n) {sample_sample_density_plot(x[,n], main)})
}



#' Generate a nice sample-sample density plot of two samples in matrix x.
#'
#' @param x Matrix of twp columns (samples) to graph.
#' @param fc Fold-change cutoff to consider (2.0 is default).
#' @param main A main title for plot
#'
#' @return Return a ggplot sample-sample density plot.
#' @export
#'
sample_sample_density_plot<-function(x, fc=2.0,main=NULL) {

  # Stats to compute
  cor_r<-stats::cor(x[,1],x[,2])
  diffs<-abs(x[,1]-x[,2])
  diffs_1.5<-length(which(diffs>log2(1.5)))
  diffs_2<-length(which(diffs>log2(2)))

  # Define the x and y mappings as the colnames of the input.
  x_aes<-as.name(colnames(x)[1])
  y_aes<-as.name(colnames(x)[2])
  # ggplot wants a data frame
  df<-as.data.frame(x)
  # Build the scatter plot as a series of hex shapes (but very small, so it looks like points).
  splot<-ggplot2::ggplot(df, ggplot2::aes_(x=x_aes, y=y_aes)) +
    ggplot2::geom_hex(binwidth=c(0.05,0.05)) +   # This specifies hexagons at 0.1x0.1 (x,y) widths.
    ggplot2::scale_fill_continuous(type="viridis",trans="log10") +   # Viridis is blue to yellow range.
    ggplot2::theme_bw() +
    # Set the space between the legend and the plot to a very small amount (save space).
    ggplot2::theme(legend.box.spacing=grid::unit(c(0.01), "npc")) +
    ggplot2::ggtitle(main) +
    # Ensure the plot is square
    ggplot2::coord_fixed(ratio=1) +

    # Add the two-fold lines and the y=x for comparisons.
    ggplot2::geom_abline(slope=1, intercept=0, color="blue", linetype=2,alpha=0.4) +
    ggplot2::geom_abline(slope=1, intercept=1, color="red", linetype=2) +
    ggplot2::geom_abline(slope=1, intercept=-1, color="red", linetype=2)


  # The annotation of correlation, etc. is defined here.
  stats<-sprintf("R=%7.5f\nNum > 1.5-fold: %5d\nNum > 2.0-fold: %5d",cor_r, diffs_1.5, diffs_2)

  # We combine the plot with a caption that includes the stats. Use the font size of the plot.
  splot + ggplot2::labs(caption = stats) #+ theme(plot.caption=element_text(size=splot$theme$text$size+2))
}
