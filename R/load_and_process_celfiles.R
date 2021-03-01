#' Load and process a celfile dataset.
#'
#' It requires a
#' table/tibble of samples with at least one column called
#' "filenames". These filenames will be read in as CEL files.
#' Any associated data in the table will be added as pData to
#' the AffyBatch.
#'
#' This does quite a bit (by calling other functions). First, it
#' loads the data into an AffyBatch and includes metadata in the
#' object (associated with the input parameter). Then, RMA is used
#' on the sample set. RSI is calculated. And the resulting AffyBatch and
#' ExpressionSet are returned as an object from this function.
#'
#' @param data A tibble of samples
#'
#' @return An list of (celfiles, rma_exprs).
#' @export
#'
load_and_process_celfiles<-function(data) {
  process_celfiles(load_celfiles(data))
}


#' Load an experimental CEL file dataset.
#'
#' This is custom for the RSI_CLIA project. It requires a
#' table/tibble of samples with at least one column called
#' "filenames". These filenames will be read in as CEL files,
#' any associated data in the table will be added as pData to
#' the AffyBatch.
#'
#' @param data A tibble of filenames (and perhaps other metadata).
#'
#' @return An AffyBatch object.
#' @export
load_celfiles<-function(data) {

  assertthat::assert_that(!missing(data), msg="Parameter data is missing.")
  assertthat::assert_that(assertthat::has_name(data, "filename"), msg="The field 'filename' is not in input data.")
  assertthat::assert_that(nrow(data)>0, msg="Data is empty.")

  # This is an optional thing, add in the sample name to ensure it exists.
  if (! "sample_name" %in% colnames(data)) data$sample_name<-basename(data$filename)

  message(sprintf("Loading %d CEL files.", nrow(data)))
  celfiles<-affy::ReadAffy(filenames=data$filename, sampleNames=data$sample_name)
  Biobase::pData(celfiles)<-data.frame(data, row.names=affy::sampleNames(celfiles))

  celfiles
}

#' Process CEL files by RMA normalization, then compute RSI.
#'
#' @param celfiles An AffyBatch of CEL files.
#'
#' @return A list of celfiles, rma_exprs corresponding to the input. The pData of both include a new variable, called rsi.
#' @export
#'
process_celfiles<-function(celfiles) {

  assertthat::assert_that(!missing(celfiles), msg="Parameter celfiles is missing.")
  assertthat::assert_that(methods::is(celfiles, "AffyBatch"), msg="Celfiles is not an AffyBatch.")


  message("Performing RMA.")
  rma_exprs<-affy::rma(celfiles)

  rma_exprs$rsi<-rsi::rsi(rma_exprs)
  # Distribute back to the AffyBatch
  celfiles$rsi<-rma_exprs$rsi

  list(celfiles=celfiles, rma_exprs=rma_exprs)
}
