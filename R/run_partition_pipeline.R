#' Run Partition Pipeline on Neuroimaging Data
#'
#' This function initializes and executes a partitioning pipeline designed for
#' processing neuroimaging data. It handles tasks such as image processing,
#' super partition analysis, mapping, and combining of data based on specified
#' thresholds and parameters.
#'
#' @param tind Index or identifier for the type of tissue under analysis.
#' @param nfl List of file names (full paths) that need to be processed.
#' @param main_dir Main directory where outputs and intermediate results will be saved.
#' @param tissue_type Type of tissue for segmentation.
#' @param ICC_thresh_vec A vector of threshold values for Intraclass Correlation Coefficient used in the Partition Algorithm.
#' @param num_cores Number of cores to use for parallel processing. Default to 1.
#' @param suppar_thresh_vec Optional; a sequence of threshold values used in Super Partitioning.
#'        Default is a sequence from 0.7 to 1 by 0.01.
#' @param B Optional; the maximum size of modules to be considered in partitioning. Default is 2000.
#' @param outp_volume Optional; a logical indicating whether volume outputs should be generated. Default is TRUE.
#'
#' @details
#' The function configures and runs a series of operations that are typical in
#' neuroimage analysis, especially focusing on ROI-based transformations. Each step
#' of the pipeline, from initial image processing (`iproc`) to the final combination
#' of independent variables with reduced variables (`Cmb_indep_with_dep`), is executed
#' in sequence. Adjustments to the pipeline's behavior can be made by changing the
#' function parameters, allowing for custom analysis flows.
#'
#' @return The function does not return a value but will output results directly to
#'         the specified `main_dir` as side effects of the processing steps.
#'
#' @export
run_partition_pipeline <- function(tind, nfl, main_dir, tissue_type,ICC_thresh_vec, num_cores = 1, suppar_thresh_vec = seq(0.7, 1, 0.01),B = 2000, outp_volume = TRUE) {
  pipeline <- PartitionPipeline$new(
    tind = tind,
    nfl = nfl,
    main_dir = main_dir,
    tissue_type = tissue_type,
    ICC_thresh_vec = ICC_thresh_vec,
    num_cores = num_cores,
    suppar_thresh_vec = suppar_thresh_vec,
    B = B,
    outp_volume = outp_volume
  )

  pipeline$iproc()
  pipeline$supparfun()
  pipeline$map_suppar_roi()
  pipeline$partition_intensities()
  pipeline$tissue_segment()
  pipeline$Cmb_tissue_type()
  pipeline$process_indep_variables()
  pipeline$Cmb_indep_with_dep()
}
