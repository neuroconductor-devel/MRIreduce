#' Process Tissue Segmentation and Compute Volumes
#'
#' This function performs tissue segmentation for specified regions of interest (ROIs),
#' computes intensity and volume metrics, and saves the results for different thresholds.
#' It supports different tissue types and uses parallel processing for efficiency.
#'
#' @param slist A list containing details of the segmentation task including the ROI and module index.
#' @param thresh_vec A numeric vector of threshold values to be used for processing intensity maps.
#' @param tissue_type An integer indicating the type of tissue to process:
#'   - `1` for CFS
#'   - `2` for GM (Gray Matter)
#'   - `3` for WM (White Matter)
#' @param brain_volume_path A string indicating the path to the brain volume RDS file which contains unit voxel information.
#' @param indir A string, the directory path where tissue and intensity data are stored.
#' @param indir.par A string, the directory path where partition data are stored.
#' @param store.dir A string, the directory path where results should be stored.
#'
#' @details The function reads the intensity and tissue data from the specified directories,
#' adjusts the data by removing the first three rows, and computes the mean intensity
#' and volume metrics based on the specified tissue type. Results are saved in RDS format.
#'
#' @return Returns a vector of messages indicating the completion status for each threshold processed.
#' These messages are also printed to the console.
#'
#' @examples
#' \dontrun{
#'   slist <- list(ROI = "hippocampus", Index = 1, Variables = list(c('V1', 'V2', 'V3'),c('V4', 'V5', 'V6')))
#'   thresh_vec <- seq(0.8, 1.0, by = 0.05)
#'   tissue_segment(slist, thresh_vec, 3, "/path/to/brain_volume.rds",
#'                  "/path/to/data", "/path/to/partitions", "/path/to/store")
#' }
#'
#' @export
tissue_segment <- function(slist, thresh_vec, tissue_type, brain_volume_path, indir, indir.par, store.dir) {
  require(dplyr)
  if (tissue_type == 3){
    t.name = 'WM'
  } else{
    if (tissue_type == 2){
      t.name = 'GM'
    }else{
      t.name = 'CFS'
    }
  }

  ROI = slist$ROI #name of interest position

  par.path = file.path(indir.par, ROI) #Path for partition directory per ROI

  if (!dir.exists(file.path(store.dir,ROI))) {
    # Directory does not exist, so create it
    dir.create(file.path(store.dir,ROI), recursive = TRUE)
  }
  output_path = file.path(store.dir,ROI,t.name)
  if (!dir.exists(output_path)) {
    # Directory does not exist, so create it
    dir.create(output_path, recursive = TRUE)
  }

  #Read in tissue data
  tissue.path = file.path(indir,paste0('tissues_', ROI, '.rds'))
  fl.tissue <- readRDS(file = tissue.path)
  fl.tissue <- fl.tissue[-c(1:3),] #Very important line
  indices <- which(fl.tissue == tissue_type, arr.ind = TRUE)

  #Read in intensity data
  intensity.path = file.path(indir, paste0('intensities_', ROI, '.rds'))
  fl.intensity = readRDS(file = intensity.path)
  fl.intensity <- fl.intensity[-c(1:3),] #Very important line

  #Get module number
  module = slist$Index

  # Extract relevant rows from df_intensity based on the index
  row <- unique(indices[,1])
  data_tissue <- fl.tissue[row, ]
  rowname <- rownames(data_tissue) #images
  data <- fl.intensity[rowname, ]

  #Brain unit voxel
  brainv_df <- readRDS(file = brain_volume_path)

  process_threshold <- function(thred) {
    # Define the output file paths
    intensity_file_path <- file.path(output_path, paste0('intensity_', module, '_', thred, '_.rds'))
    volume_file_path <- file.path(output_path, paste0('volume_', module, '_', thred, '_.rds'))

    # Check if both files exist
    if (file.exists(intensity_file_path) && file.exists(volume_file_path)) {
      return(paste("No processing needed; output files already exist for threshold", thred))
    }

    # File processing operations if files do not exist
    file_path <- file.path(par.path, paste0('Par_intensities_map_', module, '_', thred, '_', ROI, '_.rds'))
    par_map <- readRDS(file_path)

    # Initialize matrices to store results
    count_voxel_matrix <- matrix(0, nrow = length(rowname), ncol = 1)
    intensity_mean_matrix <- matrix(NA, nrow = length(rowname), ncol = nrow(par_map))

    results <- lapply(1:nrow(par_map), function(i) {
      mapping <- par_map[i, 'mapping'][[1]]
      condition <- data_tissue[, mapping, drop = FALSE] == tissue_type

      intensity_mean <- rowMeans(data[, mapping, drop = FALSE] * condition, na.rm = TRUE)
      count_voxel <- rowSums(condition)

      return(list(intensity_mean, count_voxel))
    })

    for (i in 1:length(results)) {
      intensity_mean_matrix[, i] <- results[[i]][[1]]
      count_voxel_matrix[, 1] = count_voxel_matrix[, 1] + results[[i]][[2]]
    }

    colnames(intensity_mean_matrix) <- par_map[, 'variable']
    colnames(count_voxel_matrix) <- 'count_voxel'

    intensity_df <- data.frame(intensity_mean_matrix, fname = rowname)
    volume_df <- data.frame(count_voxel_matrix, fname = rowname)

    # Perform the inner join and calculations as before
    joined_df <- inner_join(volume_df, brainv_df, by = "fname")
    joined_df <- joined_df %>%
      mutate(volume = count_voxel * unit_voxel) %>%
      select(fname, volume)

    # Save the output files
    saveRDS(intensity_df, file = intensity_file_path)
    saveRDS(joined_df, file = volume_file_path)

    return(paste("Processed and saved data for threshold", thred))
  }

  rslts <- lapply(thresh_vec, process_threshold)
  # Optionally, print results
  print(rslts)
}
