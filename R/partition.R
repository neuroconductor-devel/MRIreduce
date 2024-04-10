#' Perform Partition Analysis on ROI Intensities
#'
#' This function conducts a partition analysis on a set of regions of interest (ROIs) based on intensity measures.
#' It reads in intensity data for each ROI, filters by specified features, and performs partition analysis at different information thresholds.
#' The results, including reduced data and mapping keys, are saved to specified output directories.
#'
#' @param liste A list containing three elements:
#'   1. The name of the ROI as a string.
#'   2. The super partition element index (integer) within the ROI.
#'   3. A nested list of length 1 containing a vector of feature names (strings) identifying features in the super partition.
#' @param thresh_vec A numeric vector of thresholds to use for partition analysis.
#' @param indir A string specifying the directory path where the intensity files (`intensities_ROI.rds`) are located.
#' @param outdir A string specifying the directory path where the output files should be saved.
#'   The function creates a subdirectory for each ROI within this directory.
#'
#' @return Returns a vector containing the ROI name and the super partition element index.
#'
#' @examples
#' # Assuming appropriate structures and data exist in 'indir' and 'outdir' paths
#' result <- parfun(liste = my_list, thresh_vec = c(0.1, 0.2), indir = "path/to/indir", outdir = "path/to/outdir")
#'
#' @import Partition
#' @export
parfun = function(liste, thresh_vec, indir, outdir){
  roi = liste[[1]]
  mydir = paste0(outdir, roi, "/")
  if(!dir.exists(mydir)) dir.create(mydir)

  module = liste[[2]]
  features = liste[[3]][[module]]

  # read in intensities
  fl = paste0(indir, "intensities_", roi,".rds")
  dati = readRDS(fl) # first 3 rows of dati are location coordinates, intensities are below
  dati = dati[-1:-3,features] # we don't need coordinates for this Partition analysis

  # conduct partition for each information threshold
  for(tind in 1:length(thresh_vec)){
    thresh = thresh_vec[tind]
    tmp = Partition::partition(dati, threshold=thresh) #
    tmpdat = as.data.frame(tmp$reduced_data)
    tmpmap = as.data.frame(mapping_key(tmp))

    fl.dat = paste0(mydir, "Par_intensities_", module, "_", thresh, "_", roi, "_.rds")
    fl.map = paste0(mydir, "Par_intensities_map_", module, "_", thresh, "_", roi, "_.rds")

    saveRDS(tmpdat, file=fl.dat)
    saveRDS(tmpmap, file=fl.map)
  } # end thresh loop
  outp = c(roi, module)
  return(outp)
}
