#' Process Intensity and Tissue Data for Specified ROI
#'
#' This function loads ROI labels and descriptions, processes a set of image data files
#' to extract intensity and tissue information for a specified Region of Interest (ROI),
#' and optionally aggregates brain volume data. Outputs are saved as RDS files.
#'
#' @param tind An integer indicating the target ROI based on the integer labels in `eve_label_array.RData`.
#' @param nfl A character vector of paths to RDS files, each containing a list with "intensities", "tissues", and "brain_volume_cm3" for different images.
#' @param outpath_intensity The directory path where intensity data RDS files will be saved.
#' @param outpath_tissue The directory path where tissue data RDS files will be saved.
#' @param outpath_volume Optional. The directory path where brain volume data RDS files will be saved, applicable only when `tind` is 1. If `NULL`, volume data is not processed.
#'
#' @details The function starts by loading the `eve_label_array.RData` and `eve_label_info_dataframe.RData`
#' files to access ROI labels and descriptions. For each file in `nfl`, it extracts and processes intensity
#' and tissue data corresponding to the specified `tind`. If `outpath_volume` is provided and `tind` equals 1,
#' brain volume data is also processed and saved. Finally, the function saves the processed data as RDS files
#' in the specified output paths.
#'
#' @return The text label of the processed ROI.
#' @export
#'
#' @examples
#' \dontrun{
#'   iproc(tind = 1,
#'         nfl = c("/path/to/image1.rds", "/path/to/image2.rds"),
#'         outpath_intensity = "/path/to/output/intensity",
#'         outpath_tissue = "/path/to/output/tissue",
#'         outpath_volume = "/path/to/output/volume")
#' }
iproc <- function(tind, nfl, outpath_intensity, outpath_tissue, outpath_volume = NULL) {
  # Load Eve labels
  fl = system.file("extdata", "eve_label_array.RData", package = "WHIMs")
  load(file=fl) # dat_eve, array of integer labels, contains integer labels for each voxel in the image data.

  # Load text_label and structure descriptions of ROI's
  fl = system.file("extdata", "eve_label_info_dataframe.RData", package = "WHIMs")
  load(file=fl) # lab_df, DataFrame, contains mappings from integer labels to text labels and structures.

  ROI = lab_df[as.character(tind), "text_label"]
  print(paste("Start", ROI, Sys.time()))
  loc = which(dat_eve == tind,arr.ind = TRUE)
  loci = as.data.frame(t(which(dat_eve == tind,arr.ind = TRUE)))
  names(loci) = paste0("V", 1:ncol(loci))
  loct = loci
  if(!is.null(outpath_volume) && tind == 1){
    vvec = as.data.frame(matrix(NA, nrow=0, ncol=2))
    names(vvec) = c("fname", "brainVolume")
  }
  t1 = Sys.time()

  for(find in 1:length(nfl)){
    print(c(find, Sys.time(), nfl[find]))
    fpath = nfl[find]
    load(fpath) # outp, list with "intenisites", "tissues", "brain_volume_cm3"
    idat = outp[[1]]
    tdat = outp[[2]]
    if(tind == 1) vdat = outp[[3]]
    rm(outp)
    dind = find + 3
    loci[dind, ] = idat[loc]
    tmp = strsplit(nfl[find], "/", fixed = TRUE)[[1]]
    tmp = tmp[length(tmp)]
    tmp = strsplit(tmp, ".", fixed = TRUE)[[1]][1]
    rownames(loci)[dind] = tmp
    loct[dind, ] = tdat[loc]
    rownames(loct)[dind] = tmp
    rm(idat, tdat, tmp)
    if(!is.null(outpath_volume) && tind == 1){
      vvec[find, "brainVolume"] = vdat
      vvec[find, "fname"] = tmp
    }
  } # end find loop

  # output intensity and tissue files for ROI
  inm = paste0(outpath_intensity, "/intensities_", ROI, tind, ".rds")
  saveRDS(loci, file=inm)
  tnm = paste0(outpath_tissue, "/tissues_", ROI, tind, ".rds")
  saveRDS(loct, file=tnm)
  if(!is.null(outpath_volume) && tind == 1){
    vnm = paste0(outpath_volume, "/Bvolumes.rds")
    saveRDS(vvec, file=vnm)
  }
  print(paste("End", ROI, Sys.time()))
  return(ROI)
}
