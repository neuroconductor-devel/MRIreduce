#' Process and Extract ROI Intensity and Tissue Type Information
#'
#' This function processes image data for a specified tissue index by loading
#' pre-defined region of interest (ROI) labels and associated information. It
#' extracts regions of interest based on these labels and concatenates intensity
#' and tissue type information across multiple images.
#'
#' @param tind Integer, specifying the tissue index corresponding to the ROI.
#' @param nfl Character vector, specifying the file paths to the RData files containing image data.
#' @param outpath_intensity Character, specifying the output path for intensity data.
#' @param outpath_tissue Character, specifying the output path for tissue data.
#' @param outpath_volume Character, specifying the output path for volume data.
#'
#' @return Character, the text label of the processed ROI.
#'
#' @examples
#' \dontrun{
#'   nfl <- list.files(pattern = "\\.RData$")
#'   iproc(tind = 1, nfl = nfl,
#'         outpath_intensity = "path/to/intensity",
#'         outpath_tissue = "path/to/tissue",
#'         outpath_volume = "path/to/volume")
#' }
#'
#' @export
iproc = function(tind, nfl, outpath_intensity, outpath_tissue, outpath_volume){
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
  if(tind == 1){
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
    if(tind == 1){
      vvec[find, "brainVolume"] = vdat
      vvec[find, "fname"] = tmp
    }
  } # end find loop

  # output intensity and tissue files for ROI
  inm = paste0(outpath_intensity, "/intensities_", ROI, tind, ".rds")
  saveRDS(loci, file=inm)
  tnm = paste0(outpath_tissue, "/tissues_", ROI, tind, ".rds")
  saveRDS(loct, file=tnm)
  if(tind == 1){
    vnm = paste0(outpath_volume, "/Bvolumes.rds")
    saveRDS(vvec, file=vnm)
  }
  print(paste("End", ROI, Sys.time()))
  return(ROI)
}
