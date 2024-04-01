#' Perform T1 Image Processing and Analysis
#'
#' This function performs a series of operations on T1-weighted MRI images,
#' including reading the image, reorienting, bias correction, brain extraction,
#' segmentation, registration to the Eve template, and extraction of intensity
#' and tissue data. The brain volume is calculated, and the output includes
#' arrays of intensities, tissue segmentation, and brain volume measurements.
#'
#' @param fpath Character string specifying the file path of the T1-weighted image.
#' @param eve_brain File path to the Eve template brain image for registration.
#' @param eve_brain_mask File path to the Eve template brain mask for segmentation.
#' @param outpath Character string specifying the file path of processed image.
#' @param fsl_path File path to FSL installation on your computer.
#' @param fsl_outputtype Specify the output type of fsl. Default is "NIFTI_GZ".
#' @return A named list with elements "intensities" for the image intensities,
#'         "tissues" for the segmented tissue data, and "brain_volume_cm3" as a vector for
#'         unit of volume (cm) for one voxel and the total intracranial volume in cubic centimeters. Additionally,
#'         an .Rdata file is saved to a specified directory with this information.
#' @examples
#' \dontrun{
#' # Example usage:
#' result <- eve_T1("path/to/T1/image.nii.gz", "path/to/eve_brain.nii.gz", "path/to/eve_brain_mask.nii.gz")
#' }
#' @importFrom neurobase readnii writenii
#' @importFrom fslr fslreorient2std fslbet fsl_biascorrect fslstats
#' @importFrom ANTsR flirt fast
#' @importFrom dplyr %>%
#' @export
eve_T1 <- function(fpath, eve_brain, eve_brain_mask, outpath, fsl_path, fsl_outputtype = "NIFTI_GZ") {
  require(neurobase)
  require(ANTsR)
  require(fslr)
  require(EveTemplate)
  require(MNITemplate)
  require(dplyr)

  options(fsl.path = fsl_path)
  options(fsl.outputtype = fsl_outputtype)

  # get file name
  tmp <- strsplit(fpath, "/", fixed = TRUE)
  fnm <- tmp[[1]][length(tmp[[1]])]

  # read in an image
  print(paste(Sys.time(), "Reading image:", fpath))
  t1 <- neurobase::readnii(fpath)

  # unit of volume (cm) for one voxel
  vres <- oro.nifti::voxres(t1, units = "cm")

  # reorient to standard
  print(paste(Sys.time(), "Reorienting image:", fnm))
  t1_ro <- fslreorient2std(t1)

  vres <- oro.nifti::voxres(t1_ro, units = "cm")

  print(paste(Sys.time(), "Bias correct:", fnm))
  bc_t1 <- fsl_biascorrect(file = t1_ro)

  temp_dir <- tempdir()
  fl <- tempfile(pattern = fnm, tmpdir = temp_dir, fileext = ".nii.gz")
  writenii(nim = bc_t1, filename = fl)

  # FSL’s Brain Extraction Tool (BET)
  # To avoid warnings:
  # export LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH
  print(paste(Sys.time(), "Brain extraction:", fnm))
  bc_bet <- fslbet(
    infile = fl,
    opts = "-B -f 0.1 -v", # from Popescu et al.
    betcmd = "bet",
    intern = FALSE
  )


  # Calculate brain volume
  # segmentation of image into white matter (class = 3), grey matter (class = 2),
  # and cerebrospinal fluid (CFS) (class = 1)
  print(paste(Sys.time(), "Brain volume segmentation:", fnm))
  writenii(nim = bc_bet, filename = fl)
  msk_fast <- fast(fl, type = "T1", retimg = TRUE, opts = "-N", reorient = FALSE)

  # Initialize a variable to store the total intracranial volume
  total_icv <- 0

  # Loop through the tissue types (1: CSF, 2: Grey Matter, 3: White Matter)
  for (tissue_type in 1:3) {
    # For each tissue type, use fslstats to calculate the volume.
    # Specify the lower (-l) and upper (-u) threshold to isolate each tissue type.
    # The '-V' option returns the volume of the specified tissue type.
    tissue_volume_info <- fslstats(file = msk_fast, opts = paste0("-l ", tissue_type - 0.5, " -u ", tissue_type + 0.5, " -V"))

    # The output is a character string that includes the voxel count and the total volume.
    # Here, we split the string and convert the second part (the volume in mm^3) to numeric.
    tissue_volume <- as.numeric(strsplit(tissue_volume_info, " ")[[1]][2])

    # Sum up the volume to calculate total ICV
    total_icv <- total_icv + tissue_volume
  }

  # register images to Eve templet
  print(paste(Sys.time(), "Register to Eve:", fnm))
  bc_bet <- flirt(infile = bc_bet, reffile = eve_brain)

  # extract intensity data from the image
  print(paste(Sys.time(), "Extract intensities:", fnm))
  adat <- oro.nifti::img_data(bc_bet) # array

  # segmentation of image into white matter (class = 3), grey matter (class = 2),
  # and cerebrospinal fluid (CFS) (class = 1)
  print(paste(Sys.time(), "Segmentation:", fnm))
  # write and read to fix file-not-found error
  writenii(nim = bc_bet, filename = fl)
  msk_fast <- fast(fl, type = "T1", retimg = TRUE, opts = "-N", reorient = FALSE) # -N means no inhomogeneity correction

  print(paste(Sys.time(), "Tissue array:", fnm))
  adat_fast <- oro.nifti::img_data(msk_fast) # array

  # output array intensities and tissues registered to Eve
  outp <- vector("list", 3)
  names(outp) <- c("intenisites", "tissues", "brain_volume_cm3")
  outp[[1]] <- adat
  outp[[2]] <- adat_fast
  outp[[3]] <- c(vres, total_icv)

  outfile <- paste0(outpath, fnm, ".Rdata")
  save(outp, file = outfile)
  unlink(fl)
  return(outp)
} # end eve_T1()
