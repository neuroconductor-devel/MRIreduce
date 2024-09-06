#' Process T1-weighted Brain MRI Data with FSL and Register to EVE Atlas
#'
#' This function takes a T1-weighted brain MRI image, performs bias correction,
#' reorientation to standard space, brain extraction using FSL's BET, and registration
#' to the EVE brain template. It then segments the brain volume into different tissue
#' types, calculates intracranial volume and outputs the results as an Rdata file.
#'
#' @param fpath Character string specifying the path to one T1-weighted MRI file. The file should be in NIFTI file format (.nii.gz).
#' @param outpath Character string specifying the output directory where the results
#' will be saved as an Rdata file.
#' @param fsl_path Character string specifying the directory where FSL is installed.
#' @param fsl_outputtype Character string specifying the FSL output file type.
#' Defaults to "NIFTI_GZ".
#'
#' @return Saves the image data after processing to the specified output path.
#' Outputs an Rdata file containing three components: image intensities array,
#' segmented tissues array, and brain volume metrics.
#'
#' @details
#' The function begins by loading the EVE brain template for image registration. It then reads the
#' T1-weighted MRI file and reorients it to standard space using FSL's `fslreorient2std`.
#' Following reorientation, it applies bias correction with FSL's `fsl_biascorrect`, which
#' is necessary to correct for field inhomogeneities that can affect quantitative analysis.
#' The next step involves using FSL’s Brain Extraction Tool (BET) to isolate the brain from
#' non-brain tissue which is crucial for accurate subsequent analysis.
#' After brain extraction, the image is registered to the EVE brain atlas using FLIRT,
#' ensuring that the image is aligned with a standard coordinate space for comparable
#' anatomical analysis. Subsequent to registration, the function uses FSL's FAST tool to
#' segment the brain into white matter, grey matter, and cerebrospinal fluid, which
#' are essential for studying brain structure and pathology. Finally, it calculates the
#' volume of these tissues, providing key data points for clinical and research applications.
#' Each step logs its progress with timestamps, aiding in debugging and optimization of processing times.
#'
#' @examples
#' eve_T1("path/to/your/image.nii.gz", "path/to/output", "/usr/local/fsl",  "NIFTI_GZ")
#'
#' @importFrom neurobase readnii writenii
#' @importFrom fslr fslreorient2std fsl_biascorrect fslbet fast fslstats
#' @importFrom dplyr %>%
#' @importFrom EveTemplate getEvePath
#' @export
eve_T1 <- function(fpath, outpath, fsl_path = '/Users/jinyaotian/fsl', fsl_outputtype = "NIFTI_GZ") {
  #Brain template
  eve_brain_fname = getEvePath("Brain")
  eve_brain = readnii(eve_brain_fname)

  options(fsl.path = fsl_path)
  options(fsl.outputtype = fsl_outputtype)

  # get file name
  tmp <- strsplit(fpath, "/", fixed = TRUE)
  fnm <- tmp[[1]][length(tmp[[1]])]

  # read in an image
  print(paste(Sys.time(), "Reading image:", fpath))
  t1 <- neurobase::readnii(fpath)

  # reorient to standard
  print(paste(Sys.time(), "Reorienting image:", fnm))
  t1_ro <- fslreorient2std(t1)

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

  # register images to Eve templet
  print(paste(Sys.time(), "Register to Eve:", fnm))
  bc_bet <- flirt(infile = bc_bet, reffile = eve_brain)

  # extract intensity data from the image
  print(paste(Sys.time(), "Extract intensities:", fnm))
  adat <- oro.nifti::img_data(bc_bet) # array

  # segmentation of image into white matter (class = 3), grey matter (class = 2),
  # and cerebrospinal fluid (CFS) (class = 1)
  print(paste(Sys.time(), "Brain volume segmentation:", fnm))
  writenii(nim = bc_bet, filename = fl)
  msk_fast <- fast(fl, type = "T1", retimg = TRUE, opts = "-N", reorient = FALSE)

  #Calculate Brain Volume
  print(paste(Sys.time(), "Intracranial volume calculation:", fnm))
  vres <- oro.nifti::voxres(bc_bet, units = "cm")
  # Initialize a variable to store the total intracranial volume
  total_icv <- 0
  # Initialize a variable to store the brain volume = WM + GM
  bv <- 0
  # Loop through the tissue types (1: CSF, 2: Grey Matter, 3: White Matter)
  for (tissue_type in 1:3) {
    # For each tissue type, use fslstats to calculate the volume.
    # Specify the lower (-l) and upper (-u) threshold to isolate each tissue type.
    # The '-V' option returns the volume of the specified tissue type.
    tissue_volume_info <- fslstats(file = msk_fast, opts = paste0("-l ", tissue_type - 0.5, " -u ", tissue_type + 0.5, " -V"))

    # The output is a character string that includes the voxel count and the total volume.
    # Here, we split the string and convert the second part (the volume in cm^3) to numeric.
    tissue_volume <- as.numeric(strsplit(tissue_volume_info, " ")[[1]][2])

    # Sum up the volume to calculate total ICV
    total_icv <- total_icv + tissue_volume * vres
    if(tissue_type != 1){
      bv <- bv + tissue_volume * vres
    }
  }
  print(paste(Sys.time(), "Tissue array:", fnm))
  adat_fast <- oro.nifti::img_data(msk_fast) # array

  # output array intensities and tissues registered to Eve
  outp <- vector("list", 3)
  names(outp) <- c("intenisites", "tissues", "brain_volume_cm3")
  outp[[1]] <- adat
  outp[[2]] <- adat_fast
  outp[[3]] <- c(vres,bv, total_icv)

  outfile <- file.path(outpath,paste0(fnm, ".Rdata"))
  save(outp, file = outfile)
  unlink(fl)
} # end eve_T1()

################################################################################
# #test
# fpath = "/Users/jinyaotian/Downloads/whims_test/raw_image"
# opath = "/Users/jinyaotian/Downloads/whims_test/eve_T1"
# files = list.files(fpath, full.names = TRUE)
# # Apply the eve_T1 function to each file with additional parameters
# results <- lapply(files, function(x) eve_T1(fpath = x, outpath = opath))

