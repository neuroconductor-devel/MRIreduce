#' Process FLAIR Neuroimages and Register to EVE Brain Template
#'
#' This function processes FLAIR (Fluid-Attenuated Inversion Recovery) neuroimages for better lesion detection,
#' especially in the periventricular area by suppressing the CSF signal. It involves steps such as reading the image,
#' reorienting, bias correction using N4, brain extraction, registration to the EVE template, and tissue segmentation.
#' It finally calculates the intracranial volume and outputs the data.
#'
#' @param fpath Character string specifying the path to the FLAIR image file. The file should be in NIFTI file format (.nii.gz).
#' @param outpath Character string specifying the output directory where the processed data is saved.
#' @param fsl_path Character string specifying the path to the FSL software on the system.
#' @param fsl_outputtype Character string specifying the type of output file format for FSL; defaults to "NIFTI_GZ".
#'
#' @return Returns a list containing three elements: `intensities`, `tissues`, and `brain_volume_cm3`.
#' Each element corresponds to the array of intensities, the segmented tissue data, and the calculated brain volumes,
#' respectively. The function also saves these results as an .Rdata file at the specified output path.
#'
#' @details The function uses specific FSL tools for image processing steps such as reorientation to standard space,
#' bias correction with N4 method from ANTsR, and brain extraction using a robust method from `extrantsr`.
#' Segmentation into different tissue types (CSF, grey matter, and white matter) is performed using FSL's FAST tool.
#' Volumes are calculated based on the segmented tissues.
#'
#' @examples
#' \dontrun{
#'    eve_Fl("path/to/your/flair/image.nii.gz",
#'          "path/to/output/",
#'          "/usr/local/fsl",
#'          "NIFTI_GZ")
#' }
#'
#' @importFrom neurobase readnii writenii
#' @importFrom extrantsr bias_correct fslbet_robust
#' @importFrom fslr fslreorient2std fslstats
#' @importFrom oro.nifti img_data voxres
#' @export

eve_Fl <- function(fpath,outpath, fsl_path, fsl_outputtype = "NIFTI_GZ") {
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
  fl <- neurobase::readnii(fpath)

  # reorient to standard
  print(paste(Sys.time(), "Reorienting image:", fnm))
  #This step is the exact same as for T1 format
  fl_ro <- fslreorient2std(fl)

  # Here we will use the bias_correct function in extrantsr,
  # which calls n4BiasFieldCorrection from ANTsR.
  print(paste(Sys.time(), "Bias correct:", fnm))
  bc_fl = bias_correct(fl_ro, correction = "N4")

  temp_dir <- tempdir()
  fl_name <- tempfile(pattern = fnm, tmpdir = temp_dir, fileext = ".nii.gz")
  writenii(nim = bc_fl, filename = fl_name)

  # FSL’s Brain Extraction Tool (BET)
  #Differ from eve_T1
  print(paste(Sys.time(), "Brain extraction:", fnm))
  bc_bet <- extrantsr::fslbet_robust(infile = fl_name,
                                     correct = FALSE,
                                     reorient = FALSE,
                           remover = "double_remove_neck")

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
  writenii(nim = bc_bet, filename = fl_name)
  msk_fast <- fast(fl_name, retimg = TRUE, opts = "-N", reorient = FALSE) # -N means no inhomogeneity correction

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
    # Here, we split the string and convert the second part (the volume in mm^3) to numeric.
    tissue_volume <- as.numeric(strsplit(tissue_volume_info, " ")[[1]][2])

    # Sum up the volume to calculate total ICV
    total_icv <- total_icv + tissue_volume* vres
    if(tissue_type != 1){
      bv <- bv + tissue_volume* vres
    }
  }
  print(paste(Sys.time(), "Tissue array:", fnm))
  adat_fast <- oro.nifti::img_data(msk_fast) # array

  # output array intensities and tissues registered to Eve
  outp <- vector("list", 3)
  names(outp) <- c("intensities", "tissues", "brain_volume_cm3")
  outp[[1]] <- adat
  outp[[2]] <- adat_fast
  outp[[3]] <- c(vres, bv, total_icv)

  outfile <- paste0(outpath, fnm, ".Rdata")
  save(outp, file = outfile)
  unlink(fl_name)
  return(outp)
} # end eve_Fl()
