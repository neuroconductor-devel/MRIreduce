#Processing FLAIR (Fluid-Attenuated Inversion Recovery) neuroimages involves a pipeline somewhat
#similar to that used for T1-weighted images, with a few adjustments tailored to the specific characteristics of FLAIR images.
#FLAIR images provide high contrast for detecting lesions, especially in the periventricular area, by suppressing the cerebrospinal fluid (CSF) signal.
#This property makes FLAIR particularly useful in the identification and study of white matter lesions and other pathologies.

#Programmer: Jinyao Tian
#Date: 04/08/2024

eve_Fl <- function(fpath, eve_brain, eve_brain_mask, outpath, fsl_path, fsl_outputtype = "NIFTI_GZ") {
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

  print(paste(Sys.time(), "Tissue array:", fnm))
  adat_fast <- oro.nifti::img_data(msk_fast) # array

  # output array intensities and tissues registered to Eve
  outp <- vector("list", 3)
  names(outp) <- c("intenisites", "tissues")
  outp[[1]] <- adat
  outp[[2]] <- adat_fast

  outfile <- paste0(outpath, fnm, ".Rdata")
  save(outp, file = outfile)
  unlink(fl_name)
  return(outp)
} # end eve_Fl()
