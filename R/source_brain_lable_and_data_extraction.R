#Program: brain lable and data extraction
# Programmer: Jinyao Tian
# Date: 03/20/2024
library(neurobase)
library(ANTsR)
library(fslr)
library(EveTemplate)
library(MNITemplate)
library(dplyr)

################################################################################
eve_T1 = function(fpath, eve_brain, eve_brain_mask){
  require(neurobase)
  require(ANTsR)
  require(fslr)
  require(EveTemplate)
  require(MNITemplate)
  require(dplyr)

  options(fsl.path = "/Users/jinyaotian/fsl")
  options(fsl.outputtype = "NIFTI_GZ")

  # get file name
  tmp = strsplit(fpath, "/", fixed=TRUE)
  fnm = tmp[[1]][length(tmp[[1]])]

  # read in an image
  print(paste(Sys.time(), "Reading image:", fpath))
  t1 = neurobase::readnii(fpath)

  # unit of volume (cm) for one voxel
  vres = oro.nifti::voxres(t1, units = "cm")

  # reorient to standard
  print(paste(Sys.time(), "Reorienting image:", fnm))
  t1_ro = fslreorient2std(t1)

  vres = oro.nifti::voxres(t1_ro, units = "cm")

  print(paste(Sys.time(), "Bias correct:", fnm))
  bc_t1 = fsl_biascorrect(file = t1_ro)

  fl = paste0("/Users/jinyaotian/Downloads/baseline_tmp/", fnm)
  writenii(nim=bc_t1, filename=fl)

  # FSL’s Brain Extraction Tool (BET)
  # To avoid warnings:
  # export LD_LIBRARY_PATH=/usr/lib/fsl/5.0:$LD_LIBRARY_PATH
  print(paste(Sys.time(), "Brain extraction:", fnm))
  bc_bet = fslbet(infile = fl,
                  opts = "-B -f 0.1 -v",  # from Popescu et al.
                  betcmd = "bet",
                  intern = FALSE)


  # Calculate brain volume
  # segmentation of image into white matter (class = 3), grey matter (class = 2),
  # and cerebrospinal fluid (CFS) (class = 1)
  print(paste(Sys.time(), "Brain volume segmentation:", fnm))
  fl = paste0("/Users/jinyaotian/Downloads/baseline_tmp/", fnm)
  writenii(nim=bc_bet, filename=fl)
  msk_fast = fast(fl, type="T1", retimg=TRUE, opts='-N', reorient=FALSE)

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
  bc_bet = flirt(infile=bc_bet, reffile=eve_brain)

  # label image by eve parcellation
  # print(paste(Sys.time(), "Label to Eve:", fnm))
  # bc_bet = mask_img(bc_bet, mask = eve_brain_mask)

  # extract intensity data from the image
  print(paste(Sys.time(), "Extract intensities:", fnm))
  adat = oro.nifti::img_data(bc_bet) # array

  # segmentation of image into white matter (class = 3), grey matter (class = 2),
  # and cerebrospinal fluid (CFS) (class = 1)
  print(paste(Sys.time(), "Segmentation:", fnm))
  # write and read to fix file-not-found error
  fl = paste0("/Users/jinyaotian/Downloads/baseline/", fnm)
  writenii(nim=bc_bet, filename=fl)
  msk_fast = fast(fl, type="T1", retimg=TRUE, opts='-N', reorient=FALSE) # -N means no inhomogeneity correction

  print(paste(Sys.time(), "Tissue array:", fnm))
  adat_fast = oro.nifti::img_data(msk_fast) # array

  # output array intensities and tissues registered to Eve
  outp = vector('list', 3)
  names(outp) = c("intenisites", "tissues", "brain_volume_cm3")
  outp[[1]] = adat
  outp[[2]] = adat_fast
  outp[[3]] = c(vres,total_icv)

  outfile = paste0("/Users/jinyaotian/Downloads/baseline_Rdata/", fnm, ".Rdata")
  save(outp, file=outfile)
  #return(outfile)
} # end eve_T1()

################################################################################
#Parallel running
if(1){
  #Data Preprocess

  #Exclude images
  id_excluded <- read.csv(file = '/Users/jinyaotian/Desktop/Bio/Project_783/Data/ID_excluded.csv')%>%
    select(Subject)

  #T1 and Fl files from baseline folder
  baseline_path <- "/Users/jinyaotian/Downloads/baseline"
  # List T1 files
  t1_files <- list.files(path = baseline_path, pattern = "T1", full.names = TRUE) #1319 T1 files in baseline
  # List FL files
  fl_files <- list.files(path = baseline_path, pattern = "FL", full.names = TRUE) #1262 FL files in baseline


  filter_files <- function(files) {
    # Filter and return file names not in the id_excluded
    filtered_files <- Filter(function(file) {
      # Extract the numeric ID from the file name
      id <- sub("-.*", "", basename(file))
      # Check if the ID is not in the id_excluded vector
      !id %in% id_excluded$Subject
    }, files)

    return(filtered_files)
  }

  # Apply the function to T1 and FL files
  clean_t1_files <- filter_files(t1_files) #1279
  clean_fl_files <- filter_files(fl_files) #1222
  #40 images are removed from the baseline


  library(foreach)
  library(doParallel)
  # Start time measurement
  start_time <- Sys.time()

  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  # Parallel processing of clean_t1_files using foreach
  results <- foreach(t1_file = clean_t1_files[63:length(clean_t1_files)], .packages = c("dplyr", "ANTsR", "fslr","neurobase", "EveTemplate", "MNITemplate")) %dopar% {
    # Assuming eve_T1 is your function that processes each file
    # You might need to make sure eve_brain and eve_brain_mask are available for each parallel process
    eve_brain_fname = getEvePath("Brain")
    eve_brain = readnii(eve_brain_fname) # read in brain-extracted Eve T1 image
    eve_brain_mask = readEve(what = "Brain_Mask")
    # This can involve loading or passing them within the foreach loop if necessary
    eve_T1(t1_file, eve_brain, eve_brain_mask)
  }

  # End time measurement
  end_time <- Sys.time()

  # Print execution time
  print(end_time - start_time)

  # Deregister the parallel backend and stop the cluster
  stopCluster(cl)
  registerDoSEQ()
}



