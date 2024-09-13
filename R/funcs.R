# Process Intensity and Tissue Data for Specified ROI
#
# This function loads ROI labels and descriptions, processes a set of image data files
# to extract intensity and tissue information for a specified Region of Interest (ROI),
# and optionally aggregates brain volume data. Outputs are saved as RDS files.
#
# Output:
#   intensities file structure: main_dir/intensities/"intensities_ROI_tind.rds"
#   tissue file structure: main_dir/tissues/"tissues_ROI_tind.rds"
load_required_data <- function(filename, package_name, envir = .GlobalEnv) {
  filepath <- system.file("extdata", filename, package = package_name)
  if (file.exists(filepath)) {
    load(filepath, envir = envir)
  } else {
    stop("File does not exist:", filepath)
  }
}
iproc <- function(tind, nfl, main_dir, outp_volume = TRUE) {
  # Load Eve labels
  # Load necessary data
  load_required_data("eve_label_array.RData", "MRIreduce")
  load_required_data("eve_label_info_dataframe.RData", "MRIreduce")
  # Initial checks and setup
  if (!dir.exists(main_dir)) {
    stop("main_dir does not exist: ", main_dir)
  }

  roi = lab_df[as.character(tind), "text_label"]

  print(paste("Start", roi, Sys.time()))

  loc = which(dat_eve == tind,arr.ind = TRUE)
  loci = as.data.frame(t(which(dat_eve == tind,arr.ind = TRUE)))
  names(loci) = paste0("V", 1:ncol(loci))
  loct = loci

  outpath_volume = file.path(main_dir,"volumes")
  if (dir.exists(outpath_volume)){
    outp_volume = FALSE
  }else{
    dir.create(outpath_volume, recursive = TRUE)
  }
  # Pre-allocate memory for performance and to avoid growing objects inside the loop
  num_files <- length(nfl)
  if (outp_volume) {
    vvec <- data.frame(fname=character(num_files), unit_voxel=numeric(num_files),
                       brainVolume=numeric(num_files), ICV=numeric(num_files), stringsAsFactors=FALSE)
  }
  # Pre-compute file name processing to reduce redundancy
  file_names <- sapply(nfl, function(path) {
    basename(path)
  })
  # Processing each file
  for (find in seq_len(num_files)) {
    print(c(find, Sys.time(), nfl[find]))

    # Load and process the data
    load(nfl[find])  # assuming outp contains the list with "intensities", "tissues", "brain_volume_cm3"

    # Assign data
    loci[find + 3, ] <- outp[[1]][loc]  # Using pre-determined locations
    loct[find + 3, ] <- outp[[2]][loc]

    # Set row names based on pre-computed file names
    rownames(loci)[find + 3] <- file_names[find]
    rownames(loct)[find + 3] <- file_names[find]

    if (outp_volume) {
      vvec[find, ] <- c(file_names[find], outp[[3]])
    }

    # Free up memory
    rm(outp)
  }

  outpath_intensity = file.path(main_dir, "intensities")
  dir.create(outpath_intensity, recursive = TRUE)
  outpath_tissue = file.path(main_dir, "tissues")
  dir.create(outpath_tissue, recursive = TRUE)

  # output intensity and tissue files for ROI
  inm = file.path(outpath_intensity,paste0("intensities_", roi,'_',tind, ".rds"))
  saveRDS(loci, file=inm)
  tnm = file.path(outpath_tissue,paste0("tissues_", roi, '_',tind, ".rds"))
  saveRDS(loct, file=tnm)
  if(outp_volume){
    vnm = file.path(outpath_volume, "Bvolumes.rds")
    saveRDS(vvec, file=vnm)
  }
  message(roi, " : ", "data extraction has been completed")
  return(roi)
}

# Apply suppar() to intensities
# thresh_vec: Numeric vector specifying correlation thresholds.
# B: Integer, the maximum group size.
# output file structure: main_dir/suppar/roi.rds
supparfun <- function(tind, roi,thresh_vec = seq(0.7, 1, 0.01), B = 2000, main_dir){
  inm = file.path(main_dir, "intensities",paste0("intensities_", roi,'_',tind, ".rds"))
  dati = readRDS(inm)
  # dati = dati[4:nrow(dati), ] # remove coordinates
  k = ceiling(ncol(dati) / B)
  dist.thresh = 5 # we are only interested in combining features no more than 5 mm apart
  dir.tmp = paste0(getwd(), "/", "temp_corr_", roi)
  mygrps = suppar(tmp=dati, thresh=thresh_vec, n.chunkf=10000, B=B, dist.thresh=dist.thresh, dir.tmp=dir.tmp)

  pdir = file.path(main_dir, "suppar")
  if(!dir.exists(pdir)){dir.create(pdir, recursive = TRUE)}
  fl = file.path(pdir,paste0(roi, ".rds"))
  saveRDS(mygrps, file=fl)
  message(roi," ", "Super partition has been completed" )
}

# The function generates an R list object that is a map of voxels to superPartitions.
# output:
#   dependent list: main_dir/dep_list/roi.rds
#   independent list: main_dir/indep_list/roi.rds
map_suppar_roi <- function(roi, main_dir) {
  # Function to read RDS files and extract information
  read_rds_file <- function(filename) {
    list_data <- readRDS(filename)
    # Removing prefix and suffix from the filename to extract the ROI name
    roi_name <- gsub("^(.*)\\.rds$", "\\1", basename(filename))
    # Return both components and ROI name
    list(variables_group = list_data[[1]], ind_names_vector = list_data[[2]], roi_name = roi_name)
  }


  # Initialize lists
  dependent_list <- list()
  independent_list <- list()

  file_path = file.path(main_dir, "suppar", paste0(roi, ".rds"))
  data <- read_rds_file(file_path)
  ROI = data$roi_name
  # Building the independent list
  independent_list[[length(independent_list) + 1]] <- list(ROI = data$roi_name, Variables = data$ind_names_vector)

  # Extract and build dependent list
  for (i in seq_along(data$variables_group)) {
    dependent_list[[length(dependent_list) + 1]] <- list(ROI = data$roi_name, Index = i, Variables = data$variables_group)
  }

  if(!dir.exists(file.path(main_dir, "dep_list"))){dir.create(file.path(main_dir, "dep_list"), recursive = TRUE)}
  if(!dir.exists(file.path(main_dir, "indep_list"))){dir.create(file.path(main_dir, "indep_list"), recursive = TRUE)}
  output_dep_path = file.path(main_dir, "dep_list", paste0(ROI,'.rds'))
  output_indep_path = file.path(main_dir, "indep_list", paste0(ROI,'.rds'))
  # Save the extracted lists to RDS files
  saveRDS(independent_list, file = output_indep_path)
  message(data$roi_name, " : ", "independent features list has been created.")
  saveRDS(dependent_list, file = output_dep_path)
  message(data$roi_name, " : ", "reduced features list has been created")
}


#Read ROI summary chunks with T1 data and apply partition to each one.
## Read in super partitions
# slist: super partition element, length(slist) = number of super partition elements
# 1nd level, first element: name of ROI
# 2nd level, second element: super partition element index, 1 to number of super partitions in the ROI
# 3nd level, third element: list of length 1 with vector of feature names identifying features in the super partition

# For each ROI, there will be an "independent_list" for partition algorithm
# liste: sublist of independent_list
# output:
#   intensities file structure: /main_dir/partition/roi/thresh/intensities/par_intensities_module
#   intensities file map structure: /main_dir/partition/roi/thresh/map/par_intensities_module
parfun = function(liste, tind, thresh_vec, main_dir){
  roi = liste$ROI
  mydir = file.path(main_dir,"partition" ,roi)
  if(!dir.exists(mydir)) dir.create(mydir, recursive = TRUE)

  module = liste$Index
  features = liste$Variables[[module]]

  indir = file.path(main_dir, "intensities")
  # read in intensities
  fl = file.path(indir,paste0("intensities_", roi,'_',tind, ".rds"))
  dati = readRDS(fl) # first 3 rows of dati are location coordinates, intensities are below
  dati = dati[-1:-3,features, drop = FALSE] # we don't need coordinates for this Partition analysis
  row_index = row.names(dati) #image name

  # Assuming thresh_vec and other variables are appropriately defined
  results <- lapply(thresh_vec, function(thresh) {
    if(dim(dati)[2] > 1){
      tmp = partition(dati, threshold=thresh)
      tmpdat = as.data.frame(tmp$reduced_data)
      tmpmap = as.data.frame(mapping_key(tmp))
    } else {
      tmpdat = dati
      tmpmap = data.frame(
        variable = features,
        mapping = features,
        information = 1.0,
        indices = "1:1"
      )
    }

    row.names(tmpdat) = row_index
    tmp_path_intensity = file.path(mydir,thresh, "intensities")
    if(!dir.exists(tmp_path_intensity)) dir.create(tmp_path_intensity, recursive = TRUE)
    #directory structure: /roi/thresh/intensities/par_intensities_module
    fl.dat = file.path(tmp_path_intensity, paste0("Par_intensities_", module, ".rds"))

    tmp_path_map = file.path(mydir, thresh, "map")
    if(!dir.exists(tmp_path_map)) dir.create(tmp_path_map, recursive = TRUE)
    #directory structure: /roi/thresh/map/par_intensities_module
    fl.map = file.path(tmp_path_map, paste0("Par_intensities_map_", module,".rds"))

    saveRDS(tmpdat, file=fl.dat)
    saveRDS(tmpmap, file=fl.map)
    message(roi, " module: ", module, " has been partitioned at information threshold: ", thresh)
  })
}

# Tissue segment
# output:
#   Intensities file structure: /main_dir/partition/roi/thresh/tissue/intensities_module.rds
#   Volume file structure: /main_dir/partition/roi/thresh/tissue/volume_module.rds
tissue_segment <- function(liste, thresh_vec, tind, tissue_type, main_dir) {
  if (tissue_type == 3){
    t.name = 'WM'
  } else{
    if (tissue_type == 2){
      t.name = 'GM'
    }else{
      t.name = 'CFS'
    }
  }

  ROI = liste$ROI #name of interest position

  par.path = file.path(main_dir, "partition", ROI)
  output_path = file.path(main_dir, "partition", ROI)

  #Read in tissue data
  tissue.path = file.path(main_dir, "tissues", paste0("tissues_", ROI, '_',tind, ".rds"))
  fl.tissue <- readRDS(file = tissue.path)
  data_tissue <- fl.tissue[-c(1:3),] #Very important line

  #Read in intensity data
  intensity.path = file.path(main_dir, "intensities", paste0("intensities_", ROI,'_',tind, ".rds"))
  fl.intensity = readRDS(file = intensity.path)
  data <- fl.intensity[-c(1:3),] #Very important line
  img_name <- rownames(data)

  #Get module number
  module = liste$Index

  process_threshold <- function(thred) {
    # Define the output file paths
    # Intensities file structure: /main_dir/partition/roi/thresh/tissue/intensities_module.rds
    # Volume file structure: /main_dir/partition/roi/thresh/tissue/volume_module.rds
    if(!dir.exists(file.path(output_path,thred,t.name))){dir.create(file.path(output_path,thred,t.name), recursive = TRUE)}

    intensity_file_path <- file.path(output_path,thred,t.name, paste0('intensities_', module, '.rds'))
    volume_file_path <- file.path(output_path,thred,t.name, paste0('volume_', module, '.rds'))

    # File processing operations if files do not exist
    file_path <- file.path(par.path, thred, "map",paste0("Par_intensities_map_", module, ".rds"))
    par_map <- readRDS(file_path)

    # Initialize matrices to store results
    count_voxel_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(par_map))
    intensity_mean_matrix <- matrix(NA, nrow = nrow(data), ncol = nrow(par_map))

    results <- lapply(1:nrow(par_map), function(i) {
      mapping <- par_map[i, 'mapping'][[1]]
      condition <- data_tissue[, mapping, drop = FALSE] == tissue_type

      intensity_mean <- rowMeans(data[, mapping, drop = FALSE] * condition, na.rm = TRUE)
      count_voxel <- rowSums(condition)

      return(list(intensity_mean, count_voxel))
    })

    for (i in 1:length(results)) {
      intensity_mean_matrix[, i] <- results[[i]][[1]]
      count_voxel_matrix[, i] <- results[[i]][[2]]
    }

    colnames(intensity_mean_matrix) <- par_map[, 'variable']
    colnames(count_voxel_matrix) <- par_map[, 'variable']

    intensity_df <- as.data.frame(intensity_mean_matrix)
    row.names(intensity_df)<- img_name
    volume_df <- as.data.frame(count_voxel_matrix)
    volume_df$fname <- img_name

    #Be careful here:the order of img_name in volume_df may change as a result of the merge operation. This is because the merge function
    #in R can reorder rows based on the common key column used for merging (in this case, the fname column).
    brain_volume_path = file.path(main_dir,"volumes","Bvolumes.rds")
    if (file.exists(brain_volume_path)) {
      brainv_df <- readRDS(file = brain_volume_path)
      # First, merge the data frames on the "fname" column to align the "unit_voxel" values with the corresponding rows in volume_df
      #merged_df <- merge(volume_df, brainv_df, by = "fname")
      merge_df <- inner_join(volume_df, brainv_df, by = "fname")
      # Multiply each column in volume_df (except "fname") by the "unit_voxel" column
      cols_to_multiply <- setdiff(names(volume_df), "fname")
      merged_df[cols_to_multiply] <- merged_df[cols_to_multiply] * merged_df$unit_voxel
      # Resulting data frame has the product, and you might want to remove the extra unit_voxel column if not needed anymore
      final_df <- merged_df[, !names(merged_df) %in% c("unit_voxel", "brainVolume")]
      row.names(final_df) = final_df$fname
      final_df$fname = NULL
      saveRDS(final_df, file = volume_file_path)
    } else {
      row.names(volume_df) = volume_df$fname
      volume_df$fname = NULL
      saveRDS(volume_df, file = volume_file_path)
    }

    saveRDS(intensity_df, file = intensity_file_path)
    message(paste0(ROI, "......Processed and saved data for threshold: ", thred, "......", "tissue: ", t.name))
  }
  rslts <- lapply(thresh_vec, process_threshold)
}

##Test
# main_dir = '/Users/jinyaotian/Downloads/whims_test'
# roi = 'inferior_frontal_gyrus_left'
# dep_list_path = file.path(main_dir, "dep_list", paste0(roi,'.rds'))
# listes <- readRDS(dep_list_path)
# cat("tissue segmentation for all sublists based on super partitions.\n")
# lapply(listes, function(sub_list) {
#   tissue_segment(liste = sub_list,thresh_vec = 0.8,tind = 5, tissue_type =2, main_dir = main_dir )
# })

#Combine by tissue type for each threshold and roi
Cmb_tissue_type <- function(thresh, roi, tissue_type, main_dir){
  modify_data <- function(data, roi_name) {
    colnames(data) <- paste(roi_name, colnames(data), sep = "_")
    return(data)
  }
  read_rds_file <- function(file_path) {
    # Extract the file name
    doc_name <- basename(file_path)
    # Extract the number from the file name using regex
    number <- gsub("\\D", "", doc_name)  # Remove non-numeric characters
    # Read the .rds file
    df <- readRDS(file_path)
    # Update the column names based on whether they contain "reduced_var"
    colnames(df) <- sapply(colnames(df), function(col_name) {
      if (grepl("reduced_var", col_name)) {
        return(paste0("module", number, "_", col_name))
      } else {
        return(col_name)
      }
    })
    return(df)  # Return the modified dataframe
  }

  if (tissue_type == 3){
    tissue.prefix = 'WM'
  } else{
    if (tissue_type == 2){
      tissue.prefix = 'GM'
    }else{
      tissue.prefix  = 'CFS'
    }
  }
  print(paste("Processing for threshold:", thresh, "and tissue type:", tissue.prefix))

  directory = file.path(main_dir, "partition", roi, thresh, tissue.prefix)

  intensity_files_with_prefix <- list.files(path = directory, pattern = paste0("^intensities_"), full.names = TRUE)
  volume_files_with_prefix <- list.files(path = directory, pattern = paste0("^volume_"), full.names = TRUE)

  data_frames_intensity <- lapply(intensity_files_with_prefix, read_rds_file)
  data_frames_volume <- lapply(volume_files_with_prefix, read_rds_file)

  processed_data_frames_intensity <- lapply(data_frames_intensity, function(ls_) {
    df <- ls_
    modify_data(df, roi_name = roi)
  })

  # Combine all processed data frames by columns
  combined_data_frame_intensity <- do.call(cbind, processed_data_frames_intensity)

  processed_data_frames_volume <- lapply(data_frames_volume, function(ls_) {
    df <- ls_
    row.names(df) <- df$fname
    df$fname <- NULL
    modify_data(df, roi_name = roi)
  })

  # Combine all processed data frames by columns
  combined_data_frame_volume <- do.call(cbind, processed_data_frames_volume)
  save.dir = file.path(main_dir, "partition", roi, thresh, tissue.prefix, "cmb")
  if(!dir.exists(save.dir)){dir.create(save.dir, recursive = TRUE)}
  saveRDS(combined_data_frame_intensity, file = file.path(save.dir,"intensities.rds"))
  saveRDS(combined_data_frame_volume, file = file.path(save.dir,"volume.rds"))
}

##Test
#Cmb_tissue_type(thresh = 0.8, roi = roi, tissue_type = 2, main_dir = main_dir)

#Process independent variables from Super-Partition
process_indep_variables <- function(indep_list, tissue_type, tind, roi, main_dir) {
  if (tissue_type == 3){
    t.name = 'WM'
  } else{
    if (tissue_type == 2){
      t.name = 'GM'
    }else{
      t.name = 'CFS'
    }
  }

  # Initialize empty lists for intensity and volume data
  intense_list <- list()
  volume_list <- list()

  # Iterate over each sublist in data
  for (i in seq_along(indep_list)) {
    sublist <- indep_list[[i]]
    ROI <- sublist$ROI
    Variables <- sublist$Variables

    # Skip if Variables is empty
    if (length(Variables) == 0) next

    # Process intensity data
    intensity_file <- file.path(main_dir, "intensities",paste0("intensities_", ROI,'_',tind, ".rds"))
    if (file.exists(intensity_file)) {
      intensity_data <- readRDS(intensity_file)
      selected_intensity_data <- intensity_data[,Variables,drop = FALSE]
      selected_intensity_data <- selected_intensity_data[-c(1:3), , drop = FALSE]
      colnames(selected_intensity_data) <- paste0(ROI, "_", colnames(selected_intensity_data))
      intense_list[[length(intense_list) + 1]] <- selected_intensity_data
      rm(selected_intensity_data)
    }

    # Process volume data
    volume_file <- file.path(main_dir, "tissues",paste0("tissues_", ROI,'_',tind, ".rds"))
    if (file.exists(volume_file)) {
      volume_data <- readRDS(volume_file)
      selected_volume_data <- volume_data[,Variables,drop = FALSE]
      selected_volume_data <- selected_volume_data[-c(1:3), , drop = FALSE]
      condition <- selected_volume_data == tissue_type
      if(file.exists(file.path(main_dir,"volumes", "Bvolumes.rds"))){
        volumes = readRDS(file.path(main_dir,"volumes", "Bvolumes.rds"))
        unit_ = volumes$unit_voxel
      }else{
        unit_ = 0.001
      }
      selected_volume_data <- selected_volume_data*condition*unit_ #cm^3
      colnames(selected_volume_data) <- paste0(ROI, "_", colnames(selected_volume_data))
      volume_list[[length(volume_list) + 1]] <- selected_volume_data
      rm(selected_volume_data)
    }
  }

  # Combine lists into data frames
  intense_dat <- bind_cols(intense_list)
  volume_dat <- bind_cols(volume_list)
  save_dir = file.path(main_dir, "indep_variables", roi,t.name )
  if(!dir.exists(save_dir)){dir.create(save_dir, recursive = TRUE)}
  saveRDS(intense_dat, file = file.path(save_dir,paste0("intensities.rds")))
  saveRDS(volume_dat, file = file.path(save_dir, paste0("volume.rds")))
}

#Combine independent variables with reduced variables by tissue type by roi
Cmb_indep_with_dep <- function(tissue_type, roi, thresh, main_dir){
  if (tissue_type == 3){
    t.name = 'WM'
  } else{
    if (tissue_type == 2){
      t.name = 'GM'
    }else{
      t.name = 'CFS'
    }
  }
  iten_path_reduced = file.path(main_dir, "partition", roi, thresh,t.name, "cmb", "intensities.rds")
  vol_path_reduced = file.path(main_dir, "partition", roi, thresh,t.name, "cmb", "volume.rds")
  iten_path_indep = file.path(main_dir, "indep_variables",roi, t.name, "intensities.rds")
  vol_path_indep = file.path(main_dir, "indep_variables",roi, t.name, "volume.rds")

  # Check if any of the files do not exist
  if (!all(file.exists(iten_path_reduced, vol_path_reduced, iten_path_indep, vol_path_indep))) {
    stop("Error: Missing combining steps")
  }else{
    iten_reduced = readRDS(iten_path_reduced)
    iten_indep = readRDS(iten_path_indep)
    iten_cmb = do.call(cbind, list(iten_reduced,iten_indep))
    saveRDS(iten_cmb, file = file.path(main_dir, "partition", roi, thresh,t.name, "cmb", "intensities_whole.rds"))

    vol_reduced = readRDS(vol_path_reduced)
    vol_indep = readRDS(vol_path_indep)
    vol_cmb = do.call(cbind, list(vol_reduced, vol_indep))
    saveRDS(vol_cmb, file = file.path(main_dir, "partition", roi, thresh,t.name, "cmb", "volume_whole.rds"))
  }
}

##Test
#Cmb_indep_with_dep(tissue_type = 2, roi = "inferior_frontal_gyrus_left", thresh = 0.8, main_dir = '/Users/jinyaotian/Downloads/whims_test')

#Combine by ROI
concat_reduced_var <- function(roi, thresh, main_dir){
  #Structure of dir: /roi/thresh
  dir = file.path(main_dir, roi, thresh)
  files_intensity <- list.files(path = file.path(dir, "intensities"))
  files_map <- list.files(path = file.path(dir, "map"))
  data_frames_intensity <- lapply(files, function(file) {
    module <- str_extract(file, "(?<=Par_intensities_).*(?=_.rds)")
    df <- readRDS(file)
    cols_to_rename <- names(df)[grepl("^reduced", names(df))]
    new_names <- paste(module, cols_to_rename, sep = "_")
    names(df)[grepl("^reduced", names(df))] <- new_names
    return(df)
  })
  combined_df_intensity <- do.call(cbind, data_frames_intensity)
  saveRDS(combined_df_intensity, file = file.path(dir, "reduced_intensities.rds"))

  data_frames_map <- lapply(files, function(file) {
    module <- str_extract(file, "(?<=Par_intensities_map_).*(?=_.rds)")
    df <- readRDS(file)
    rows_to_rename <- df$variable[grepl("^reduced", df$variable)]
    new_names <- paste(module, rows_to_rename, sep = "_")
    df$variable[grepl("^reduced", df$variable)] <- new_names
    return(df)
  })
  combined_df_map <- do.call(rbind, data_frames_map)
  saveRDS(combined_df_map, file = file.path(dir, "reduced_intensities_map.rds"))
}
