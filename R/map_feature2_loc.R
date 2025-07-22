#' Map Reduced Feature Back to Brain Image Voxel Locations
#'
#' This function maps a given feature name to voxel locations by identifying the pattern in the feature name.
#' If the feature name contains "reduced_var", it follows the format "Roi_module_reduced_var" and this indicates that it is a reduced feature. Otherwise, it follows
#' the format "Roi_Vnumber", which is not a reduced feature by Partition. Based on the pattern, the function extracts the appropriate region of interest (ROI) and
#' retrieves voxel locations from intensity data.
#'
#' @param feature_name String. A feature name.
#' @param threshold a numeric between 0 and 1. The Partition threshold value applied to the data.
#' @param main_dir String. The main directory containing the data files.
#'
#' @return A data frame containing the voxel locations (x,y,z Coordinates) corresponding to the extracted mapping column.
#'
#' @details
#' Each voxel in a brain image corresponds to a specific feature, establishing a one-to-one mapping between a voxel and a feature. This relationship allows us
#' to localize particular brain features to specific brain areas at the voxel level, enabling visualization of these features. However, after applying Super-Partition and Partition techniques,
#' multiple brain features are aggregated into a single reduced feature as part of a data reduction process. The goal of the following function is to relate the reduced feature back to its component features,
#' and subsequently identify the brain image voxels to which those component features are mapped.
#' This function performs the following steps:
#'
#' 1. **Check Feature Name Format:** It first checks if the `feature_name` contains "reduced_var". If it does, the function assumes the format "Roi_module_reduced_var", otherwise it assumes the format "Roi_Vnumber".
#'
#' 2. **Extract ROI and Mapping Information:**
#'    - For the "Roi_module_reduced_var" format, the ROI and module number are extracted, and the corresponding partition file is read to retrieve the mapping vector.
#'    - For the "Roi_Vnumber" format, the `V` number is extracted from the feature name and used as the mapping vector.
#'
#' 3. **Extract Voxel Locations:** The mapping vector is then used to extract voxel locations from the intensity data.
#'
#' @examples
#' \dontrun{
#' loc_df <- map_feature2_loc(feature_name = "inferior_frontal_gyrus_left_module4_reduced_var_13",
#'   threshold = 0.8, main_dir = "/path/to/data")
#' }
#'
#' @export
map_feature2_loc <- function(feature_name, threshold,main_dir) {
  lab_df <- NULL
  if (grepl("reduced_var", feature_name)) {
    # Extract the ROI for the format "Roi_module_reduced_var"
    roi <- sub("_module.*", "", feature_name)  # Remove everything after and including "_module"
    module = sub(".*module(\\d+).*", "\\1", feature_name)
    variable_name = sub(".*(reduced_var_\\d+).*", "\\1", feature_name)
    file_path = file.path(main_dir, "partition", roi, threshold, "map", paste0("Par_intensities_map_", module, ".rds"))
    df = readRDS(file_path)
    vec <- df[df$variable == variable_name, "mapping"][[1]]
  } else {
    # Extract the ROI for the format "Roi_Vnumber"
    roi <- sub("_V\\d+.*", "", feature_name)  # Remove everything after and including "_V<number>"
    vec <- sub(".*(V\\d+).*", "\\1", feature_name)
  }

  # Extract voxel locations
  extract_mapping_voxel <- function(tind, roi, Cols) {
    file_name <- paste0("intensities_", roi, "_", tind,".rds")
    df <- readRDS(file.path(main_dir, "intensities", file_name))
    locs <- df[1:3, Cols, drop = FALSE]
    return(locs)
  }
  load_required_data("eve_label_info_dataframe.RData", "whims")
  tind = lab_df[lab_df$text_label == roi, "integer_label"]
  # Extract the voxel locations
  loc_df <- extract_mapping_voxel(tind, roi, vec)
  return(loc_df)
}

##Test
# feature_name = "inferior_frontal_gyrus_left_module11_reduced_var_5"
# threshold = 0.8
# main_dir = '/Users/jinyaotian/Downloads/whims_test'
# loc_df = map_feature2_loc(feature_name = feature_name, threshold = threshold, main_dir = main_dir)
