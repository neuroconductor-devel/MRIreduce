#' Plot Mask on EVE Template using Python and Nilearn
#'
#' This function plots a mask image on the EVE template using Python's Nilearn library and optionally saves the plot as an image.
#' It ensures the required Python libraries are installed, and it handles the appropriate Conda environment setup.
#'
#' @param mask_img_path A string representing the file path to the mask NIfTI file.
#' @param cmap The colormap to use. Either a string (name of a matplotlib colormap) or a matplotlib colormap object. Default is 'bwr_r'.
#' @param alpha Transparency level for the overlay. Default is 1.
#' @param save_path A string representing the file path where the plot image should be saved (e.g., "output.png"). If NULL, the plot will be shown interactively instead of saved. Default is NULL.
#' @param ... Additional arguments to customize the anatomical plot. These arguments are passed
#'        directly to the Python function `nilearn.plotting.plot_anat`. You can specify parameters
#'        such as `title`, `display_mode`, `cut_coords`, `dim`, etc.
#'        For more details on available options, refer to the official documentation at:
#'        https://nilearn.github.io/stable/modules/generated/nilearn.plotting.plot_anat.html
#' @details
#' The function first detects the system architecture (ARM or x86) and ensures the appropriate Conda environment is set up (Miniforge for ARM or Miniconda for x86).
#' It then checks and installs the required Python libraries: `nilearn`, `nibabel`, and `matplotlib`.
#' The function loads the Python environment and calls the `plot_mask_on_eve` function, which plots the mask on the EVE template. If `save_path` is provided, the plot is saved as an image at the specified location.
#'
#' @return None
#' @importFrom reticulate source_python
#' @examples
#' \dontrun{
#' map2_eve(
#'   mask_img_path = "/path/to/mask_nifti_GM_Volume.nii.gz",
#'   save_path = "/path/to/save_output.png",
#'   cmap = "bwr_r",
#'   alpha = 0.8,
#'   title = "Mask on EVE Template"
#' )
#' }
#'
#' @export
map2_eve <- function(mask_img_path, cmap = 'bwr_r', alpha = 1, save_path = NULL, ...) {
  # Function to install and use the correct Conda environment (Miniforge for ARM, Miniconda for x86)
  install_and_use_conda <- function() {
    # Detect system architecture (ARM vs. x86)
    system_arch <- Sys.info()["machine"]
    conda_path <- NULL  # Initialize conda path

    if (grepl("arm64", system_arch) || grepl("aarch64", system_arch)) {
      # ARM-based systems (e.g., Apple M1/M2, Raspberry Pi)
      conda_path <- "~/Library/r-miniconda-arm64"

      # Check if Miniforge is already installed and available
      if (!reticulate::py_available(initialize = FALSE) || !grepl("arm64|aarch64", reticulate::py_config()$version)) {
        message("ARM-compatible Python environment not found. Installing Miniforge (ARM-compatible Python)...")
        reticulate::install_miniconda(path = conda_path, update = TRUE, force = TRUE)
      } else {
        message("ARM-compatible Python environment is already available.")
      }

    } else if (grepl("x86_64", system_arch)) {
      # x86_64-based systems (Intel/AMD)
      conda_path <- "~/Library/r-miniconda-x86"

      # Check if Miniconda is already installed and available
      if (!reticulate::py_available(initialize = FALSE)) {
        message("x86-compatible Python environment not found. Installing Miniconda...")
        reticulate::install_miniconda(path = conda_path, update = TRUE, force = TRUE)
      } else {
        message("x86-compatible Python environment is already available.")
      }

    } else if (grepl("ppc64le", system_arch)) {
      # PowerPC 64-bit Little Endian systems
      stop("PowerPC architecture detected. Please manually configure your environment as Conda support is limited.")
    } else {
      stop("Unsupported system architecture: ", system_arch)
    }

    # Use the correct conda environment
    reticulate::use_condaenv(condaenv = conda_path, required = TRUE)
  }

  # Ensure that Python libraries are installed
  install_python_dependencies <- function() {
    # Check if 'nilearn' is available
    if (!reticulate::py_module_available("nilearn")) {
      message("Installing 'nilearn' Python package...")
      reticulate::py_install("nilearn")
    } else {
      message("'nilearn' is already installed.")
    }

    # Check if 'nibabel' is available
    if (!reticulate::py_module_available("nibabel")) {
      message("Installing 'nibabel' Python package...")
      reticulate::py_install("nibabel")
    } else {
      message("'nibabel' is already installed.")
    }
    # Check if 'matplotlib' is available (required for nilearn plotting)
    if (!reticulate::py_module_available("matplotlib")) {
      message("Installing 'matplotlib' Python package...")
      reticulate::py_install("matplotlib")
    } else {
      message("'matplotlib' is already installed.")
    }
  }

  # Install and set up appropriate Conda environment based on architecture
  install_and_use_conda()
  # Ensure Python packages are installed
  install_python_dependencies()

  # Get the path to the template image (eve_t1.nii.gz)
  template_img_path <- system.file("extdata", "eve_t1.nii.gz", package = "MRIreduce")

  if (template_img_path == "") {
    stop("The template image 'eve_t1.nii.gz' could not be found.")
  }

  # Load the Python environment
  reticulate::source_python(system.file("python/plot_mask_on_eve.py", package = "MRIreduce"))

  # Call the Python function with the template image path
  plot_mask_on_eve(mask_img_path = mask_img_path, template_img_path = template_img_path, cmap = cmap, alpha = alpha,save_path = save_path, ...)
}

##test
#map2_eve(mask_img_path = "/Users/jinyaotian/Downloads/mask_nifti_GM_Volume_pm25_test1_avg_red.nii.gz",save_path = '/Users/jinyaotian/Desktop/test.png',title = "Mask on EVE Template" )
