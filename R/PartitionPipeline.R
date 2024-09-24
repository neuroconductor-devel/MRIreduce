#' Partition Pipeline for Image Analysis
#'
#' This R6 class is designed to streamline the processing pipeline for image analysis,
#' including steps from initial processing to combining independent variables with
#' reduced variables by tissue type by ROI.
#'
#' @import R6
#' @importFrom parallel makeCluster parLapply stopCluster clusterExport clusterEvalQ
#' @export
PartitionPipeline <- R6Class("PartitionPipeline",
                           public = list(
                             tind = NULL,
                             nfl = NULL,
                             main_dir = NULL,
                             tissue_type = NULL,
                             outp_volume = TRUE,
                             ICC_thresh_vec = NULL,
                             suppar_thresh_vec = seq(0.7, 1, 0.01),
                             B = 2000,
                             roi = NULL,
                             num_cores = NULL,
                             # Constructor
                             initialize = function(tind = NULL, nfl = NULL, main_dir = NULL, tissue_type = NULL,
                                                   outp_volume = TRUE, ICC_thresh_vec = NULL, suppar_thresh_vec = seq(0.7, 1, 0.01), B = 2000, roi = NULL, num_cores = NULL) {
                               self$tind <- tind
                               self$nfl <- nfl
                               self$main_dir <- main_dir
                               self$tissue_type <- tissue_type
                               self$outp_volume <- outp_volume
                               self$ICC_thresh_vec <- ICC_thresh_vec
                               self$suppar_thresh_vec <- suppar_thresh_vec
                               self$B <- B
                               self$roi <- roi
                               self$num_cores = num_cores
                             },

                             # Step1: Process images and extract intensities and tissues
                             iproc = function() {
                               cat("Starting image processing for ROI:", self$tind, "\n")
                               result = iproc(tind = self$tind, nfl = self$nfl, main_dir = self$main_dir, outp_volume = self$outp_volume)
                               self$roi = result
                             },

                             # Step2: Apply super partition analysis
                             supparfun = function() {
                               cat("Starting applying Super Partition to ROI:", self$roi, "\n")
                               supparfun(
                                  tind = self$tind,
                                  roi = self$roi,
                                  thresh_vec = self$suppar_thresh_vec,
                                  B = self$B,
                                  main_dir = self$main_dir)
                             },

                             # Step3: Map super partitions to ROI
                             map_suppar_roi = function() {
                               map_suppar_roi(roi = self$roi, main_dir = self$main_dir)
                             },

                             # Step4: Partition intensities based on super partitions using lapply
                             partition_intensities = function() {
                               dep_list_path = file.path(self$main_dir, "dep_list", paste0(self$roi,'.rds'))
                               listes <- readRDS(dep_list_path)
                               cat("Partitioning intensities for all sublists based on super partitions.\n")
                               cl <- makeCluster(self$num_cores)
                               # Export necessary functions and variables
                               clusterExport(cl, varlist = c("parfun", "self"), envir = environment())
                               # Load required packages on worker nodes
                               clusterEvalQ(cl, {library(MRIreduce); library(partition)})
                               parLapply(cl, listes, function(sub_list) {
                                 parfun(liste = sub_list,tind = self$tind, thresh_vec = self$ICC_thresh_vec, main_dir = self$main_dir)
                               })
                               # Stop the cluster after use
                               on.exit(stopCluster(cl))
                             },

                             # Step5: Segment tissue types
                             tissue_segment = function() {
                               dep_list_path = file.path(self$main_dir, "dep_list", paste0(self$roi,'.rds'))
                               listes <- readRDS(dep_list_path)
                               cat("tissue segmentation for all sublists based on super partitions.\n")
                               lapply(listes, function(sub_list) {
                                 tissue_segment(liste = sub_list,thresh_vec = self$ICC_thresh_vec,tind = self$tind, tissue_type =self$tissue_type, main_dir = self$main_dir)
                               })
                             },

                             # Step6: Combine by tissue type for each threshold and roi
                             Cmb_tissue_type = function(){
                               lapply(as.list(self$ICC_thresh_vec), function(x){
                                 Cmb_tissue_type(thresh = x, roi = self$roi, tissue_type = self$tissue_type, main_dir = self$main_dir)
                               })
                             },

                             # Step7: Process independent variables from Super-Partition
                             process_indep_variables = function(){
                               indep_list_path = file.path(self$main_dir, "indep_list", paste0(self$roi,'.rds'))
                               indep_list <- readRDS(indep_list_path)
                               process_indep_variables(indep_list = indep_list, tissue_type = self$tissue_type, tind = self$tind,
                                                       roi = self$roi, main_dir = self$main_dir)
                             },

                             # Step8: Combine independent variables with reduced variables by tissue type by roi
                             Cmb_indep_with_dep = function(){
                               lapply(as.list(self$ICC_thresh_vec), function(x){
                                 Cmb_indep_with_dep(thresh = x, roi = self$roi, tissue_type = self$tissue_type, main_dir = self$main_dir)
                               })
                             }
                           )
)

#test
# pipeline <- PartitionPipeline$new(tind = 5, nfl = list.files('/Users/jinyaotian/Downloads/pipeline_test/eve_t1', full.names = TRUE),
#                                   main_dir = "/Users/jinyaotian/Downloads/pipeline_test",
#                                   tissue_type = 2,
#                                   ICC_thresh_vec = c(0.8, 0.9),
#                                   roi = "inferior_frontal_gyrus_left",
#                                   num_cores = 5)
#
# # Step 1
# pipeline$iproc()
# # Step 2
# pipeline$supparfun()
# # Step 3
# pipeline$map_suppar_roi()
# # Step 4
# pipeline$parfun()
# # Step 5
# pipeline$tissue_segment()
# # Step 6
# pipeline$Cmb_tissue_type()
# # Step 7
# pipeline$process_indep_variables()
# # Step 8
# pipeline$Cmb_indep_with_dep()

