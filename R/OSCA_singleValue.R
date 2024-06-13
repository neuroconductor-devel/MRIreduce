#New version of covariate analysis
#Programmer: Jinyao Tian
#Date: 04/18/2024

options(stringsAsFactors=FALSE)
# .libPaths("/project/millstei_783/millstei/Rlibs")
library(haven)
library(dplyr)
library(modACDC)

coVar_analysis <- function(data.path,perm = FALSE,covar.file = '/users/jinyaotian/Downloads/vars_df.rds', link.path = '/users/jinyaotian/Downloads/whimsID_link.csv'){
  #permutation (perm): No: NUll, Yes: 'perm'
  df <- readRDS(data.path)
  #covariate data
  vars_df <- readRDS(covar.file)

  df$img <- rownames(df)
  df <- df %>%
    mutate(ID = as.numeric(substr(sub("-.*", "", img),1,7)), Date = as.Date(sub(".*-(.*)_(.*)", "\\1", img), format = "%Y%m%d"))
  #Linking WHI_ID (WHI common ID in all WHI data sets) & UPenn_ID (ID shown in scan file/folder names)
  link <- read.csv(link.path)
  #left join
  result <- merge(df, link, by.x = 'ID', by.y = 'whims_id', all.x = TRUE)
  Age_df <- vars_df[,!colnames(vars_df) %in% c('MRI1_whi_day')]
  Age_df$age <- Age_df$age_at_MRI1
  #Inner join
  result <- merge(result, Age_df, by.x = 'WHI_ID', by.y = 'ID', all.x = FALSE)
  #Order the result by Whi ID and then by Date
  result_ordered <- result[order(result$WHI_ID, result$Date), ]
  for (j in 2:nrow(result_ordered)){
    if (result_ordered[j,'WHI_ID'] == result_ordered[j-1,'WHI_ID']){
      #Calculate the age at MRI2 (follow up visit)
      time_difference = as.numeric(result_ordered[j,'Date'] - result_ordered[j-1,'Date'])/365.0
      result_ordered[j,'age'] = result_ordered[j,'age'] + time_difference
    }
  }
  co_vars <- c('WHI_ID', 'ID', 'img', 'Date', colnames(Age_df))
  Full_Data <- result_ordered[, !names(result_ordered) %in% co_vars]
  NumCovars <- result_ordered[, c("PSHTDEP_logit", "nSES_MRI1_on_CT")]
  CatCovars <- result_ordered[, c('BMI4_MRI1', 'SMOKING_MRI1', 'ALCOHOL_cat_MRI1', 'MODSTRENEX_MRI1', 'HORMSTAT2', 'REGION',
                                  'HRTARM', 'MRI_clinic', 'EMPLOY2', 'ApoE_imputed','CVD_risk_MRI1_2', 'stroke_by_MRI1', 'Scanner_Type', 'race_eth_NIH',
                                  'EDUCG3', 'INCOME_CAT','ApoE_e4_positive_new', 'PD_bf_MRI1', 'MCI_PD_bf_MRI1')]
  if(perm == 'perm'){ExternalVar <- sample(result_ordered[, 'age'])}else{
    ExternalVar <- result_ordered[, 'age']
  }

  require(modACDC)
  t_out <- OSCA_singleValue(df = Full_Data,
                            externalVar = ExternalVar,
                            # numCovars = NumCovars,
                            # catCovars= CatCovars,
                            oscaPath = '/Users/jinyaotian/Downloads/osca_Mac_v0.45/osca_Mac')

  variance_explained<- t_out$Variance
  message(basename(data.path),'......Covariate Analysis of Aging has been successfully done......')
  return(list(file = basename(data.path), Variance_explained_by_age = variance_explained))
}

if(0){
  directory = '/project/millstei_783/jinyaoti/image_preprocess/whims1/cmb'
  files_with_pattern <- list.files(path = directory,
                                   pattern = "^combined_volume_[0-9]+(\\.[0-9]+)?_GM\\.rds$",
                                   full.names = TRUE)
  rslts <- lapply(as.list(files_with_pattern), function(path){
    coVar_analysis(data.path = path)
  })
}
#Test
if(1){
  rslt_08 <- coVar_analysis(data.path = '/Users/jinyaotian/Downloads/combined_intense_0.8_WM.rds')
}
