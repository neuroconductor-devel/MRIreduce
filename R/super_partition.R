#' Superimpose Partition
#'
#' Assign features to groups in an agglomerative approach according to correlations
#' between pairs of features with constraints on group size and feature group membership.
#'
#' @param tmp DataFrame of features.
#' @param thresh Numeric vector specifying correlation thresholds.
#' @param n.chunkf Integer, the number of features to process in each chunk.
#' @param B Integer, the maximum group size.
#' @param compute.corr Logical, whether to compute correlations.
#' @param dist.thresh Numeric, distance threshold for including correlations.
#' @param dir.tmp Character, path to a temporary directory for storing intermediate results.
#'
#' @return A list containing two elements: a list of feature groups and a vector of independent features.
#' @export
#' @examples
#' \dontrun{
#'   result <- suppar(tmp = myData, thresh = c(0.5, 0.75), n.chunkf = 100,
#'                    B = 2000, compute.corr = TRUE, dist.thresh = NULL,
#'                    dir.tmp = "temp_corr_temp")
#' }
#' @export
suppar = function(tmp, thresh=NULL, n.chunkf=10000, B = 2000, compute.corr=TRUE, dist.thresh=NULL, dir.tmp){
  cchunk.l = 10^6 # chunk processing of correlation file datasets, number of rows
  if(is.null(colnames(tmp))) colnames(tmp) = paste0("V", 1:ncol(tmp))
  fnms.indep = colnames(tmp) # all features start classified as independent
  fnms.dep = NULL
  #dir.tmp = paste(getwd(), "temp_corr_temp", sep="/") # temporary directory
  if(compute.corr){
    dir.create(dir.tmp) # create temporary directory
    print(paste("created:", dir.tmp))
    corfun2(tmp, thresh, n.chunkf, dir.tmp, dist.thresh=dist.thresh) # compute correlations and save as temporary files
  }

  grpl = vector('list', 1)
  thresh1 = getthresh(dir.tmp) # vector of correlation bounds for files, largest to smallest
  # loop across correlation files
  # read corrbin with largest correlations first
  for(fl.ind in 1:length(thresh1)){
    fnm = paste0(dir.tmp, "/corrbin_", thresh1[fl.ind], "_.Rdata") # file path for thresh1[fl.ind]
    load(fnm) # corrbin
    if(nrow(corrbin) > 0){
      corrbin = corrbin[order(corrbin[,"value"], decreasing=TRUE), ]
      corrbin = corrbin[!duplicated(corrbin[, c("Var1", "Var2")]), ]
      if(length(grpl[[1]]) < 1){
        grpl[[1]] = as.character(corrbin[1 ,c("Var1", "Var2")])
        corrbin = corrbin[-1, ]
      }
      n.cchunk = ceiling(nrow(corrbin) / cchunk.l) # chunk processing of correlation file datasets
      for(cchunk in 1:n.cchunk){
        startc = (cchunk - 1) * cchunk.l + 1
        endc = startc + cchunk.l - 1
        if(endc > nrow(corrbin)) endc = nrow(corrbin)
        cbtemp = corrbin[startc:endc, ]
        grpl = partagg(grpl, cbtemp, B)
        # check for duplicates and remove
        for(j in 1:length(grpl)){
          grpl[[j]] = grpl[[j]][!duplicated(grpl[[j]])]
        }
      }

      fnms.dep = unlist(grpl)
      fnms.indep = fnms.indep[!is.element(fnms.indep, fnms.dep)]
      print(paste("n.modules:", length(grpl), "n.indep.features:", length(fnms.indep),  "n.dep.features:", length(fnms.dep)))
      print(paste("End processing file:", fnm))
    } # if nrow corrbin
  } # end fl.ind loop

  if(compute.corr){
    unlink(dir.tmp, recursive = TRUE)
    unlink(dir.tmp)
  }
  # check for NAs and remove
  if(!is.null(grpl[[1]])){
    lg = rep(NA, length(grpl))
    for(j in 1:length(grpl)){
      grpl[[j]] = grpl[[j]][grpl[[j]] != "NA"]
      lg[j] = length(grpl[[j]])
    }
    grpl = grpl[lg > 0]
  } else {
    grpl = NULL
  } # end check grpl if null
  outp = vector('list', 2)
  outp[[1]] = grpl
  outp[[2]] = fnms.indep
  return(outp)
} # end suppar

dist2mat = function(m1, m2){
  m2 = t(m2)
  tmp = apply(m1, 2, function(v) sqrt(rowSums(m2^2)+sum(v^2)-2*(m2%*%as.matrix(v))[,1]))
  colnames(tmp) = colnames(m1)
  rownames(tmp) = rownames(m2)
  tmp = t(tmp)
  return(tmp)
}


corfun1 = function(dat, dat2=NULL, thresh, dir.tmp, dist.thresh=NULL){
  if(!is.null(dist.thresh)){
    if(is.null(dat2)){
      dmat = as.matrix(dist(t(dat[1:3,])))
    } else {
      dmat = dist2mat(dat[1:3,], dat2[1:3,])
      dat2 = dat2[-c(1:3), ]
    }
    dat = dat[-c(1:3), ]
  } # end if !is.null dist.thresh

  if(is.null(dat2)){
    tmp = cor(dat, use="pairwise.complete.obs")
  } else {
    tmp = cor(dat, dat2, use="pairwise.complete.obs")
  }
  if(!is.null(dist.thresh)){
    tmp = ifelse(dmat < dist.thresh, tmp, 0)
  }
  if(is.null(dat2)){
    tmp1.l = melt(upper.tri(tmp, diag = FALSE)) # logical for upper triangle
    tmp2 = melt(tmp) # long skinny format
    tmp2 = tmp2[tmp1.l[,"value"], ] # keep upper triangle only
  } else {
    tmp2 = melt(tmp) # long skinny format
  } # is null dat2
  rm(tmp)
  tmp2[,"Var1"] = as.character(tmp2[,"Var1"])
  tmp2[,"Var2"] = as.character(tmp2[,"Var2"])
  tmp2 = tmp2[tmp2[,"value"] > thresh[1], ]
  n.edges = nrow(tmp2)

  # save results to files.
  #Threshold Filtering and Data Saving:
  #The correlations are filtered based on the thresh values. The function iterates over the thresh vector, creating subsets of the correlation data (corrbin) that fall between consecutive threshold values.
  #These subsets are then saved into separate files within the specified dir.tmp directory. If a file for a specific threshold range already exists, the new data is appended to it.
  for(j in 1:(length(thresh) - 1)){
    aa = tmp2[,"value"] > thresh[j]
    bb = tmp2[,"value"] <= thresh[j+1]
    corrbin = tmp2[aa & bb, ]
    fnm = paste0(dir.tmp, "/corrbin_", thresh[j], "_.Rdata")
    if(file.exists(fnm)){
      tmp4 = corrbin
      rm(corrbin)
      load(fnm) # load corrbin
      corrbin = rbind(corrbin, tmp4)
      save(corrbin, file=fnm)
    } else {
      save(corrbin, file=fnm)
    }
  }
  #The function returns the number of edges (n.edges), which is essentially the number of rows in the filtered correlation data (tmp2) after applying the first threshold.
  return(n.edges)
}


corfun2 = function(dat, thresh, n.chunkf, dir.tmp, dist.thresh=NULL){
  n.f = ncol(dat)
  n.chunk = ceiling(n.f/n.chunkf) #number of chunk
  size.chunk = ceiling(n.f/n.chunk) #number of features within each chunk
  for(chnk in 1:n.chunk){
    i.start = (chnk - 1) * size.chunk + 1
    i.end = i.start + size.chunk - 1
    if(i.end > n.f) i.end = n.f
    n.edges = corfun1(dat=dat[,i.start:i.end], thresh=thresh, dir.tmp=dir.tmp, dist.thresh=dist.thresh)
    print(paste("Finishing chunk", chnk, "start.ind:", i.start, "end.ind:", i.end, "n.edges:", n.edges))
    if(n.chunk > 1){
      for(chnk2 in (chnk+1):n.chunk){
        if(chnk2 <= n.chunk){
          j.start = (chnk2 - 1) * size.chunk + 1
          j.end = j.start + size.chunk - 1
          if(j.end > n.f) j.end = n.f
          n.edges = corfun1(dat=dat[,i.start:i.end], dat2=dat[,j.start:j.end], thresh=thresh, dir.tmp=dir.tmp, dist.thresh=dist.thresh)
          print(paste("Finishing chunk2", chnk2, "start.ind:", j.start, "end.ind:", j.end, "n.edges:", n.edges))
        }
      }
    }
  }
  return(NULL)
}

# Function to retrieve correlation thresholds
getthresh = function(dir.tmp){
  fl.nms = list.files(dir.tmp)
  thresh1 = NULL
  if(length(fl.nms) > 0){
    tmpl = strsplit(fl.nms, "_", fixed=TRUE)
    thresh1 = rep(NA, length(fl.nms))
    for(j in 1:length(tmpl)) thresh1[j] = as.numeric(tmpl[[j]][2])
    thresh1 = thresh1[order(thresh1, decreasing = TRUE)]
  }
  return(thresh1)
}

# check membership across list element
efun = function(vec, mypair){
  lvec = is.element(mypair, vec)
  s = ifelse(sum(lvec) > 0, TRUE, FALSE)
  return(s)
}


