## what the heck does this do?
## I think the point is to extact/construct identifiers to use as row names.
rnames <- function(lst, item) {
  f <- function(ll, inds) {
    if ((ii <- match(item, names(ll), FALSE)))
      list(inds=c(inds, ii), len=length(ll[[ii]]))
    else if (all(is.atomic(unlist(ll, FALSE))) || !is.list(ll))
      NULL
    else
      lapply(seq_along(ll), function(i) f(ll[[i]], inds=c(inds, i)))
  }
  unlist(f(lst, NULL))
}

extractOneLGF <- function(J) {
  bands <- J$iscn2016_bands # cytoband labels/names
  ## extract the status
  OUT <- J$output
  Status <- sapply(OUT, "[", "status") # pulls out all the status information
  ## extract the (sub)clone ids
  clone <- unlist(lapply(OUT, function(x){
    if((x$status)=="Success"){
      1:length(lapply(x$parsing_result, "[[", "loss_gain_fusion_computing"))
    }
  }))
  ## extract the binary-mapped LGF data
  lgf <- lapply(OUT, function(x) { # Get the loss data from parsing_result
    if((x$status) == "Success"){
      lgf <- lapply(x$parsing_result, "[[","loss_gain_fusion_computing")
    }
  })
  ## compute row-name IDS
  df.ID <- as.data.frame(matrix(rnames(lgf, "loss"), ncol = 4, byrow = TRUE))
  id <- df.ID$V1
  if(any(df.ID$V2 != clone)) {
    warning("Disagreement among clone IDs.")
  }

  fullID <- apply(df.ID[, 1:3], 1, paste, collapse=".")
  if (any(duplicated(fullID))) {
    warning("Duplicated identifiers should not happen.")
  }
  ## Extract the loss, gain, and fusion binary data separately
  loss1 <- sapply(lgf, function(x){
    lapply(x, "[[","loss")
  })
  df_Loss <- as.data.frame(matrix(unlist(loss1), ncol = 916, byrow = TRUE),
                           stringsAsFactors = FALSE)
  colnames(df_Loss) <- paste("Loss", bands, sep = "_")

  gain1 <- sapply(lgf, function(x){
    lapply(x, "[[","gain")
  })
  df_Gain <- as.data.frame(matrix(unlist(gain1), ncol = 916, byrow = TRUE),
                           stringsAsFactors = FALSE)
  colnames(df_Gain) <- paste("Gain", bands, sep = "_")
  fusion1 <- sapply(lgf, function(x){
    lapply(x, "[[","fusion")
  })
  df_Fusion <- as.data.frame(matrix(unlist(fusion1), ncol = 916, byrow = TRUE),
                             stringsAsFactors = FALSE)
  colnames(df_Fusion) <- paste("Fusion", bands, sep = "_")
  ## Put them back together
  df.lgf <- cbind(df_Loss, df_Gain, df_Fusion)
  rownames(df.lgf) <- fullID
  df.lgf$ID = id
  df.lgf$Clones = clone
  w <- which(colnames(df.lgf) == "ID")
  df.lgf <- df.lgf[, c(w:(w + 1), 1:(w-1))]
  df.lgf
}

readLGF <- function(files = NULL, folder = NULL) {
  ## Figure out which files we want to read
  if (is.null(folder)) {
    folder  <-  "."
  }
  if (is.null(files)) {
    files  <- list.files(folder, pattern = "*.json", full.names = TRUE) # character vector , one file per entry
  }
  if (length(files) < 1) {
    stop("No JSON input files to read.")
  }
  message("Reading ", length(files), " file(s) from '", folder, "'.\n")

  ## make sure that, when done,  we leave the working directory in the same state that we found it.
  home <- getwd()
  on.exit(setwd(home))
  setwd(folder)

  ## It might get large -- not said by Davis Guggenheim
  myJSON <- lapply(files, function(x) fromJSON(file = x)) # a list, one element per JSON file
  temp  <- lapply(myJSON, extractOneLGF)
  temp
}
