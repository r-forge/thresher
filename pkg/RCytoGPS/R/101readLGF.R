## I think the point is to extract/construct identifiers to use as row names.
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

## Make up sample names that include row number, padded with leading zeros
padnames <- function(prefix, values, ndigits = NULL) {
  M <- max(values)
  if (is.null(ndigits)) {
    ndigits <- trunc(log10(M))
  }
  padlength <- sapply(values, function(n) sum(n < 10^(1:ndigits)))
  pad <- sapply(padlength, function(p) paste(rep("0", p),collapse = ""))
  paste(prefix, pad, values, sep="")
}


extractOneLGF <- function(J) {
  bands <- J$iscn2016_bands # cytoband labels/names
  ## extract the status
  OUT <- J$output
  KY <- matrix(NA, ncol = 2, nrow = length(OUT))
  colnames(KY) <- c("Status", "Karyotype")
  for (I in 1:length(OUT)) {
    x <- OUT[[I]]
    kary <- x$karyotype
    stat <- x$status
    KY[I,] <- c(stat, kary)
  }
  KY <- as.data.frame(KY, stringsAsFactors = FALSE)
  rownames(KY) <- padnames("RN", 1:length(OUT))
  KY$Status  <- factor(KY$Status)

  ## extract the (sub)clone ids
  clone <- unlist(lapply(OUT, function(x){
    if((x$status) %in% c("Success", "Fixable grammar error and success")) {
      1:length(lapply(x$parsing_result, "[[", "loss_gain_fusion_computing"))
    }
  }))
  ## extract the binary-mapped LGF data
  lgf <- lapply(OUT, function(x) { # Get the loss data from parsing_result
    if((x$status) %in% c("Success", "Fixable grammar error and success")) {
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
  list(Status = KY, LGF = df.lgf)
}

if(FALSE) { # maybe later...
setClass("LGF",
         slots = c(
           raw = "list",             # of data frame, rows = clones, columns = LGF-cytobands
           frequency = "data.frame", # rows = cytobands, columns = Loss, Gain, Fusion
           source = "character",     # file names
           size = "numeric",         # number of clones in each file
           CL = "data.frame"         # cytoband locations
         ))
}

Idioformat <- function(df, CL){
  ## element name will be the same as karyotype
  ## first example loss.1, loss.2
  Loss <- df[, grepl("Loss", names(df))]
  tags <- sub("Loss_", "", colnames(Loss))
  Gain <- df[, grepl("Gain", names(df))]
  Fusion <- df[, grepl("Fusion", names(df))]
  temp <- data.frame(Loss = colMeans(Loss > 0),
                     Gain = colMeans(Gain > 0),
                     Fusion = colMeans(Fusion > 0))
  rownames(temp) <- tags
  list(Frequency = temp[rownames(CL),], N = nrow(df))
}

readLGF <- function(files = NULL, folder = NULL, verbose = TRUE) {
  ## Figure out which files we want to read
  if (is.null(folder)) {
    folder  <-  "."
  }
  if (is.null(files)) {
    files  <- list.files(folder, pattern = "*.json",
                         full.names = FALSE) # character vector , one file per entry
  }
  if (length(files) < 1) {
    stop("No JSON input files to read.")
  }
  if (verbose) {
    message("Reading ", length(files), " file(s) from '", folder, "'.\n")
  }
  ## make sure that, when done,  we leave the working directory in the same state that we found it.
  home <- getwd()
  on.exit(setwd(home))
  setwd(folder)

  ## It might get large -- not said by Davis Guggenheim
  myJSON <- lapply(files, function(x) fromJSON(file = x)) # a list, one element per JSON file
  raw  <- lapply(myJSON, extractOneLGF)
  names(raw) <- sub(".json", "", basename(files))
  ## Get the summary statistics
  ick <- lapply(raw, function(R) {
    Idioformat(R$LGF, CL = cytobandLocations)
  })
  bundle <- do.call(cbind, lapply(ick, function(F) F$Frequency))
  size <- sapply(ick, function(F) F$N)
  list(source = files, raw = raw, frequency = bundle, size = size, CL = cytobandLocations)
}


