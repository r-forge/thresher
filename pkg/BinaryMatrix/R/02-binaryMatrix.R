setOldClass("dist")

setClass("BinaryMatrix",
         slots = c(
           binmat = "matrix",
           features = "data.frame",
           info = "list",
           history = "character")
)

BinaryMatrix <- function(binmat, features) {
  new("BinaryMatrix",
      binmat = binmat,
      features = features,
      info = list(),
      history = "Newly created.")
}

validBinaryMatrix <- function(object) {
  (ncol(object@binmat) == nrow(object@features)                # same size and
    & all(colnames(object@binmat) == rownames(object@features))) # same names
}
setValidity("BinaryMatrix", validBinaryMatrix)

setMethod("[", signature = "BinaryMatrix", function(x, i, j, ..., drop=FALSE) {
  new("BinaryMatrix",
      binmat = x@binmat[i, j, drop = drop],
      features = x@features[j, , drop = drop],
      info = x@info,
      history = c(x@history, "Subsetted."))
})

setMethod("dim", signature = "BinaryMatrix", function(x) {
  dim(x@binmat)
})

setMethod("summary", signature(object="BinaryMatrix"),
          function(object, ...) {
  cat("An object of the 'BinaryMatrix' class, of size\n")
  print(dim(object))
  cat("History:\n")
  print(object@history)
})

removeDuplicateFeatures <- function(object) {
  I <- object@info
  neverHit <- apply(object@binmat, 2, sum) == 0 
  I$notUsed <- object@features[neverHit,]
  binstring <- apply(object@binmat, 2, paste, collapse="") #vectors back to 0-1-strings
  dupstring <- duplicated(binstring)
  I$redundant <- object@features[dupstring,]
  retainBand <- !(dupstring | neverHit)
  binmat <- object@binmat[, retainBand, drop = FALSE]
  features <- object@features[retainBand,, drop = FALSE]
  history <- c(object@history, "Duplicate features removed.")
  new("BinaryMatrix",
      binmat = binmat,
      features = features,
      info = I,
      history = history)
}

