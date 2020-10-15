setOldClass("dist")

setClass("BinaryMatrix",
         slots = c(
           binmat = "matrix",
           columnInfo = "data.frame",
           rowInfo = "data.frame",
           info = "list",
           history = "character")
)

BinaryMatrix <- function(binmat, columnInfo, rowInfo) {
  if (missing(rowInfo)) {
    if (is.null(rownames(binmat))) {
      rownames(binmat) <- paste("R", 1:nrow(binmat), sep = "")
    }
    rowInfo <- data.frame(Names = rownames(binmat))
    rownames(rowInfo) <- rownames(binmat)
  }
  if (missing(columnInfo)) {
    if (is.null(colnames(binmat))) {
      colnames(binmat) <- paste("X", 1:ncol(binmat), sep = "")
    }
    columnInfo <- data.frame(Names = colnames(binmat))
    rownames(columnInfo) <- colnames(binmat)
  }
  new("BinaryMatrix",
      binmat = binmat,
      columnInfo = columnInfo,
      rowInfo = rowInfo,
      info = list(),
      history = "Newly created.")
}

validBinaryMatrix <- function(object) {
  (ncol(object@binmat) == nrow(object@columnInfo)                # same size and
    & all(colnames(object@binmat) == rownames(object@columnInfo)) # same names
    & (nrow(object@binmat) == nrow(object@rowInfo))                # same size and
    & all(rownames(object@binmat) == rownames(object@rowInfo))) # same names
}
setValidity("BinaryMatrix", validBinaryMatrix)

setMethod("[", signature = "BinaryMatrix", function(x, i, j, ..., drop=FALSE) {
  new("BinaryMatrix",
      binmat = x@binmat[i, j, drop = drop],
      columnInfo = x@columnInfo[j, , drop = drop],
      rowInfo  = x@rowInfo[i, , drop = drop],
      info = x@info,
      history = c(x@history, "Subsetted."))
})

setMethod("print", signature = "BinaryMatrix", function(x, ...) {
  show(x)
})

setMethod("show", signature = "BinaryMatrix", function(object) {
  M <- min(R <- nrow(object), 10)
  N <- min(C <- ncol(object), 10)
  B <- object@binmat[1:M, 1:N]
  if (M < R | N < C) {
    cat("The full matrix contains", R, "rows and", C,
        "columns.\nTo see the whole thing, use 'x@binmat'",
        "to extract the underlying matrix.\n Here is it the upper-left corner:\n")
  }
  print(B)
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

setMethod("t",  signature = "BinaryMatrix",
          function(x) {
  new("BinaryMatrix",
      binmat = t(x@binmat),
      columnInfo = x@rowInfo,
      rowInfo = x@columnInfo,
      info = x@info,
      history = c(x@history, "transposed")
      )
})

removeDuplicateFeatures <- function(object) {
  I <- object@info
  neverHit <- apply(object@binmat, 2, sum) == 0 
  I$notUsed <- object@columnInfo[neverHit,]
  binstring <- apply(object@binmat, 2, paste, collapse="") #vectors back to 0-1-strings
  dupstring <- duplicated(binstring)
  I$redundant <- object@columnInfo[dupstring & !neverHit,]
  retainBand <- !(dupstring | neverHit)
  binmat <- object@binmat[, retainBand, drop = FALSE]
  columnInfo <- object@columnInfo[retainBand,, drop = FALSE]
  history <- c(object@history, "Duplicate features removed.")
  new("BinaryMatrix",
      binmat = binmat,
      columnInfo = columnInfo,
      rowInfo = object@rowInfo,
      info = I,
      history = history)
}

