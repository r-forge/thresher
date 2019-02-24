setClass("ThreshedBinaryMatrix",
         slots = c(
           thresher = "Thresher",
           reaper = "Reaper"
           ),
         contains = "BinaryMatrix"
)

threshLGF <- function(object, cutoff = 0) {
  I <- object@info
  I$cutoff <- cutoff
  thresh <- Thresher(object@binmat, method="auer.gervini")
  reap   <- Reaper(thresh, useLoadings=TRUE,
                     cutoff = cutoff, metric="pearson")
  keep <- thresh@delta  > cutoff
  new("ThreshedBinaryMatrix",
      binmat = object@binmat[, keep],
      features = object@features[keep,],
      info = I,
      history = c(object@history, "Threshed."),
      thresher = thresh,
      reaper = reap)
}


