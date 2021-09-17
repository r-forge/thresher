setClass("CytobandData",
         slots = c(DATA = "data.frame",
                   INFO = "data.frame"
                   ))

setValidity("CytobandData", function(object) {
  msg <- ""
  if (nrow(object@INFO) != ncol(object@DATA)) {
    msg <- "Mismatched sizes."
  } else if (!all(rownames(object@INFO) == colnames(object@DATA))) {
    msg <- "Mismatched names."
  }
  if (!all(c("Chromosome", "loc.start", "loc.end") %in% colnames(object@DATA))) {
    msg <- paste(msg, "Missing required data column.")
  }
  if (!all(c("Label", "Description") %in% colnames(object@INFO))) {
    msg <- paste(msg, "Missing required info column.")
  }
  ifelse(msg == "", TRUE, msg)
})

CytobandData <- function(data, info, genome = NULL) {
  if (missing(info)) {
    info <- data.frame(Label = colnames(data),
                       Description = colnames(data))
    rownames(info) <- colnames(data)
  }
  if (!is.null(genome)) {
    ginfo <- data.frame(Label = colnames(genome),
                        Description = colnames(genome))
    rownames(ginfo) <- colnames(genome)
    data <- cbind(data, genome)
    info <- rbind(info, genome)
  }
  new("CytobandData", DATA = data, INFO = info)
}

setMethod("summary", "CytobandData", function(object, ...) {
  cat("A object of the 'CytobandData' class.\n")
  summary(object@DATA)
})

setMethod("barplot", "CytobandData",
          function(height, what, col = "blue", altcol = "#FED4C4",
                   ylab = "Percent", hline = NULL,
                   xform = function(x) x, ...) {
  if (is.list(what)) what <- unlist(what)
  if ( length(what) != 1) stop("'what' must identify exactly one data column.")
  genomeBarplot(xform(height@DATA[ , what]),
                col = col, altcol = altcol, ylab = ylab, h=hline)
})

setMethod("image", "CytobandData", function(x, chr, what,
         pal = palette(), nrows = 2, labels = NULL,
         horiz = FALSE, axes = chr != "all", debug = FALSE, legend = FALSE,
         sigcolumn = NA, sigcut = 0.01, alpha = 63, clip = FALSE) {
  if (length(chr) != 1) {
    stop("Invalid chromosome value.")
  }
  if (is.list(what)) {
    if (length(what) != 2) {
      stop("'what' must be a list of length exactly two.")
    }
    if(chr == "all") {
      if(debug) cat("biIdiogram\n", file = stderr())
      biIdiogram(x@DATA, what[[1]], what[[2]],
                 pal = pal, nrows = nrows, horiz = horiz,
                 axes = axes, legend = legend,
                 sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
    } else {
      if(debug) cat("plot2Chrom\n", file = stderr())
      plot2Chrom(x@DATA, what[[1]], what[[2]],
                 chr = chr, pal = pal,
                 horiz = !horiz, # changed our minds about what "horiz" means
                 axes = axes, legend = legend,
                 sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha,
                 clip = clip)
    }
  } else { # now we must have chr equal to a character vector
    if (chr == "all") {
      if (length(what) == 1) {
        if(debug) cat("makeIdiogram\n", file = stderr())
        makeIdiogram(x@DATA, what, color = pal, axes = axes, legend = legend,
                     sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
      } else {
        if(debug) cat("stackIdiogram\n", file = stderr())
        stackIdiogram(x@DATA, what, pal = pal,
                      horiz = !horiz, # changed our minds about "horiz"
                      axes = axes, nrows = nrows, legend = legend,
                      sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
      }
    } else {
      if(debug) cat("plot1Chrom\n", file = stderr())
      plot1Chrom(x@DATA, what, chr = chr,
                 labels = labels, pal = pal, axes = axes,
                 horiz = !horiz, # changed our minds about what "horiz" means
                 legend = legend,
                 sigcolumn = sigcolumn, sigcut = sigcut, alpha = alpha)
    }
  }
  invisible(x)
})
