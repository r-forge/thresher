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
          function(height, what, col = "blue", altcol = "#FE4C4",
                   ylab = "Percent", h = NULL, ...) {
  if (is.list(what)) what <- unlist(what)
  if ( length(what) != 1) stop("'what' must identify exactly one data column.")
  genomeBarplot(height@DATA[ , what], col = col, altcol = altcol, ylab = ylab, h=h)
})

setMethod("image", "CytobandData", function(x, chr, what,
           pal = palette(), horiz = FALSE, nrows = 2, labels = NULL) {
  if (length(chr) != 1) {
    stop("Invalid chromosome value.")
  }
  if (is.list(what)) {
    if (length(what) != 2) {
      stop("'what' must be a list o length exactly two.")
    }
    if(chr == "all") {
      biIdiogram(x@DATA, what[[1]], what[[2]],
                 pal = pal, horiz = horiz, nrows = nrows)
    } else {
      plot2Chrom(x@DATA, what[[1]], what[[2]],
                 chr = chr, pal = pal,
                 horiz = !horiz) # changed our minds about what "horiz" means
    }
  } else { # now we must have chr equal to a character vector
    if (chr == "all") {
      if (length(what) == 1) {
        makeIdiogram(x@DATA, what, color = pal)
      } else {
        stackIdiogram(x@DATA, what, pal = pal,
                      horiz = !horiz, # changed our minds about "horiz"
                      nrows = nrows)
      }
    } else {
      plot1Chrom(x@DATA, what, chr = chr,
                 labels = labels, pal = pal,
                 horiz = !horiz) # changed our minds about what "horiz" means
    }
  }
  invisible(x)
})
