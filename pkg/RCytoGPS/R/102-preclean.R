preclean <- function(x, targetColumns, dirt) {
  if (is.numeric(targetColumns)) {
    if (any(targetColumns < 1)) {
      stop("Only positive column numbers are allowed.\n")
    }
    if (any(targetColumns > dim(x)[2])) {
      nx <- deparse(substitute(x))
      stop("Column numbers exceed dimension of x (", nx, ").\n", sep="")
    }
    targetColumns <- colnames(x)[targetColumns]
  }
  if (!all(targetColumns %in% colnames(x))) {
    warning("Will skip unknown column names.\n")
  }
  for (target in targetColumns) {
    if (!(target %in% colnames(x))) {
      cat("\tSkipping '", target, "'\n", sep="", file=stderr())
    }
    vec <- x[,target]
    for (smudge in dirt) {
      vec <- gsub(smudge, "", vec, fixed = TRUE)
    }
    x[,target] <- vec
  }
  x
}
