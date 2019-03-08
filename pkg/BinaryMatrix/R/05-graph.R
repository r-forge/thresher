### Based on the code from SPADE by Peng Qiu, which was implemented
### specifically for masa cyteometry data stored in FCS files.
downsample <- function(target, distanceMat, cutoff) {
  localDensity <- function(distanceMat, epsilon) {
    apply(distanceMat, 1, function(x) mean(x < epsilon))
  }
  ldens <- localDensity(distanceMat, cutoff)
  scalefactor <- sum(1/ldens)/length(ldens)
  P <- rbinom(length(ldens), 1, (1/ldens)/length(ldens) * target/scalefactor)
  P == 1
}
