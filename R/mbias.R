calc_mbias <- function(bam) {
  tags <- c("XR", "XG", "XM")
  param <- ScanBamParam(tag = tags)
  dat <- scanBam(bam, param = param)
}