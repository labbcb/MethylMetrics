#' Calculate coverage of a given matrix
#' 
#' This function expects integer matrix where columns are samples and rows are
#' bases. It normally inputs [bsseq::getCoverage()] result object.
#'
#' @param cov integer matrix
#' @param min_cov minimum value to count base as covered
#'
#' @return data frame with the following columns:
#' * `sample` sample names
#' * `min_cov` minimum coverage threshold (same as `min_cov` parameter)
#' * `covered` number of covered bases
#' * `total` number of bases
#' 
#' @export
#' @md
#'
#' @examples
#' library(bsseqData)
#' data(BS.cancer.ex)
#' getCoverage(BS.cancer.ex) %>%
#'   calc_coverage(min_cov = 10)
calc_coverage <- function(cov, min_cov = 1) {
  total <- nrow(cov)
  covered <- apply((cov), 2, function(x) sum(x >= min_cov))
  
  tibble::tibble(sample=names(covered), min_cov, covered, total) %>%
    dplyr::mutate(freq = covered / total, perc = freq * 100)
}

#' Calculate coverage of a BSseq object where bases overlap given regions
#' 
#' This function inputs a [bsseq::BSseq] object and a set of regions
#' (data frame or [GRanges]). Bases that overlaps one of these regions
#' are considered for coverage calculation.
#' 
#' Regions as data frame must have the columns: `seqnames`, `start` and `end`.
#'
#' @param regions data frame or GenomicRanges
#' @param type region type (exon, intron, gene symbol, etc)
#' @param bs BSseq object
#' @param min_cov minimum value to count base as covered
#'
#' @return data frame with the following columns:
#' * `sample` sample names
#' * `type` region type
#' * `min_cov` minimum coverage threshold (same as `min_cov` parameter)
#' * `covered` number of covered bases
#' * `total` number of bases
#' 
#' @export
#' @md
#'
#' @examples
#' library(bsseqData)
#' data(BS.cancer.ex)
#' calc_coverage_regions(
#'   regions = GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 2*10^7)),
#'   type = "custom",
#'   bs = BS.cancer.ex,
#'   min_cov = 10)
calc_coverage_regions <- function(regions, type, bs, min_cov = 1) {
  if (!is(regions, "GenomicRanges")) {
    regions <- makeGRangesFromDataFrame(regions)
  }
  
  subsetByOverlaps(bs, regions) %>%
    getCoverage() %>%
    calc_coverage(min_cov = min_cov) %>%
    tibble::add_column(type, .before = "min_cov")
}

#' Calculate coverage of a BSseq object for a given chromosome
#' 
#' This function calculates coverage of bases of a single chromosome
#' (reference name, `seqnames`).
#'
#' @param chr single chromosome
#' @param bs BSseq object
#' @param min_cov minimum value to count base as covered
#'
#' @return data frame with the following columns:
#' * `sample` sample names
#' * `chr` chromosome
#' * `min_cov` minimum coverage threshold (same as `min_cov` parameter)
#' * `covered` number of covered bases
#' * `total` number of bases
#' 
#' @export
#' @md
#'
#' @examples
#' library(bsseqData)
#' data(BS.cancer.ex)
#' calc_coverage_chr("chr21", bs = BS.cancer.ex, min_cov = 10)
calc_coverage_chr <- function(chr, bs, min_cov = 1) {
  getCoverage(chrSelectBSseq(bs, chr)) %>%
    calc_coverage(min_cov = min_cov) %>%
    tibble::add_column(chr, .before = "min_cov")
}

#' Calculate coverage metrics of a given BSseq object
#' 
#' This is the main function to calculate coverage metrics.
#' It inputs a [bsseq::BSseq] object, zero or more chromosomes
#' (reference names, `seqnames`), and zero or more regions.
#' 
#' The `regions_list` must be a names list of data frames and 
#' [GenomicRanges::GRanges] objects. Names of list elements are used to
#' differentiate metrics from one set of regions to others (column `type`).
#' 
#' If `chrs` or `regions_list` is not provided it will not calculate for these
#' genomic intervals.
#'
#' @param bs BSseq object
#' @param chrs zero or more chromosomes
#' @param regions_list named list with zero or more regions
#' @param min_cov 
#'
#' @return named list of the following data frames:
#' * `cov` genome-wide coverage metrics (see [calc_coverage()])
#' * `cov_chr` metrics summarized by chromosomes (see [calc_coverage_chr()])
#' * `cov_regions` metrics summarized by regions (see [calc_coverage_regions()])
#' 
#' @export
#' @md
#'
#' @examples
#' library(bsseqData)
#' data(BS.cancer.ex)
#' calc_coverage_metrics(BS.cancer.ex,
#'   chrs = c("chr21","chr22"),
#'   regions_list = list(custom = GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 2*10^7))),
#'   min_cov = 10)
calc_coverage_metrics <- function(bs, chrs, regions_list, min_cov = 1) {
  cov <- calc_coverage(getCoverage(bs), min_cov=min_cov)
  cov_chr <- purrr::map_df(chrs, calc_coverage_chr, bs=bs, min_cov=min_cov)
  cov_regions <- purrr::map2_df(regions_list, names(regions_list),
    calc_coverage_regions, bs = bs, min_cov = min_cov)
  
  list(cov=cov, cov_chr=cov_chr, cov_regions=cov_regions)
}

#' Write coverage metrics data as tabular files
#' 
#' This function inputs the result object of [calc_coverage_metrics()].
#' Each data frame in the named list is written as tabular file (with header).
#' The base name of each file is defined by name of element plus `ext`.
#' Some of files may not create if there are missing metrics results,
#' for instance no regions were passed to `calc_coverage_metrics` function.
#' 
#' The `prefix_path` are added before base name of each file.
#' It can be the cytosine context for example (CpG, CHG, CHH).
#' It also can be use to change destination directory (example: `metrics/`).
#' 
#'
#' @param cov_metrics named list of metrics data (see [calc_coverage_metrics()])
#' @param prefix_path single string to add to base file name as prefix
#' @param delim delimiter used to separate values. Default is ",".
#' @param ext file extension. Default is ".csv".
#'
#' @return file paths to created files.
#' Possible files are (prefix and ext omitted):
#' * `cov` genome-wide coverage metrics (see [calc_coverage()])
#' * `cov_chr` metrics summarized by chromosomes (see [calc_coverage_chr()])
#' * `cov_regions` metrics summarized by regions (see [calc_coverage_regions()])
#' 
#' @export
#' @md
#'
#' @examples
#' \dontrun{
#' library(bsseqData)
#' data(BS.cancer.ex)
#' calc_coverage_metrics(BS.cancer.ex,
#'   chrs = c("chr21","chr22"),
#'   regions_list = list(custom = GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 2*10^7))),
#'   min_cov = 10) %>%
#'   write_coverage_metrics(prefix = "metrics/CpG_")
#'}
write_coverage_metrics <- function(cov_metrics, prefix_path = "", delim = ",", ext = ".csv") {
  for (name in names(cov_metrics)) {
    readr::write_csv(cov_metrics[[name]], paste0(prefix_path, name, ext))
  }
}

