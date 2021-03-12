#' Read methylcx M-Bias file
#' 
#' @param file path M-Bias file
#'
#' @return data frame with the following columns:
#' * Context: CpG, CHG or CHH
#' * Strand: Forward or Reverse
#' * Cycle: Sequencing cycle
#' * Count.Methylated:
#' * Count.Unmethylated:
#' * Percent.Methylation:
#' * Coverage: 
#' 
#' @export
#' @md
read_methycx_mbias <- function(file) {
  readr::read_csv(file) %>%
    rename(Count.Methylated = Methylated, Count.Unmethylated = Unmethylated) %>%
    mutate(Percent.Methylation = Count.Methylated / Coverage * 100)
}