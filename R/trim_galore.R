#' Read TrimGalore stats file
#'
#' @param file stats file
#'
#' @return
#' @export
#' @md
read_trim_galore_stats <- function(file) {
  tibble(text = read_file(file)) %>%
    mutate(Raw.Reads = str_remove_all(str_match(., "reads processed:\\s+([\\d,]+)\\n")[2], ",")) %>%
    select(-text) %>%
    mutate(across(everything(), as.numeric))
}