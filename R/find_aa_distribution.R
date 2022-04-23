#'Find the distribution of amino acids in the polypeptide.
#'
#'Create a histogram of the distribution
#' @param polypep
#' @importFrom magrittr %>%
#' @return polypptide_plot
#' @author Miss Oriade Latifah Simpson, \email{s172084@@dtu.dk}
#' @export
#'
#' @examples
#' find_aa_distribution("HEGGHNRSRSSRTTSLKEET")
find_aa_distribution <- function(polypep){
  # split up the polypeptide and create distinct amino acids.
  uniq_amino <- polypep %>%
    stringr::str_split(pattern = stringr::boundary("character"), simplify = TRUE) %>%
    as.character() %>%
    unique()
  # count how many amino acids are in the polypeptide.
  counts <- sapply(uniq_amino, function(amino_acid) stringr::str_count(string = polypep, pattern =  amino_acid)) %>%
    as.data.frame()

  colnames(counts) <- c("counts")
  counts[["polypep"]] <- rownames(counts)
  # create a plot of amino acid counts.
  polypeptide_plot <- counts %>%
    ggplot2::ggplot(ggplot2::aes(x = polypep, y = counts, fill = polypep)) +
    ggplot2::geom_col(alpha = 0.6) +
    ggplot2::labs(x = "The Amino Acid",
                  y = "The counts",
                  subtitle = "The distribution of amino acids",
                  title = "Histogram",
                  caption = "Figure ") +
    ggthemes::theme_clean()+
    ggplot2::theme(legend.position = "none")

  return(polypeptide_plot)
}
