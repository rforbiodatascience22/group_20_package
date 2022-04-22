#' Title Creation of a sequence of DNA nucleotides.
#'
#'Creates different lengths of nucleotides of a DNA sequence.
#' @param Set sequence_length indicating the length of the DNA strand.
#'
#' @return Return character strings representing the DNA.
#' @author Oriade Latifah Simpson, \email{s172084@@dtu.dk}
#' @export
#'
#' @examples
#' my_dna <- make_dna(5)
#' my_dna
make_dna <- function(sequence_length){
  # sample nucleotides with replacement to create a sequence of nucleotides of a specified length.
  dna_sample <- sample(c("A", "T", "G", "C"), size = sequence_length, replace = TRUE)

  # collapse the string together.
  the_dna_sequence <- paste0(dna_sample, collapse = "")
  return(the_dna_sequence)
}
