#' Title Creation of a sequence of codons
#'
#' Creates different codons based on an RNA sequence.
#' @param any_rna_sequence
#' @param start indicating RNA strand as input to the function.
#'
#' @return Return character strings representing the codons
#' @author Oriade Latifah Simpson, \email{s172084@@dtu.dk}
#' @export
#'
#' @examples
#' my_codon <- make_codon(my_rna)
#' my_codon
make_codon <- function(any_rna_sequence, start = 1){
  # count the number of nucleotides
  number_of_nucleotides <- nchar(any_rna_sequence)

  # create sub strings of 3
  codons <- substring(any_rna_sequence,
                      first = seq(from = start, to = number_of_nucleotides-3+1, by = 3),
                      last = seq(from = 3+start-1, to = number_of_nucleotides, by = 3))
  return(codons)

}
